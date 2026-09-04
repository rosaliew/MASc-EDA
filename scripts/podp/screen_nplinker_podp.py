#!/usr/bin/env python3
"""Cheap pre-flight screen for NPLinker's PoDP mode.

NPLinker's ``DatasetArranger`` in ``podp`` mode fetches the PoDP record, reads
the GNPS task id from ``metabolomics.project.molecular_network``, downloads that
task's results, and derives ``gnps/file_mappings.tsv`` -- the spectrum -> source
file join -- then folds it into ``strain_mappings.json``. That is exactly the
bridge the 48 molecular-family links are missing.

It only works if three things hold, and all three are checkable without
downloading anything:

  1. the PoDP API still serves the record  (HEAD on the project endpoint)
  2. the record carries a GNPS task id      (already in ground_truth.csv)
  3. the GNPS task is still alive and is a workflow NPLinker parses
     (GET status.jsp, read the Workflow and Status rows -- the same probe
     ``gnps_format_from_task_id`` uses)

Nothing here downloads a task archive. Writes reports/nplinker_podp_screen.csv.
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
import sys
import time
from pathlib import Path

import httpx
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[2]  # repo root (this file is scripts/podp/*)
PODP = ROOT / "results" / "podp"   # was data/PoDP, moved in the 2026-09 reorg
GT = PODP / "podp_bgc_ms2_links.csv"
OUT = PODP / "reports" / "nplinker_podp_screen.csv"

PODP_PROJECT_URL = "https://pairedomicsdata.bioinformatics.nl/api/projects/{}"
GNPS1_STATUS_URL = "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}"
GNPS2_STATUS_URL = "https://gnps2.org/status?task={}"

# workflows NPLinker's GNPSFormat enum knows how to parse
SUPPORTED = {
    "METABOLOMICS-SNETS",
    "METABOLOMICS-SNETS-V2",
    "FEATURE-BASED-MOLECULAR-NETWORKING",
}
TASK_RE = re.compile(r"task=([0-9a-f]{32})", re.I)

FIELDS = [
    "study_id", "podp_id", "podp_id_live", "gt_links", "links_needing_bridge",
    "task_id", "task_source", "podp_api", "gnps_status", "gnps_workflow",
    "nplinker_supported", "verdict", "note",
]


def load_targets(only_unresolved: bool) -> dict[str, dict]:
    """Group ground-truth rows by study, carrying every task id seen."""
    with GT.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    studies: dict[str, dict] = {}
    for r in rows:
        study = r.get("study_id") or r["podp_id"]
        s = studies.setdefault(study, {
            "study_id": study, "podp_id": r["podp_id"], "gt_links": 0,
            "links_needing_bridge": 0, "project_tasks": set(), "network_tasks": set(),
        })
        s["gt_links"] += 1
        if r.get("genome_unresolved_reason") == "gnps_resolution_not_run":
            s["links_needing_bridge"] += 1
        if r.get("gnps_task_id"):
            s["project_tasks"].add(r["gnps_task_id"].strip())
        m = TASK_RE.search(r.get("network_nodes_url") or "")
        if m:
            s["network_tasks"].add(m.group(1).lower())

    if only_unresolved:
        studies = {k: v for k, v in studies.items() if v["links_needing_bridge"]}
    return studies


def probe_podp(client: httpx.Client, podp_id: str) -> tuple[str, str]:
    """Is the PoDP record still served, and under which version?

    Returns (status, live_id). The API needs the full versioned id -- stripping
    the .N suffix 404s, which is why NPLinker's config.podp_id carries it. If
    the exact version is gone the record has been revised upstream, so walk
    neighbouring versions to find what PoDP serves now.
    """
    def get(pid: str) -> str:
        try:
            r = client.get(PODP_PROJECT_URL.format(pid), timeout=30)
        except Exception as exc:
            return f"error:{type(exc).__name__}"
        if r.status_code == 200:
            return "ok" if r.text.lstrip().startswith("{") else "non_json"
        return f"http_{r.status_code}"

    status = get(podp_id)
    if status == "ok" or "." not in podp_id:
        return status, podp_id if status == "ok" else ""

    base, _, ver = podp_id.rpartition(".")
    if not ver.isdigit():
        return status, ""
    for cand in range(1, 10):
        if cand == int(ver):
            continue
        time.sleep(0.2)
        if get(f"{base}.{cand}") == "ok":
            return "stale_version", f"{base}.{cand}"
    return status, ""


def probe_gnps(client: httpx.Client, task_id: str) -> tuple[str, str]:
    """Return (status, workflow) for a GNPS1 task, without downloading it."""
    try:
        r = client.get(GNPS1_STATUS_URL.format(task_id), timeout=30)
    except Exception as exc:
        return f"error:{type(exc).__name__}", ""
    if r.status_code != 200:
        return f"http_{r.status_code}", ""

    soup = BeautifulSoup(r.text, features="html.parser")

    def cell(label: str) -> str:
        # match on the rendered text: several of these th tags wrap their
        # label in nested markup, which a string= match silently misses
        for th in soup.find_all("th"):
            if th.get_text(strip=True) != label:
                continue
            td = th.find_next_sibling("td")
            if not td:
                return ""
            # the value can be nested (Status wraps DONE in <div><a>), so read
            # rendered text rather than contents[0], then drop the trailing
            # [Clone] / [View All Library Hits] action links
            return td.get_text(" ", strip=True).split("[")[0].strip()
        return ""

    workflow, status = cell("Workflow"), cell("Status")
    workflow = re.split(r"\s*\(version", workflow)[0].strip()
    if not workflow and not status:
        # no task table rendered at all -- the id is unknown to the server
        return "not_found", ""
    return (status or "unknown"), workflow


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all-studies", action="store_true",
                    help="screen every study, not just those needing the bridge")
    ap.add_argument("--delay", type=float, default=0.4, help="seconds between requests")
    args = ap.parse_args()

    studies = load_targets(only_unresolved=not args.all_studies)
    if not studies:
        print("nothing to screen", file=sys.stderr)
        return 1
    print(f"screening {len(studies)} studies "
          f"({sum(s['links_needing_bridge'] for s in studies.values())} links need the bridge)\n")

    out: list[dict] = []
    podp_cache: dict[str, str] = {}
    task_cache: dict[str, tuple[str, str]] = {}

    with httpx.Client(follow_redirects=True,
                      headers={"User-Agent": "podp-screen/1.0"}) as client:
        for i, (study, s) in enumerate(sorted(studies.items()), 1):
            tasks = [(t, "molecular_network") for t in sorted(s["project_tasks"])]
            tasks += [(t, "network_nodes_url") for t in sorted(s["network_tasks"])
                      if t not in s["project_tasks"]]

            if s["podp_id"] not in podp_cache:
                podp_cache[s["podp_id"]] = probe_podp(client, s["podp_id"])
                time.sleep(args.delay)
            api, live_id = podp_cache[s["podp_id"]]

            if not tasks:
                out.append({**{f: "" for f in FIELDS},
                            "study_id": study, "podp_id": s["podp_id"],
                            "gt_links": s["gt_links"],
                            "links_needing_bridge": s["links_needing_bridge"],
                            "podp_id_live": live_id,
                            "podp_api": api, "verdict": "no_task_id",
                            "note": "PoDP record names no GNPS task; the arranger "
                                    "has nothing to download for this study."})
                print(f"[{i}/{len(studies)}] {study[:8]}  no task id")
                continue

            for task, src in tasks:
                if task not in task_cache:
                    task_cache[task] = probe_gnps(client, task)
                    time.sleep(args.delay)
                status, workflow = task_cache[task]

                supported = workflow in SUPPORTED
                if status == "not_found":
                    verdict, note = "task_gone", "GNPS does not know this task id."
                elif status.startswith(("http_", "error:")):
                    verdict, note = "unreachable", f"status probe returned {status}."
                elif status.upper() != "DONE":
                    verdict, note = "not_done", f"task status is {status!r}, not DONE."
                elif not supported:
                    verdict, note = "unsupported_workflow", (
                        f"workflow {workflow!r} is outside NPLinker's GNPSFormat enum.")
                else:
                    verdict, note = "ready", "task is DONE and parseable by NPLinker."
                if verdict == "ready" and api == "stale_version":
                    verdict, note = "ready_after_version_bump", (
                        f"GNPS task is fine, but PoDP no longer serves "
                        f"{s['podp_id']}; it now serves {live_id}. Point the "
                        "NPLinker config at the live id.")
                elif verdict == "ready" and api != "ok":
                    verdict, note = "podp_api_down", (
                        f"GNPS task is fine but the PoDP API returned {api}; "
                        "the arranger fetches the record first.")

                out.append({
                    "study_id": study, "podp_id": s["podp_id"],
                    "gt_links": s["gt_links"],
                    "links_needing_bridge": s["links_needing_bridge"],
                    "podp_id_live": live_id,
                    "task_id": task, "task_source": src, "podp_api": api,
                    "gnps_status": status, "gnps_workflow": workflow,
                    "nplinker_supported": "yes" if supported else "no",
                    "verdict": verdict, "note": note,
                })
                print(f"[{i}/{len(studies)}] {study[:8]}  {task[:10]}  "
                      f"{status:12} {workflow[:34]:36} {verdict}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)

    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(out)} rows)\n")
    print("verdicts:")
    for k, v in collections.Counter(r["verdict"] for r in out).most_common():
        print(f"  {v:4}  {k}")
    print("\nPoDP API:")
    for k, v in collections.Counter(r["podp_api"] for r in out).most_common():
        print(f"  {v:4}  {k}")

    stale = {r["podp_id"]: r["podp_id_live"] for r in out
             if r["podp_api"] == "stale_version"}
    if stale:
        print("\nlocal descriptors PoDP has since revised:")
        for old, new in sorted(stale.items()):
            print(f"  {old}  ->  {new}")

    ok_verdicts = {"ready", "ready_after_version_bump"}
    ready = {r["study_id"] for r in out if r["verdict"] in ok_verdicts}
    recoverable = sum(s["links_needing_bridge"] for k, s in studies.items() if k in ready)
    total = sum(s["links_needing_bridge"] for s in studies.values())
    print(f"\nstudies with at least one ready task: {len(ready)}/{len(studies)}")
    print(f"bridge-needing links in those studies : {recoverable}/{total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
