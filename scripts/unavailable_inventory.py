#!/usr/bin/env python3
"""Inventory every ground-truth data file we cannot currently obtain, and why.

Two kinds of gap get conflated easily, so they are coded separately here:

  * ``mapping_gap``   -- the file exists on disk, but nothing joins it to the
                         link. Fixable in code; no download required.
  * ``file_gap``      -- the file is genuinely not in hand. Each carries a
                         reason and a recoverability verdict.

Recoverability is the actionable column:

  ``automatable``     -- a known route exists and is not yet implemented
  ``manual``          -- a human has to fetch or adjudicate it
  ``deferred``        -- reachable, but deliberately not being pursued
  ``unavailable``     -- no known route; the source does not hold it

Reads ground_truth.csv, refetch_plan_all.csv and mibig_migration_report.csv.
Writes data/PoDP/reports/unavailable_data_inventory.csv plus a markdown digest.
"""

from __future__ import annotations

import csv
import collections
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PODP = ROOT / "data" / "PoDP"
REPORTS = PODP / "reports"
GT = PODP / "ground_truth.csv"
PAIRINGS = PODP / "strain_pairings.csv"
PLAN = REPORTS / "refetch_plan_all.csv"
MIGRATION = REPORTS / "mibig_migration_report.csv"
JGI = REPORTS / "jgi_worklist_297c364c.csv"
OUT_CSV = REPORTS / "unavailable_data_inventory.csv"
OUT_MD = REPORTS / "unavailable_data_inventory.md"

FIELDS = [
    "gap_id", "scope", "file_kind", "gap_type", "reason_code", "recoverability",
    "n_ground_truth_rows", "n_strain_pairing_rows", "n_distinct_items", "items",
    "studies", "detail", "next_action",
]


def rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def study_of(row: dict) -> str:
    return (row.get("study_id") or row.get("podp_id") or "")[:8]


def _cap(items, n=12) -> str:
    """Pipe-join, truncating long lists so the cell stays readable."""
    items = sorted(set(items))
    if len(items) <= n:
        return "|".join(items)
    return "|".join(items[:n]) + f"|... (+{len(items) - n} more)"


def collect() -> list[dict]:
    gt, plan, mig, jgi = rows(GT), rows(PLAN), rows(MIGRATION), rows(JGI)
    sp_by_acc: dict[str, int] = collections.Counter(
        r["genome_accession"] for r in rows(PAIRINGS) if r["genome_accession"])
    out: list[dict] = []

    def add(**kw):
        kw.setdefault("n_strain_pairing_rows", "")
        kw.setdefault("n_distinct_items", "")
        kw.setdefault("items", "")
        kw.setdefault("studies", "")
        kw["gap_id"] = f"G{len(out) + 1:02d}"
        out.append({f: kw.get(f, "") for f in FIELDS})

    # ---- 1. MS2 spectra: never extracted -------------------------------
    pending = [r for r in gt if not r["file_ms2_spectrum"]]
    add(scope="ground_truth", file_kind="ms2_spectrum", gap_type="not_attempted",
        reason_code="spectrum_extraction_not_run", recoverability="automatable",
        n_ground_truth_rows=len(pending),
        studies=_cap(study_of(r) for r in pending),
        detail="No per-link MGF has been cut yet. 66 rows carry an explicit "
               "ms2_scan and 49 carry a GNPS component; both are extractable "
               "from files already downloaded.",
        next_action="Extract scans from the mzML/mzXML already on disk.")

    # ---- 2. Molecular-family links: genome + run unmapped ---------------
    unres = [r for r in gt if r["genome_unresolved_reason"] == "gnps_resolution_not_run"]
    add(scope="ground_truth", file_kind="genome", gap_type="mapping_gap",
        reason_code="gnps_component_bridge_not_run", recoverability="automatable",
        n_ground_truth_rows=len(unres),
        studies=_cap(study_of(r) for r in unres),
        detail="Every one of these is a 'GNPS molecular family' link. They carry "
               "no MS2_URL, so the ms2_url_bridge cannot fire and no genome is "
               "attached. The genome FILES are on disk -- 297c364c alone holds "
               "130 -- this is a join that was never built, not missing data. "
               "SCREENED 2026-08-27: NPLinker's own podp mode supplies this join "
               "for free. All 14 affected studies pass pre-flight (48/48 links).",
        next_action="Run NPLinker in podp mode per study rather than hand-building "
                    "a bridge; its gnps/file_mappings.tsv is a spectrum x strain "
                    "matrix. Verified on task c22f44b14a: 27/27 real strain "
                    "columns resolve to a genome accession. See "
                    "nplinker_podp_screen.csv.")

    nomap = [r for r in gt if not r["file_ms2_run"] and not r["ms2_run_url"]]
    add(scope="ground_truth", file_kind="ms2_run", gap_type="mapping_gap",
        reason_code="no_ms2_url_on_family_link", recoverability="automatable",
        n_ground_truth_rows=len(nomap),
        studies=_cap(study_of(r) for r in nomap),
        detail="Same 48 family links. PoDP records a network component, not a "
               "run file, so there is no URL to resolve until the component is "
               "expanded to its member spectra. NPLinker's podp mode does that "
               "expansion as part of arrange_gnps().",
        next_action="Same as G02 -- NPLinker podp mode, not a hand-built bridge.")

    # ---- 3. A genuine MS2 download failure ------------------------------
    failed = [r for r in gt if not r["file_ms2_run"] and r["ms2_run_url"]]
    if failed:
        add(scope="ground_truth", file_kind="ms2_run", gap_type="file_gap",
            reason_code="massive_download_failed", recoverability="manual",
            n_ground_truth_rows=len(failed), n_distinct_items=len(failed),
            items=_cap(r["ms2_run_url"].rsplit("/", 1)[-1] for r in failed),
            studies=_cap(study_of(r) for r in failed),
            detail="URL is present and well-formed but the file was not "
                   "retrieved. Only MS2 run in this state.",
            next_action="Retry from the MassIVE FTP path by hand.")

    # ---- 4. BGC absent from the MIBiG release --------------------------
    fate = {r["mibig_id"]: r for r in mig}
    nobgc = [r for r in gt if not r["file_bgc"]]
    by_fate = collections.defaultdict(list)
    for r in nobgc:
        by_fate[fate.get(r["mibig_id"], {}).get("fate", "unknown")].append(r)
    for f, rs in sorted(by_fate.items()):
        recov = "unavailable" if f == "absent_from_both" else "manual"
        detail = ("Cited by PoDP but present in neither MIBiG 3.1 nor 4.0. The "
                  "accession may never have been issued."
                  if f == "absent_from_both" else
                  "Active with a structure in MIBiG 3.1, but dropped from the "
                  "4.0 bundle. The 3.1 record still exists upstream.")
        action = ("Confirm the accession against MIBiG directly; if invalid, "
                  "the link's BGC claim cannot be grounded."
                  if f == "absent_from_both" else
                  "Decide whether to pin these to their 3.1 record or accept "
                  "the loss. Note BGC0000962 (barbamide) is NPOmix-covered.")
        add(scope="ground_truth", file_kind="bgc", gap_type="file_gap",
            reason_code=f"mibig_{f}", recoverability=recov,
            n_ground_truth_rows=len(rs),
            n_distinct_items=len({r["mibig_id"] for r in rs}),
            items=_cap(r["mibig_id"] for r in rs),
            studies=_cap(study_of(r) for r in rs),
            detail=detail, next_action=action)

    # ---- 5. Genome file gaps, from the refetch plan ---------------------
    gt_by_acc = collections.defaultdict(int)
    for r in gt:
        if r["genome_accession"]:
            gt_by_acc[r["genome_accession"]] += 1

    buckets = {
        "skipped_wrong_strain": (
            "ncbi_returned_wrong_strain", "manual",
            "NCBI's assembly index ignores strain tokens, so a strain-name "
            "query silently returns an arbitrary genome of the same genus. "
            "verify_strain() caught these; a content hash found six accessions "
            "collapsing onto one file. Quarantined, never promoted.",
            "Re-resolve via BioSample [Strain] then elink to assembly -- "
            "confirmed working on CNS237, CNH643, CNS205."),
        "skipped_metagenome_scale": (
            "metagenome_scale_assembly", "manual",
            "Sequence is present and complete, but the assembly is "
            "metagenome-scale rather than an isolate. Held back by the scale "
            "gate so it cannot be mistaken for a single-strain genome.",
            "A judgement call, not a gap: antiSMASH will run on these. Promote "
            "deliberately if metagenome-derived BGCs are acceptable."),
        "skipped_not_sequence": (
            "no_sequence_after_refetch", "manual",
            "Re-fetch still yielded a record with no ORIGIN block -- typically "
            "a WGS/CON master record that indexes contigs without carrying them.",
            "Fetch the component contigs, or the assembly FTP _genomic.gbff.gz "
            "for the corresponding GCA/GCF."),
    }
    for verdict, (code, recov, detail, action) in buckets.items():
        rs = [r for r in plan if r["promoted"] == verdict]
        if not rs:
            continue
        add(scope="genome_file", file_kind="genome", gap_type="file_gap",
            reason_code=code, recoverability=recov,
            n_ground_truth_rows=sum(gt_by_acc.get(r["genome_accession"], 0) for r in rs),
            n_strain_pairing_rows=sum(sp_by_acc.get(r["genome_accession"], 0) for r in rs),
            n_distinct_items=len(rs),
            items=_cap(r["genome_accession"] for r in rs),
            studies=_cap(r["study_id"][:8] for r in rs),
            detail=detail, next_action=action)

    noroute = [r for r in plan if r["downloaded"] != "yes" and r["route"] == "unresolved"]
    if noroute:
        add(scope="genome_file", file_kind="genome", gap_type="file_gap",
            reason_code="no_resolution_route", recoverability="manual",
            n_ground_truth_rows=sum(gt_by_acc.get(r["genome_accession"], 0) for r in noroute),
            n_strain_pairing_rows=sum(sp_by_acc.get(r["genome_accession"], 0) for r in noroute),
            n_distinct_items=len(noroute),
            items=_cap(r["genome_accession"] for r in noroute),
            studies=_cap(r["study_id"][:8] for r in noroute),
            detail="Neither the accession, the BioSample, nor the strain name "
                   "resolved to anything at NCBI.",
            next_action="Check the source publication for a deposition, or ask "
                        "the depositor. Some may be unreleased.")

    # ---- 6. JGI-only records, deferred by decision ----------------------
    todo = [r for r in jgi if r["already_have"] == "no"]
    if todo:
        add(scope="genome_file", file_kind="genome", gap_type="file_gap",
            reason_code="jgi_img_only_no_ncbi_mirror", recoverability="deferred",
            n_ground_truth_rows=sum(int(r["gt_links"]) for r in todo),
            n_strain_pairing_rows=sum(sp_by_acc.get(r["jgi_img_id"], 0) for r in todo),
            n_distinct_items=len(todo),
            items=_cap(r["jgi_img_id"] for r in todo),
            studies=_cap(r["study_id"][:8] for r in todo),
            detail="IMG Taxon OIDs with no NCBI equivalent. NPOmix used the "
                   "same JGI identifiers, so their supplement offers no "
                   "shortcut. DEFERRED 2026-08-27: study 297c364c already has "
                   "130 of its 157 genomes on disk, so these 7 do not change "
                   "rankability, and none is NPOmix-covered.",
            next_action="DEFERRED. If revisited: img.jgi.doe.gov, search by "
                        "Taxon OID, download 'Genbank (.gbk)' (not FASTA), "
                        "verify strain against genome_label. Login required.")
    return out


def main() -> int:
    inv = collect()
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(inv)

    lines = ["# Ground-truth data we do not have, and why", "",
             f"Generated from `ground_truth.csv` ({len(rows(GT))} links), "
             "`refetch_plan_all.csv` and `mibig_migration_report.csv`.", "",
             "A `mapping_gap` means the file is on disk but nothing joins it to "
             "the link. A `file_gap` means the file is not in hand.", ""]
    for group, title in (("mapping_gap", "Mapping gaps -- files exist, join missing"),
                         ("not_attempted", "Not attempted yet"),
                         ("file_gap", "File gaps -- not in hand")):
        sel = [g for g in inv if g["gap_type"] == group]
        if not sel:
            continue
        lines += [f"## {title}", ""]
        for g in sel:
            lines += [
                f"### {g['gap_id']} · {g['reason_code']} ({g['recoverability']})", "",
                f"- **file kind**: {g['file_kind']}",
                f"- **ground-truth rows affected**: {g['n_ground_truth_rows']}",
            ]
            if g["n_strain_pairing_rows"]:
                lines.append(
                    f"- **strain-pairing rows affected**: {g['n_strain_pairing_rows']} "
                    "(strain-axis depth, not links)")
            lines += [
            ]
            if g["n_distinct_items"]:
                lines.append(f"- **distinct items**: {g['n_distinct_items']} — `{g['items']}`")
            lines += [f"- **why**: {g['detail']}",
                      f"- **next**: {g['next_action']}", ""]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {OUT_CSV.relative_to(ROOT)}  ({len(inv)} gap classes)")
    print(f"wrote {OUT_MD.relative_to(ROOT)}\n")
    print(f"{'id':5} {'kind':13} {'reason':36} {'recov':12} {'gt':>4} {'pair':>5} {'items':>6}")
    for g in inv:
        print(f"{g['gap_id']:5} {g['file_kind']:13} {g['reason_code'][:35]:36} "
              f"{g['recoverability']:12} {g['n_ground_truth_rows']:>4} "
              f"{str(g['n_strain_pairing_rows']):>5} {str(g['n_distinct_items']):>6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
