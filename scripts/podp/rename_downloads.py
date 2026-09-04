#!/usr/bin/env python3
"""Rename the PoDP download tree onto keys that join to ground_truth.csv.

Target layout, every component of which is a column in ground_truth.csv::

    <study_id>/                                     # main study's podp_id (.N kept)
      STUDY_MEMBERS.txt                             # only when >1 project shares it
      genomes/<genome_accession>.gbk                # genome_accession
      ms2/<massive_id>__<genome_accession>__<original>   # gnps_massive_id,
                                                    # genome_accession, and the
                                                    # verbatim tail of ms2_run_url

The original filename is preserved verbatim on the end. The previous rename pass
sanitised punctuation (``-`` to ``_``, parentheses dropped) which made 20 of 67
files impossible to trace back to their URL; keeping the tail intact means the
join is exact in both directions and cannot silently drift again.

Folders are named for the *study*, not the project record. Where PoDP carries one
study as two project records -- same PI, identical MassIVE depositions -- both
share a folder named for whichever record holds the majority of links, and
STUDY_MEMBERS.txt inside names every record it covers. The original podp_id is
never lost; it stays in the podp_id column, with study_id alongside it.

Nothing moves without ``--apply``. A reversible manifest is always written first.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # repo root (this file is scripts/podp/*)
GROUND_TRUTH = ROOT / "results" / "podp" / "podp_bgc_ms2_links.csv"
DOWNLOAD_DIR = ROOT / "data" / "PoDP" / "ground_truth_paired_data"
MANIFEST = DOWNLOAD_DIR / "file_rename_manifest.csv"

MS2_SUFFIXES = (".mzML", ".mzXML")


def original_basename(url: str) -> str:
    """Verbatim final path segment of an MS2 URL."""
    return url.rstrip("/").split("/")[-1]


def build_plan(rows: list[dict], download_dir: Path) -> list[dict]:
    """Return [{kind, old, new}] with no duplicates and no collisions."""
    on_disk: dict[str, Path] = {}
    for path in download_dir.rglob("*"):
        if path.is_file():
            on_disk.setdefault(path.name, path)

    plan: dict[Path, dict] = {}
    for row in rows:
        accession = row["genome_accession"]
        # study_id, not podp_id: merged records share the main study's folder
        project = download_dir / (row.get("study_id") or row["podp_id"])

        # keyed on the accession, not file_genome: broken records (WGS stubs,
        # efetch error pages) have a blank file_genome but are still real files
        # that need organising -- they are the re-fetch worklist.
        if accession:
            old = on_disk.get(f"{accession}.gbk")
            if old is not None:
                new = project / "genomes" / f"{accession}.gbk"
                plan[old] = {"kind": "genome", "old": old, "new": new}

        if row["file_ms2_run"]:
            old = on_disk.get(row["file_ms2_run"])
            if old is not None:
                tail = original_basename(row["ms2_run_url"])
                stem = f"{row['gnps_massive_id']}__{accession}__{tail}"
                plan[old] = {"kind": "ms2", "old": old, "new": project / "ms2" / stem}

    # a genome record can be cited by several rows; keep one entry per source file
    entries = list(plan.values())

    targets: dict[Path, list[Path]] = defaultdict(list)
    for entry in entries:
        targets[entry["new"]].append(entry["old"])
    for target, sources in targets.items():
        if len(sources) > 1:
            raise SystemExit(f"collision: {len(sources)} files map to {target}")
    return entries


def write_study_markers(rows: list[dict], download_dir: Path) -> list[Path]:
    """Drop STUDY_MEMBERS.txt into every folder that covers more than one record."""
    seen: dict[str, str] = {}
    for row in rows:
        study_id = row.get("study_id") or row["podp_id"]
        seen[study_id] = row.get("study_members") or study_id

    written = []
    for study_id, members in seen.items():
        parts = members.split("|")
        folder = download_dir / study_id
        if len(parts) < 2 or not folder.is_dir():
            continue
        marker = folder / "STUDY_MEMBERS.txt"
        marker.write_text(
            "This folder holds one study that PoDP records under several project ids.\n"
            "They reference an identical set of MassIVE depositions.\n\n"
            f"folder / main record : {study_id}\n"
            + "".join(f"also covers          : {p}\n" for p in parts if p != study_id)
            + "\nFiles are shared across every record listed above. In ground_truth.csv\n"
              "and strain_pairings.csv the original podp_id is preserved; study_id names\n"
              "this folder and study_role says which record is the main one.\n"
        )
        written.append(marker)
    return written


def unplaced_files(download_dir: Path, plan: list[dict]) -> list[Path]:
    """Data files on disk that the plan does not account for."""
    planned = {entry["old"] for entry in plan}
    out = []
    for path in download_dir.rglob("*"):
        if not path.is_file() or path in planned:
            continue
        if path.suffix in MS2_SUFFIXES or path.suffix == ".gbk":
            out.append(path)
    return out


def write_manifest(plan: list[dict], manifest: Path) -> None:
    """Record old -> new before touching anything, so the move is reversible."""
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["kind", "old_path", "new_path", "bytes"])
        for entry in sorted(plan, key=lambda e: str(e["new"])):
            writer.writerow([
                entry["kind"], str(entry["old"]), str(entry["new"]),
                entry["old"].stat().st_size if entry["old"].exists() else "",
            ])


def apply_plan(plan: list[dict]) -> int:
    moved = 0
    for entry in plan:
        old, new = entry["old"], entry["new"]
        if not old.exists():
            continue
        if new.exists():
            print(f"  skip (target exists): {new.name}")
            continue
        new.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(old), str(new))
        moved += 1
    return moved


def prune_empty(download_dir: Path) -> list[Path]:
    """Remove directories left empty by the move, deepest first."""
    removed = []
    for path in sorted(download_dir.rglob("*"), key=lambda p: len(p.parts), reverse=True):
        if path.is_dir() and not any(path.iterdir()):
            path.rmdir()
            removed.append(path)
    return removed


def revert(manifest: Path) -> int:
    """Undo a previous run using its manifest."""
    if not manifest.exists():
        raise SystemExit(f"no manifest at {manifest}")
    moved = 0
    with manifest.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            new, old = Path(row["new_path"]), Path(row["old_path"])
            if new.exists() and not old.exists():
                old.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(new), str(old))
                moved += 1
    return moved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--apply", action="store_true",
                        help="perform the moves (default is a dry run)")
    parser.add_argument("--revert", action="store_true",
                        help="undo a previous run from its manifest")
    parser.add_argument("--ground-truth", type=Path, default=GROUND_TRUTH)
    parser.add_argument("--download-dir", type=Path, default=DOWNLOAD_DIR)
    args = parser.parse_args()

    if args.revert:
        print(f"reverted {revert(MANIFEST)} files")
        return 0

    if not args.ground_truth.exists():
        print(f"error: {args.ground_truth} not found -- build it first", file=sys.stderr)
        return 1

    with args.ground_truth.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    plan = build_plan(rows, args.download_dir)
    stray = unplaced_files(args.download_dir, plan)

    kinds = defaultdict(int)
    for entry in plan:
        kinds[entry["kind"]] += 1
    print(f"planned moves: {len(plan)}  ({dict(kinds)})")
    print(f"data files not covered by the plan: {len(stray)}")
    for path in stray[:10]:
        print(f"   {path.relative_to(args.download_dir)}")

    print("\nsample:")
    for entry in sorted(plan, key=lambda e: str(e["new"]))[:6]:
        print(f"   {entry['old'].name}")
        print(f"     -> {entry['new'].relative_to(args.download_dir)}")

    write_manifest(plan, MANIFEST)
    print(f"\nmanifest written: {MANIFEST}")

    if not args.apply:
        print("dry run -- nothing moved. re-run with --apply")
        return 0

    moved = apply_plan(plan)
    pruned = prune_empty(args.download_dir)
    markers = write_study_markers(rows, args.download_dir)
    print(f"\nmoved {moved} files; removed {len(pruned)} empty directories")
    for marker in markers:
        print(f"  study marker: {marker.relative_to(args.download_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
