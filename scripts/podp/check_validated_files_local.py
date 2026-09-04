"""Check how many validated genome-metabolome pairs have both files on disk.

"Validated" = tier-1 rows in results/podp/ground_truth.csv (experimentally
validated BGC-MS2 links: proof types "knockouts/heterologous expression" or
"NMR and/or detailed MS/MS analysis"). For each such row this checks:
  - genome file: path_genome resolves to an on-disk .gbk that actually
    contains sequence data (LOCUS header + ORIGIN block), matching
    build_ground_truth.classify_genome_file's definition of "sequence_ok".
    A .gbk that is metadata-only (no ORIGIN, the paired *.gbk.stub case) or
    absent does not count.
  - ms2 file: path_ms2_run (or path_ms2_spectrum) exists on disk.

ground_truth.csv's own path_genome/path_ms2_run columns still point at the
pre-rename tree (data/PoDP/ground_truth_paired_data); this script remaps
that prefix to the current data/processed/podp_ground_truth_paired_data
location before checking existence, and re-derives genome sequence status
from the file on disk rather than trusting the CSV's cached status_genome.

Usage:
    python check_validated_files_local.py [--ground-truth-csv PATH] [--tier N ...]
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

DEFAULT_GROUND_TRUTH_CSV = "/home/rosaliewang/projects/MASc-EDA/results/podp/ground_truth.csv"
OLD_PREFIX = "/home/rosaliewang/projects/MASc-EDA/data/PoDP/ground_truth_paired_data"
NEW_PREFIX = "/home/rosaliewang/projects/MASc-EDA/data/processed/podp_ground_truth_paired_data"


def remap_path(path):
    if not isinstance(path, str) or not path:
        return None
    return path.replace(OLD_PREFIX, NEW_PREFIX)


def classify_genome_file(path_str):
    if not isinstance(path_str, str) or not path_str:
        return "absent"
    path = Path(path_str)
    if not path.exists():
        return "absent"
    try:
        head = path.read_text(errors="replace")
    except OSError:
        return "absent"
    if not head.startswith("LOCUS"):
        return "not_genbank_error"
    if re.search(r"^ORIGIN", head, flags=re.MULTILINE):
        return "sequence_ok"
    return "metadata_only_no_contigs"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth-csv", default=DEFAULT_GROUND_TRUTH_CSV)
    parser.add_argument("--tier", type=int, nargs="*", default=[1],
                         help="Tier(s) to treat as 'validated' (default: 1, "
                              "experimentally validated only)")
    args = parser.parse_args()

    df = pd.read_csv(args.ground_truth_csv)

    df["path_genome_local"] = df["path_genome"].apply(remap_path)
    df["path_ms2_local"] = df["path_ms2_run"].apply(remap_path)
    if "path_ms2_spectrum" in df.columns:
        spectrum_local = df["path_ms2_spectrum"].apply(remap_path)
        df["path_ms2_local"] = df["path_ms2_local"].fillna(spectrum_local)

    df["genome_status_ondisk"] = df["path_genome_local"].apply(classify_genome_file)
    df["genome_present"] = df["genome_status_ondisk"] == "sequence_ok"
    df["ms2_present"] = df["path_ms2_local"].apply(
        lambda p: isinstance(p, str) and Path(p).exists()
    )
    df["both_present"] = df["genome_present"] & df["ms2_present"]

    validated = df[df["tier"].isin(args.tier)].copy()

    print(f"Ground truth: {args.ground_truth_csv}")
    print(f"Validated tier filter: {args.tier}")
    print(f"Validated links (rows): {len(validated)}")
    print(f"  genome present locally (sequence_ok .gbk): {validated['genome_present'].sum()}")
    print(f"  ms2 file present locally: {validated['ms2_present'].sum()}")
    print(f"  BOTH present (fully local pair): {validated['both_present'].sum()}")

    both = validated[validated["both_present"]]
    print(f"\nOf the fully-local rows:")
    print(f"  distinct genome files: {both['path_genome_local'].nunique()}")
    print(f"  distinct ms2 files: {both['path_ms2_local'].nunique()}")

    missing_genome = validated[~validated["genome_present"]]
    missing_ms2 = validated[~validated["ms2_present"]]
    stub_only = validated[validated["genome_status_ondisk"] == "metadata_only_no_contigs"]
    print(f"\nMissing genome (any reason): {len(missing_genome)}")
    print(f"  of which metadata-only stub (no ORIGIN sequence): {len(stub_only)}")
    print(f"  of which fully absent: {len(missing_genome) - len(stub_only)}")
    print(f"Missing ms2 file: {len(missing_ms2)}")

    out_cols = ["podp_id", "mibig_id", "genome_label", "genome_accession",
                "genome_status_ondisk", "genome_present", "ms2_present",
                "both_present", "path_genome_local", "path_ms2_local"]
    out_path = "/tmp/claude-1000/-home-rosaliewang-projects-NPLinker-Application/efd06ba5-18e3-4035-9c8f-243766959541/scratchpad/validated_local_file_status.csv"
    validated[out_cols].to_csv(out_path, index=False)
    print(f"\nWrote per-row detail to {out_path}")


if __name__ == "__main__":
    main()
