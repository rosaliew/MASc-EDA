#!/usr/bin/env python3
"""
00_inspect_columns.py
Run this FIRST, before the main analysis script.
Prints the exact column names in your GNPS output files so you can
fill in BROTH_LABEL and ISOLATE_LABELS correctly in 00_gnps_cmn_analysis.py
"""

from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1] / "data" / "MCE-GNPS"

def peek(path, label, nrows=3):
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"File: {path.name}")
    print(f"{'='*60}")
    try:
        df = pd.read_csv(path, sep="\t", nrows=nrows)
        print(f"Shape (first {nrows} rows): {df.shape}")
        print("\nAll columns:")
        for i, c in enumerate(df.columns):
            print(f"  [{i:03d}] {c}")
        print(f"\nFirst row sample:")
        print(df.iloc[0].to_string())
    except Exception as e:
        print(f"  ERROR reading file: {e}")

# Clustering
clustering_dir = BASE_DIR / "nf_output" / "clustering"
for f in clustering_dir.glob("*.tsv"):
    peek(f, "CLUSTERING / NODE TABLE")

# Networking
networking_dir = BASE_DIR / "nf_output" / "networking"
for f in networking_dir.glob("*.tsv"):
    peek(f, "NETWORKING / EDGE TABLE")

# Library summary
lib_dir = BASE_DIR / "nf_output" / "librarysummary"
for f in lib_dir.glob("*.tsv"):
    peek(f, "LIBRARY SUMMARY")

# Metadata
meta_dir = BASE_DIR / "metadata_filename"
for f in list(meta_dir.glob("*.txt")) + list(meta_dir.glob("*.tsv")):
    peek(f, "METADATA FILE")

print("\n\nDone. Use the column names above to fill in:")
print("  BROTH_LABEL     in 00_gnps_cmn_analysis.py")
print("  ISOLATE_LABELS  in 00_gnps_cmn_analysis.py")