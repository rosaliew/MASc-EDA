#!/usr/bin/env python3
"""Build the two curated PoDP datasets from the raw project records.

    curated_ground_truth_links.csv   64 links  -- the evaluation set
    curated_paired_dataset.csv     3585 pairs  -- the genome<->metabolome pairs

Run with --check to assert the headline counts instead of writing files; that
is the mode CI/the paper should use, so a change in the source records fails
loudly rather than silently moving a number in Section 3.

Why this script exists: notebooks/podp/podp.ipynb computed these numbers from
a Colab path with a genome-accession key bug (it looked for GenBank_ID /
RefSeq_ID / ENA_ID, but PoDP writes *_accession), which undercounted
retrievable genomes by ~2.5x. The notebook is fixed, but the definitions
belong somewhere importable rather than only in a notebook cell.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "podp_database"
OUT = REPO / "results" / "podp"

# Kept in sync with NPLinker-Application/nibi/bin/ground_truth_defs.py --
# that module is the single source of truth for the vocabulary; these are the
# same rules applied to the raw records rather than to ground_truth.csv.
GENOME_ID_KEYS = (
    "GenBank_accession", "RefSeq_accession", "ENA_NCBI_accession",
    "JGI_Genome_ID", "JGI_IMG_genome_ID", "JGI_ID",
)
CURATED_GENOME_TYPES = {"genome", "metagenome-assembled genome"}
CURATED_MS_FILE_TYPES = {"mzml", "mzxml"}
EXPERIMENTAL_PROOFS = {
    "Experimentally validated with knockouts, heterologous expression, "
    "or other gene cluster manipulation",
    "Experimentally validated with NMR and/or detailed MS/MS analysis",
}
MIBIG_PLACEHOLDER = "MIBiG number of a similar BGC"

EXPECTED = {
    "gt_links": 64, "gt_projects": 21, "gt_bgcs": 34,
    "gt_file_pairs": 53, "paired_file_pairs": 3585,
}


def _basename(url: str) -> str:
    return os.path.basename(urlparse(str(url)).path).lower()


def _ftype(url: str) -> str | None:
    bn = _basename(url)
    return bn.rsplit(".", 1)[-1] if "." in bn else None


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_rows, link_rows = [], []

    for fp in sorted(RAW.glob("*.json")):
        rec = json.loads(fp.read_text())
        podp_id = fp.name.rsplit(".", 1)[0]
        massive = ((rec.get("metabolomics") or {}).get("project") or {}).get("GNPSMassIVE_ID")
        exp = rec.get("experimental") or {}
        declared_extraction = {m.get("extraction_method") for m in (exp.get("extraction_methods") or [])}
        declared_instrument = {m.get("instrumentation_method") for m in (exp.get("instrumentation_methods") or [])}
        declared_prep = {m.get("sample_preparation_method") for m in (exp.get("sample_preparation") or [])}

        genomes = {}
        for g in rec.get("genomes") or []:
            gid = g.get("genome_ID") or {}
            genomes[g.get("genome_label")] = {
                "has_accession": any(gid.get(k) for k in GENOME_ID_KEYS),
                "genome_type": gid.get("genome_type"),
            }

        for link in rec.get("genome_metabolome_links") or []:
            label = link.get("genome_label")
            meta_file = link.get("metabolomics_file") or ""
            info = genomes.get(label)
            if not (info and label and meta_file):
                continue
            if info["genome_type"] not in CURATED_GENOME_TYPES:
                continue
            if not info["has_accession"]:
                continue
            if _ftype(meta_file) not in CURATED_MS_FILE_TYPES:
                continue
            # same source: the MS file must live under this study's own dataset
            if not (massive and massive in meta_file):
                continue
            pair_rows.append({
                "podp_id": podp_id,
                "genome_label": label,
                "genome_type": info["genome_type"],
                "metabolomics_file": meta_file,
                "gnps_massive_id": massive,
                # same experiment: each label resolves within this same study
                "same_experiment": (
                    link.get("extraction_method_label") in declared_extraction
                    and link.get("instrumentation_method_label") in declared_instrument
                    and link.get("sample_preparation_label") in declared_prep
                ),
            })

        study_genome_ok = any(
            v["has_accession"] and v["genome_type"] in CURATED_GENOME_TYPES
            for v in genomes.values()
        )
        for bl in rec.get("BGC_MS2_links") or []:
            proofs = bl.get("verification") or []
            if isinstance(proofs, str):
                proofs = [proofs]
            bgc_id = bl.get("BGC_ID") or {}
            mibig = bgc_id.get("MIBiG_number")
            if mibig == MIBIG_PLACEHOLDER:
                mibig = None
            # per-link pointer: URL *and* scan, so it resolves to ONE spectrum
            has_pointer = bool((bl.get("MS2_URL") or "").strip()) and bool(
                str(bl.get("MS2_scan") or "").strip())
            genome_ok = bool(mibig) or bool(bgc_id.get("ncbi_genome_accession")) or study_genome_ok
            if not (set(proofs) & EXPERIMENTAL_PROOFS and has_pointer and genome_ok):
                continue
            link_rows.append({
                "podp_id": podp_id,
                "mibig_id": mibig,
                "link_type": bl.get("link"),
                "ms2_url": bl.get("MS2_URL"),
                "ms2_scan": bl.get("MS2_scan"),
                "gnps_massive_id": massive,
            })

    pairs = pd.DataFrame(pair_rows).drop_duplicates(
        subset=["podp_id", "genome_label", "metabolomics_file"])
    links = pd.DataFrame(link_rows)

    # attach the backing file pair(s) so the link table carries its own unit
    by_study = pairs.groupby("podp_id")["metabolomics_file"].apply(
        lambda s: {_basename(x) for x in s}).to_dict()
    links["ms2_file_in_curated_pairs"] = [
        _basename(u) in by_study.get(pid, set())
        for pid, u in zip(links["podp_id"], links["ms2_url"])
    ]

    # Label every pair with its slice so the vocabulary travels with the data
    # and a downstream reader cannot mix the two up. "validated" means a
    # curated (experimentally-validated) BGC<->MS2 link sits on this pair.
    link_keys = set(zip(links["podp_id"], (_basename(u) for u in links["ms2_url"])))
    pairs["slice"] = [
        "validated_pairs" if (pid, _basename(mf)) in link_keys else "unvalidated_pairs"
        for pid, mf in zip(pairs["podp_id"], pairs["metabolomics_file"])
    ]
    return pairs, links


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="assert headline counts, write nothing")
    a = ap.parse_args()

    pairs, links = build()
    # The file pairs BACKING the curated links: distinct (study, genome,
    # metabolome file) rows, which is fewer than the link count because
    # several links can cite the same spectrum file.
    backing = pairs.merge(
        pd.DataFrame({
            "podp_id": links["podp_id"],
            "_ms2": [_basename(u) for u in links["ms2_url"]],
        }).drop_duplicates(),
        left_on=["podp_id", pairs["metabolomics_file"].map(_basename)],
        right_on=["podp_id", "_ms2"],
        how="inner",
    ).drop_duplicates(subset=["podp_id", "genome_label", "metabolomics_file"])

    got = {
        "gt_links": len(links),
        "gt_projects": links["podp_id"].nunique(),
        "gt_bgcs": links["mibig_id"].nunique(),
        "gt_file_pairs": len(backing),
        "paired_file_pairs": len(pairs),
    }

    print("CURATED GROUND TRUTH LINKS (curated experimentally-validated links)")
    print(f"  links                 : {got['gt_links']}")
    print(f"  projects              : {got['gt_projects']}")
    print(f"  distinct MIBiG BGCs   : {got['gt_bgcs']}")
    print(f"    with MIBiG          : {links['mibig_id'].notna().sum()}")
    print(f"    without MIBiG       : {links['mibig_id'].isna().sum()}  (kept by design)")
    print(f"  containment in pairs  : {int(links['ms2_file_in_curated_pairs'].sum())}/{got['gt_links']}")
    print(f"  backing file pairs    : {got['gt_file_pairs']}")
    print(f"  backing isolates      : {backing['genome_label'].nunique()}")
    print()
    # The three slices, named per the vocabulary in
    # NPLinker-Application/nibi/bin/ground_truth_defs.py. validated_pairs and
    # unvalidated_pairs differ ONLY by whether a validated link sits on top,
    # so the differentiator goes first in both names.
    n_val = got["gt_file_pairs"]
    n_unval = got["paired_file_pairs"] - n_val
    print("THE THREE DATA SLICES")
    print(f"  curated_links         : {got['gt_links']:5d} links")
    print(f"  validated_pairs       : {n_val:5d} pairs   backing those links")
    print(f"  unvalidated_pairs     : {n_unval:5d} pairs   no validated link "
          f"-- training pool, never ground truth")
    print(f"  discovery_pairs       :    (MCE, separate cohort -- not built here)")
    print()
    print("CURATED GENOME-METABOLOME PAIRED DATASET (curated paired dataset)")
    print(f"  data file pairs       : {got['paired_file_pairs']}"
          f"   = {n_val} validated + {n_unval} unvalidated")
    print(f"    genome-typed        : {(pairs.genome_type == 'genome').sum()}")
    print(f"    MAG-typed           : {(pairs.genome_type != 'genome').sum()}")
    print(f"  studies               : {pairs['podp_id'].nunique()}")
    print(f"  isolates              : {pairs['genome_label'].nunique()}")
    print(f"  metabolome files      : {pairs['metabolomics_file'].nunique()}")
    print(f"  same-experiment       : {int(pairs['same_experiment'].sum())}/{len(pairs)}")

    bad = [f"{k}: expected {EXPECTED[k]}, got {got[k]}"
           for k in EXPECTED if EXPECTED[k] != got[k]]
    if not links["ms2_file_in_curated_pairs"].all():
        bad.append("containment broken: some curated links are not curated pairs")
    if not pairs["same_experiment"].all():
        bad.append("some curated pairs are not same-experiment")

    if bad:
        print("\nFAILED:")
        for b in bad:
            print("  " + b)
        return 1
    print("\nOK: all headline counts and invariants match.")

    if not a.check:
        OUT.mkdir(parents=True, exist_ok=True)
        links.to_csv(OUT / "curated_ground_truth_links.csv", index=False)
        pairs.to_csv(OUT / "curated_paired_dataset.csv", index=False)
        print(f"wrote {OUT/'curated_ground_truth_links.csv'}")
        print(f"wrote {OUT/'curated_paired_dataset.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
