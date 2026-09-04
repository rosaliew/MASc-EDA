"""Check overlap between genome-metabolome pairs and BGC-MS2 pairs in the PoDP database.

Genome-metabolome pairs are keyed on (genome_accession, mzml/mzxml URL), where
genome_accession is resolved per-study from `genome_label` -> the best available
accession in that study's `genomes` list (GenBank > RefSeq > ENA/NCBI > JGI).

BGC-MS2 links in the raw JSON never carry a genome_label or NCBI accession
(verified across all 76 study files) -- they only carry an MS2_URL and a MIBiG
BGC_ID. So a true (genome_accession, url) key can't be built on that side.
Overlap is therefore checked on the mzml/mzxml URL alone: any URL used both as
a genome-metabolome metabolomics_file and as a BGC_MS2_links MS2_URL is
reported as a duplicate, with genome accession context pulled in from the
genome-metabolome side.

Usage:
    python check_pair_overlap.py [--database-dir DIR] [--out-csv PATH]
"""

import argparse
import csv
import glob
import json
import os
import sys

ACCESSION_FIELDS = ["GenBank_accession", "RefSeq_accession", "ENA_NCBI_accession", "JGI_Genome_ID"]

DEFAULT_DATABASE_DIR = "/home/rosaliewang/projects/MASc-EDA/data/raw/podp_database"


def best_accession(genome_id_info):
    for field in ACCESSION_FIELDS:
        value = genome_id_info.get(field)
        if value:
            return value, field
    return None, None


def load_study(file_path):
    with open(file_path) as f:
        return json.load(f)


def build_genome_lookup(study):
    lookup = {}
    for genome in study.get("genomes") or []:
        label = genome.get("genome_label")
        if not label:
            continue
        accession, accession_field = best_accession(genome.get("genome_ID") or {})
        lookup[label] = {
            "biosample_accession": genome.get("BioSample_accession"),
            "accession": accession,
            "accession_field": accession_field,
        }
    return lookup


def collect_genome_metabolome_pairs(study, study_id, genome_lookup):
    pairs = []
    for link in study.get("genome_metabolome_links") or []:
        genome_label = link.get("genome_label")
        url = link.get("metabolomics_file")
        if not url:
            continue
        genome_info = genome_lookup.get(genome_label, {})
        pairs.append({
            "study_identifier": study_id,
            "genome_label": genome_label,
            "biosample_accession": genome_info.get("biosample_accession"),
            "genome_accession": genome_info.get("accession"),
            "genome_accession_field": genome_info.get("accession_field"),
            "mzml_url": url,
        })
    return pairs


def collect_bgc_ms2_urls(study, study_id):
    entries = []
    for link in study.get("BGC_MS2_links") or []:
        url = link.get("MS2_URL")
        if not url:
            continue
        mibig_id = (link.get("BGC_ID") or {}).get("MIBiG_number")
        entries.append({
            "study_identifier": study_id,
            "mibig_id": mibig_id,
            "mzml_url": url,
        })
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-dir", default=DEFAULT_DATABASE_DIR,
                         help="Directory of PoDP study JSON files")
    parser.add_argument("--out-csv", default=None,
                         help="Optional path to write overlap rows as CSV")
    args = parser.parse_args()

    json_files = sorted(glob.glob(os.path.join(args.database_dir, "*.json")))
    if not json_files:
        print(f"No JSON files found in {args.database_dir}", file=sys.stderr)
        sys.exit(1)

    all_gm_pairs = []
    all_bgc_ms2 = []

    for study_idx, file_path in enumerate(json_files):
        study_id = f"Study_{study_idx}"
        study = load_study(file_path)
        genome_lookup = build_genome_lookup(study)
        all_gm_pairs.extend(collect_genome_metabolome_pairs(study, study_id, genome_lookup))
        all_bgc_ms2.extend(collect_bgc_ms2_urls(study, study_id))

    bgc_ms2_urls = {}
    for entry in all_bgc_ms2:
        bgc_ms2_urls.setdefault(entry["mzml_url"], []).append(entry)

    overlap_rows = []
    for pair in all_gm_pairs:
        matches = bgc_ms2_urls.get(pair["mzml_url"])
        if not matches:
            continue
        for match in matches:
            overlap_rows.append({
                "mzml_url": pair["mzml_url"],
                "gm_study_identifier": pair["study_identifier"],
                "genome_label": pair["genome_label"],
                "biosample_accession": pair["biosample_accession"],
                "genome_accession": pair["genome_accession"],
                "genome_accession_field": pair["genome_accession_field"],
                "bgc_study_identifier": match["study_identifier"],
                "mibig_id": match["mibig_id"],
                "same_study": pair["study_identifier"] == match["study_identifier"],
            })

    total_gm_pairs = len(all_gm_pairs)
    distinct_gm_urls = len({p["mzml_url"] for p in all_gm_pairs})
    total_bgc_ms2 = len(all_bgc_ms2)
    distinct_bgc_urls = len(bgc_ms2_urls)
    distinct_overlap_urls = len({r["mzml_url"] for r in overlap_rows})

    print(f"Studies scanned: {len(json_files)}")
    print(f"Genome-metabolome pairs: {total_gm_pairs} ({distinct_gm_urls} distinct URLs)")
    print(f"BGC-MS2 links with a URL: {total_bgc_ms2} ({distinct_bgc_urls} distinct URLs)")
    print(f"Overlapping URLs (used in both): {distinct_overlap_urls}")
    print(f"Overlap rows (pair x link matches): {len(overlap_rows)}")

    if overlap_rows:
        same_study = sum(1 for r in overlap_rows if r["same_study"])
        cross_study = len(overlap_rows) - same_study
        print(f"  same-study matches: {same_study}, cross-study matches: {cross_study}")
        print()
        print("Sample overlap rows:")
        for row in overlap_rows[:10]:
            print(f"  [{row['gm_study_identifier']} / {row['bgc_study_identifier']}] "
                  f"genome={row['genome_label']!r} acc={row['genome_accession']} "
                  f"mibig={row['mibig_id']} url={row['mzml_url']}")

    if args.out_csv:
        fieldnames = ["mzml_url", "gm_study_identifier", "genome_label", "biosample_accession",
                      "genome_accession", "genome_accession_field", "bgc_study_identifier",
                      "mibig_id", "same_study"]
        with open(args.out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(overlap_rows)
        print(f"\nWrote {len(overlap_rows)} overlap rows to {args.out_csv}")


if __name__ == "__main__":
    main()
