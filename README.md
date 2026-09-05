# MASc-EDA

Curation of the publicly available data for a machine-learning problem: linking
bacterial biosynthetic gene clusters (BGCs) to the metabolites they produce.

**This repository is the dataset, not an evaluation of any one tool.** Its job is to
survey what exists across the published record and is readily machine-usable for that
problem — curated pairs of **GBK + mzML/mzXML**, together with the scripts that
locate, download and verify them — drawn from the
[Paired omics Data Platform](https://pairedomicsdata.bioinformatics.nl/) (PoDP).

The curated pairs are intended to outlive any single experiment: they are used to
benchmark NPLinker2 (in the sibling `NPLinker-Benchmarking` repo), will be used to
benchmark **other linking models**, and are the training data for models developed
later in this thesis. That is why curation and its download scripts live here, on
their own, rather than inside a tool evaluation.

Because curation is about what is *currently known*, filtering runs against the most
up-to-date BGC records — **MIBiG 4.0**. A tool being benchmarked may prepare its own
references at a different version from inside its own pipeline; that is the tool's
business and is not reconciled here.

## Contents

```
data/raw/podp_database/                   76 PoDP JSON descriptors -- the source of truth
data/processed/podp_ground_truth_paired_data/   downloaded genomes + MS2 runs, one dir per study
data/external/                            MiBIG_4.0, NPOmix_inputs (gitignored, re-downloadable)
data/interim/                             MCE GNPS + antiSMASH working data
notebooks/{podp,mce}/                     analysis notebooks
scripts/{podp,mce}/                       reproducible pipeline steps
results/{podp,mce}/                       the exported tables and figures
tests/                                    unit tests
```

Restructured 2026-09 from a single `data/PoDP/` tree. Scripts resolve data
through candidate lists that try the current path first and the pre-reorg path
behind it, so both layouts work; when moving something, add the new path to the
front of the list rather than replacing the entry.

### The main artifact

`results/podp/curated_ground_truth_links.csv` — **the 64 curated experimentally-validated
links**. This is the evaluation set the benchmark scores against, and it is consumed
downstream by the sibling `NPLinker-Benchmarking` repo, which vendors a copy at a pinned
commit of this one (see its `data/external/podp-ground-truth/SOURCE.md`).

Its unfiltered source `results/podp/podp_bgc_ms2_links.csv` — **115 BGC↔MS2 links across
63 columns**, built from 76 PoDP JSON descriptors. Links are tiered by evidence strength, with InChIKey comparison
(not string matching) as the authoritative structural check.

Its companion `results/podp/podp_genome_metabolome_pairs.csv` holds 4,966 genome↔metabolome file-pair rows. These record
*declared co-origin only* and are **not** ground truth — see `data/PoDP/README.md`.

## Datasets

| directory | what it is |
|---|---|
| `data/raw/podp_database/` | the 76 PoDP JSON descriptors, verbatim |
| `data/processed/podp_ground_truth_paired_data/` | downloaded genomes + MS2 runs |
| `data/interim/MCE_*/` | in-house MicroChemEco GNPS and antiSMASH runs |
| `data/external/NPOmix_inputs/` | NPOmix validation set, for comparison against a published method |
| `data/external/MiBIG_4.0/` | MIBiG 4.0 reference bundle (gitignored; re-downloaded on demand) |

## Pipeline

```bash
# Run from the repo root, as modules -- the scripts import each other by package
# path, so `python scripts/podp/build_ground_truth.py` fails on ModuleNotFoundError.

# Build the ground-truth tables from the PoDP descriptors
python -m scripts.podp.build_ground_truth

# Preview without writing, or diff MIBiG 3.1 -> 4.0
python -m scripts.podp.build_ground_truth --dry-run
python -m scripts.podp.build_ground_truth --migration-report

# Build the two curated datasets and assert every headline count
python -m scripts.podp.build_curated_datasets --check

# Genome recovery: plan routes, fetch, then promote verified records
python -m scripts.podp.refetch_genomes --scope all --plan
python -m scripts.podp.refetch_genomes --scope all --download
python -m scripts.podp.refetch_genomes --promote

# Diagnostics
python -m scripts.podp.screen_nplinker_podp   # pre-flight screen for NPLinker's PoDP mode
python -m scripts.podp.unavailable_inventory  # what we cannot obtain, and why
```

The `msms*.py` scripts are a separate, numbered MS/MS cleaning and EDA sequence; run
`msms00_inspect_columns.py` first to discover the column labels the later steps need.

## Requirements

Python 3.11+. `pandas`, `numpy`, `matplotlib`, `networkx`, `httpx`, `beautifulsoup4`,
plus `rdkit` for the structural comparison in `build_ground_truth.py`.

```bash
pytest tests/
```

## What is not in this repository

Raw instrument data, genome sequences and pipeline dumps stay on disk but out of git —
the working tree is ~13 GB against ~10 MB tracked. The rules in `.gitignore` are keyed to
file type and directory name at any depth, so new download folders are covered without
editing it.

Everything excluded is either re-downloadable from a recorded accession or rebuildable by
the scripts above from a recorded input and tool version. `data/PoDP/reports/refetch_plan_all.csv`
is the provenance record for every genome: which route resolved it, whether the strain
verified, and why anything was held back.


# data/PoDP layout

```
podp_bgc_ms2_links.csv            115 BGC<->MS2 links, 63 cols. UNFILTERED source (gold+silver+bronze).
podp_genome_metabolome_pairs.csv  4,966 declared genome<->metabolome pairs. Co-origin only, NOT validated.
curated_ground_truth_links.csv      64 curated experimentally-validated links (the evaluation set).
curated_paired_dataset.csv       3,585 curated pairs (53 validated_pairs + 3,532 unvalidated_pairs).

json_descriptors/         77 PoDP descriptors, verbatim. The source of truth.
ground_truth_paired_data/ downloaded genomes + MS2 runs, one folder per study_id.
                          *.gbk.stub = a broken file kept beside its replacement.
gnps_cache/               15 GNPS tasks: clusterinfo.tsv, params.tsv, summary.tsv.
antismash/                antiSMASH output (one study so far).

reports/                  derived diagnostics -- regenerable, safe to delete.
legacy/                   superseded artifacts kept for provenance.
```

## reports/

| file | written by |
|---|---|
| `mibig_migration_report.csv` | `build_ground_truth.py --mibig-migration` |
| `refetch_plan_all.csv` | `refetch_genomes.py --scope all --plan` |
| `jgi_worklist_297c364c.csv` | ad hoc; the JGI-only genomes |
| `nplinker_podp_screen.csv` | `screen_nplinker_podp.py` |
| `unavailable_data_inventory.csv` / `.md` | `unavailable_inventory.py` |

`refetch_plan_all.csv` is the **provenance record for every genome**: which route
resolved it, whether the strain verified, and why anything was held back. The
`genome_refetch/` staging tree it describes was deleted on 2026-08-27 once
everything was promoted -- 226 of its 250 files were byte-identical duplicates of
promoted files, and the other 24 were quarantined (wrong strain, metagenome-scale,
no sequence). All 24 remain listed in this CSV and are re-downloadable with
`refetch_genomes.py --scope all --download`.

## legacy/

`verified_bgc_ms2_manifest.json` and `verified_bgc_ms2_genome_manifest.json`,
superseded by `podp_bgc_ms2_links.csv`. Verified 2026-08-27: all 67 links in both
appear in `podp_bgc_ms2_links.csv`, matched on project + MS2 scan + run URL. The two
scripts that wrote them are in `scripts/legacy/`.

## gnps_cache/ is more useful than it looks

`gnps_cache/<task>/summary.tsv` is **byte-identical** to the `file_mappings.tsv`
that NPLinker's `arrange_gnps()` extracts from a GNPS task archive (verified by
md5 on task `c22f44b14a3d450eb836d607cb9521bb`). Its columns after the first six
are strain names -- it is a spectrum x strain occurrence matrix. All 15 cached
tasks are exactly the 15 that `nplinker_podp_screen.csv` covers, so the GNPS
download step is already done locally.
