# MASc-EDA

Exploratory data analysis for a MASc thesis on linking bacterial biosynthetic gene
clusters (BGCs) to the metabolites they produce.

The work here builds and interrogates the **evaluation data** that a genome–metabolome
linking tool is scored against — chiefly a ground-truth table of BGC↔MS2 links derived
from the [Paired omics Data Platform](https://pairedomicsdata.bioinformatics.nl/) (PoDP).

## Contents

```
data/        inputs and derived tables (see data/PoDP/README.md)
notebooks/   analysis notebooks, named <topic>-<YYYY-MM>-<subject>.ipynb
scripts/     reproducible pipeline steps (see scripts/legacy/README.md)
results/     figures and exported tables
tests/       unit tests
```

### The main artifact

`data/PoDP/ground_truth.csv` — **115 BGC↔MS2 links across 63 columns**, built from 76
PoDP JSON descriptors. Links are tiered by evidence strength, with InChIKey comparison
(not string matching) as the authoritative structural check.

Its companion `data/PoDP/strain_pairings.csv` holds 4,966 strain↔file rows. These record
*declared co-origin only* and are **not** ground truth — see `data/PoDP/README.md`.

## Datasets

| directory | what it is |
|---|---|
| `data/PoDP/` | PoDP descriptors, downloaded paired data, and the ground-truth tables |
| `data/MicroChemEco/` | in-house MicroChemEco GNPS and antiSMASH runs |
| `data/NPOmix_inputs/` | NPOmix validation set, for comparison against a published method |
| `data/MiBIG_4.0/` | MIBiG 4.0 reference bundle (gitignored; re-downloaded on demand) |

## Pipeline

```bash
# Build the ground-truth tables from the PoDP descriptors
python scripts/build_ground_truth.py

# Preview without writing, or diff MIBiG 3.1 -> 4.0
python scripts/build_ground_truth.py --dry-run
python scripts/build_ground_truth.py --migration-report

# Genome recovery: plan routes, fetch, then promote verified records
python scripts/refetch_genomes.py --scope all --plan
python scripts/refetch_genomes.py --scope all --download
python scripts/refetch_genomes.py --promote

# Diagnostics
python scripts/screen_nplinker_podp.py     # pre-flight screen for NPLinker's PoDP mode
python scripts/unavailable_inventory.py    # what we cannot obtain, and why
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
