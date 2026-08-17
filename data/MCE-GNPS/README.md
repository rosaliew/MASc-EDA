# masc-utils: GNPS CMN Analysis Pipeline

**MicroChemEco Lab | R. Wang**

---

## Setup

```bash
# 1. Clone/copy these scripts into your GNPS output directory (GNPS_CMN_ALL/)
# 2. Create the conda environment
bash setup_env.sh

# 3. Activate it (every new terminal session)
conda activate masc-utils
```

---

## Directory structure expected

```
GNPS_CMN_ALL/
├── nf_output/
│   ├── clustering/         ← clusterinfo_summary.tsv  (NODE table)
│   ├── networking/         ← networking_pairs.tsv      (EDGE table)
│   ├── librarysummary/     ← library match hits
│   ├── library/
│   ├── library_intermediate/
│   ├── metadata/
│   └── temp_pairs/
├── input_spectra/
├── input_libraries/
├── metadata_filename/      ← MSV000090080-metadata.txt
├── setup_env.sh
├── 00_inspect_columns.py   ← run FIRST
├── 00_gnps_cmn_analysis.py ← main pipeline
└── analysis_output/        ← created automatically
```

---

## Execution order

### Step 1 — Inspect your columns
```bash
python 00_inspect_columns.py 2>&1 | tee inspect_output.txt
```
This prints every column in every relevant file. Use it to find:
- The exact name of your broth sample column → set `BROTH_LABEL`
- The exact names of all 11 isolate columns → set `ISOLATE_LABELS`

### Step 2 — Edit the config block in `00_gnps_cmn_analysis.py`
```python
BROTH_LABEL = "marine_broth"   # exact column name stem from Step 1
ISOLATE_LABELS = [
    "Shewanella_B2R08",        # all 11 isolate column name stems
    "Vibrio_1A01",
    # ...
]
```

### Step 3 — Run the main pipeline
```bash
python 00_gnps_cmn_analysis.py 2>&1 | tee analysis_log.txt
```

### Outputs (in `analysis_output/`)
| File | Contents |
|------|----------|
| `clusterinfo_filtered.tsv` | Node table with broth-only nodes removed + `origin` column |
| `networking_pairs_filtered.tsv` | Edge list for retained nodes |
| `library_hits_validated.tsv` | Library matches on filtered nodes |
| `cmn_piechart_network.png` | Network with 11-color isolate pie charts |
| `cmn_piechart_network.svg` | Same, vector format (for Illustrator/Inkscape) |
| `annotation_status_pie.png` | Annotated vs. unannotated fraction |

---

## Important parameter notes (your job)

| Parameter | Value | Implication |
|-----------|-------|-------------|
| `networking_min_cosine` | 0.7 | Conservative — edges are reliable |
| `networking_min_matched_peaks` | 6 | Conservative — reduces spurious edges |
| `min_cluster_size` | 2 | Singletons excluded from network |
| `topology_maxcomponent` | 100 | Components capped at 100 nodes — flag any exactly-100-node clusters |
| `library_min_cosine` | 0.7 | Matches ≥ 0.7 cosine only |
| `clustering_tool` | falcon | Fast density-based consensus spectra |
| `pm_tolerance` | 2.0 Da / 20 ppm | Standard for Q-TOF/Orbitrap |

---

## Next steps after this pipeline

1. **SIRIUS + CSI:FingerID** — run on the filtered `.mgf` to annotate unannotated nodes by chemical class (CANOPUS)
2. **microbeMASST** — search spectral features against the microbial spectral database
3. **Cytoscape** — import filtered `.graphml` for publication-quality layout; use `EnrichmentMap` for pie charts
4. **FBMN re-run** — when you have intensity data (feature table from MZmine/XCMS), redo as FBMN to get quantitative isolate-level abundance per node
5. **NPOmix/NPLinker2** — feed filtered node table + BiG-SCAPE GCF table into linking pipeline
