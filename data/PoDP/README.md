# data/PoDP layout

```
ground_truth.csv          115 BGC<->MS2 links, 63 columns. The evaluation table.
strain_pairings.csv       4,966 strain<->file rows. Declared co-origin, NOT ground truth.

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
superseded by `ground_truth.csv`. Verified 2026-08-27: all 67 links in both
appear in `ground_truth.csv`, matched on project + MS2 scan + run URL. The two
scripts that wrote them are in `scripts/legacy/`.

## gnps_cache/ is more useful than it looks

`gnps_cache/<task>/summary.tsv` is **byte-identical** to the `file_mappings.tsv`
that NPLinker's `arrange_gnps()` extracts from a GNPS task archive (verified by
md5 on task `c22f44b14a3d450eb836d607cb9521bb`). Its columns after the first six
are strain names -- it is a spectrum x strain occurrence matrix. All 15 cached
tasks are exactly the 15 that `nplinker_podp_screen.csv` covers, so the GNPS
download step is already done locally.
