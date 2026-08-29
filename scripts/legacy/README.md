# Retired PoDP scripts

Superseded by `scripts/build_ground_truth.py`, which produces
`data/PoDP/ground_truth.csv` (115 links, 63 columns).

These two wrote the flat manifests now archived in `data/PoDP/legacy/`:

* `download_verified_bgc_ms2.py`        -> `verified_bgc_ms2_manifest.json`
* `download_verified_bgc_ms2_genomes.py` -> `verified_bgc_ms2_genome_manifest.json`

Verified 2026-08-27: all 67 links in both manifests are present in
`ground_truth.csv`, matched on project + MS2 scan + run URL. Nothing was lost.

They also relied on `is_verified_bgc_ms2_link()`, the selection gate that was
replaced -- it encoded link *type* rather than evidence strength, discarding all
48 `GNPS molecular family` links while admitting links whose only evidence was
"Evidence as indicated in MIBiG".

`scripts/podp_online_file_locator.py` was deliberately NOT retired: it stays in
`scripts/` because `build_ground_truth.py` imports its genome-bridge helpers.

To run one of these, point it at `../data/PoDP/legacy/`.
