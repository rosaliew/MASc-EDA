# Ground-truth data we do not have, and why

Generated from `ground_truth.csv` (115 links), `refetch_plan_all.csv` and `mibig_migration_report.csv`.

A `mapping_gap` means the file is on disk but nothing joins it to the link. A `file_gap` means the file is not in hand.

## Mapping gaps -- files exist, join missing

### G02 · gnps_component_bridge_not_run (automatable)

- **file kind**: genome
- **ground-truth rows affected**: 48
- **why**: Every one of these is a 'GNPS molecular family' link. They carry no MS2_URL, so the ms2_url_bridge cannot fire and no genome is attached. The genome FILES are on disk -- 297c364c alone holds 130 -- this is a join that was never built, not missing data. SCREENED 2026-08-27: NPLinker's own podp mode supplies this join for free. All 14 affected studies pass pre-flight (48/48 links).
- **next**: Run NPLinker in podp mode per study rather than hand-building a bridge; its gnps/file_mappings.tsv is a spectrum x strain matrix. Verified on task c22f44b14a: 27/27 real strain columns resolve to a genome accession. See nplinker_podp_screen.csv.

### G03 · no_ms2_url_on_family_link (automatable)

- **file kind**: ms2_run
- **ground-truth rows affected**: 48
- **why**: Same 48 family links. PoDP records a network component, not a run file, so there is no URL to resolve until the component is expanded to its member spectra. NPLinker's podp mode does that expansion as part of arrange_gnps().
- **next**: Same as G02 -- NPLinker podp mode, not a hand-built bridge.

## Not attempted yet

### G01 · spectrum_extraction_not_run (automatable)

- **file kind**: ms2_spectrum
- **ground-truth rows affected**: 115
- **why**: No per-link MGF has been cut yet. 66 rows carry an explicit ms2_scan and 49 carry a GNPS component; both are extractable from files already downloaded.
- **next**: Extract scans from the mzML/mzXML already on disk.

## File gaps -- not in hand

### G04 · massive_download_failed (manual)

- **file kind**: ms2_run
- **ground-truth rows affected**: 1
- **distinct items**: 1 — `20170919_NM6_glyc2.1.mzXML`
- **why**: URL is present and well-formed but the file was not retrieved. Only MS2 run in this state.
- **next**: Retry from the MassIVE FTP path by hand.

### G05 · mibig_absent_from_both (unavailable)

- **file kind**: bgc
- **ground-truth rows affected**: 1
- **distinct items**: 1 — `BGC0002077`
- **why**: Cited by PoDP but present in neither MIBiG 3.1 nor 4.0. The accession may never have been issued.
- **next**: Confirm the accession against MIBiG directly; if invalid, the link's BGC claim cannot be grounded.

### G06 · mibig_dropped_from_4.0 (manual)

- **file kind**: bgc
- **ground-truth rows affected**: 5
- **distinct items**: 4 — `BGC0000962|BGC0001476|BGC0001599|BGC0002061`
- **why**: Active with a structure in MIBiG 3.1, but dropped from the 4.0 bundle. The 3.1 record still exists upstream.
- **next**: Decide whether to pin these to their 3.1 record or accept the loss. Note BGC0000962 (barbamide) is NPOmix-covered.

### G07 · ncbi_returned_wrong_strain (manual)

- **file kind**: genome
- **ground-truth rows affected**: 0
- **strain-pairing rows affected**: 166 (strain-axis depth, not links)
- **distinct items**: 21 — `2515154124|2515154177|2516143022|2516493032|2516653042|2517287019|2517287023|2517434008|2517572194|2518285561|2524614807|2528311034|... (+9 more)`
- **why**: NCBI's assembly index ignores strain tokens, so a strain-name query silently returns an arbitrary genome of the same genus. verify_strain() caught these; a content hash found six accessions collapsing onto one file. Quarantined, never promoted.
- **next**: Re-resolve via BioSample [Strain] then elink to assembly -- confirmed working on CNS237, CNH643, CNS205.

### G08 · metagenome_scale_assembly (manual)

- **file kind**: genome
- **ground-truth rows affected**: 3
- **strain-pairing rows affected**: 2 (strain-axis depth, not links)
- **distinct items**: 2 — `JAAHTG010000000|JAAHTH010000000`
- **why**: Sequence is present and complete, but the assembly is metagenome-scale rather than an isolate. Held back by the scale gate so it cannot be mistaken for a single-strain genome.
- **next**: A judgement call, not a gap: antiSMASH will run on these. Promote deliberately if metagenome-derived BGCs are acceptable.

### G09 · no_sequence_after_refetch (manual)

- **file kind**: genome
- **ground-truth rows affected**: 2
- **strain-pairing rows affected**: 91 (strain-axis depth, not links)
- **distinct items**: 10 — `2515154126|2515154180|2561511034|2561511103|2563366517|2563366531|2563366533|AZHW00000000.1|GCF_000514575.1|GCF_000514975.1`
- **why**: Re-fetch still yielded a record with no ORIGIN block -- typically a WGS/CON master record that indexes contigs without carrying them.
- **next**: Fetch the component contigs, or the assembly FTP _genomic.gbff.gz for the corresponding GCA/GCF.

### G10 · no_resolution_route (manual)

- **file kind**: genome
- **ground-truth rows affected**: 0
- **strain-pairing rows affected**: 76 (strain-axis depth, not links)
- **distinct items**: 9 — `2515154126|2515154180|2561511034|2561511103|2563366517|2563366531|2563366533|GCF_000514575.1|GCF_000514975.1`
- **why**: Neither the accession, the BioSample, nor the strain name resolved to anything at NCBI.
- **next**: Check the source publication for a deposition, or ask the depositor. Some may be unreleased.

### G11 · jgi_img_only_no_ncbi_mirror (deferred)

- **file kind**: genome
- **ground-truth rows affected**: 0
- **strain-pairing rows affected**: 72 (strain-axis depth, not links)
- **distinct items**: 7 — `2515154126|2515154180|2561511034|2561511103|2563366517|2563366531|2563366533`
- **why**: IMG Taxon OIDs with no NCBI equivalent. NPOmix used the same JGI identifiers, so their supplement offers no shortcut. DEFERRED 2026-08-27: study 297c364c already has 130 of its 157 genomes on disk, so these 7 do not change rankability, and none is NPOmix-covered.
- **next**: DEFERRED. If revisited: img.jgi.doe.gov, search by Taxon OID, download 'Genbank (.gbk)' (not FASTA), verify strain against genome_label. Login required.
