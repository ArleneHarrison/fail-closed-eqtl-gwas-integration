# External cardiovascular validation data preparation audit

Generated: 2026-09-01T17:40:13Z

## Outcome

The data layer is technically ready for held-out heart/vascular eQTL--CAD/HF workflow validation. Four GTEx v8 cardiovascular tissues, the original CAD reference, an additional 2017 CAD resource, HERMES 2.0 European-ancestry HF, and a 503-person 1000 Genomes EUR LD panel have been inventoried and checked. No whole-genome 1000 Genomes VCF was downloaded.

## Verified resources

- QTD000131 (aorta), QTD000136 (coronary artery), QTD000251 (atrial appendage), and QTD000256 (left ventricle): exact official byte sizes, complete gzip streams, readable official tabix indexes, 23 contigs each, and successful regional queries.
- HERMES 2.0 EUR: exact outer-archive size of 1,569,449,258 bytes. The nested archive contains four HF phenotype files with matching embedded MD5 values. The README declares b37/hg19 coordinates and 15 columns. Pheno1 (all clinical HF) was extracted as the primary outcome and passed gzip and tabix checks.
- GCST005195: the official 2017 CAD/UK Biobank gzip and README were downloaded. The stream is intact and contains 7,947,837 data rows and 13 columns. The GWAS Catalog API reports `snpCount=7,947,838`, a difference of one; because the downloaded bytes and gzip stream are intact, this is retained as a metadata-count note rather than treated as file corruption. Locked analysis metadata are 34,541 cases / 296,525 total, with case fraction `s=0.116485962397774`.
- Locked HERMES2 Pheno1 analysis metadata are 139,533 cases / 1,708,342 total, with case fraction `s=0.0816774392949421`.
- 1000 Genomes: the official population panel contains 503 EUR samples. One-byte range requests returned HTTP 206 for chromosomes 1--22, and the 22 official tabix indexes were cached with SHA-256 values. Full chromosome VCFs were not downloaded.
- Coordinate handling: eQTL inputs are GRCh38 while both GWAS inputs and 1000 Genomes LD are GRCh37. Both official directional chain files are present. Regional GRCh37 windows must be mapped to GRCh38 before eQTL tabix queries, and mapped boundaries/variants must pass unique exact reverse-roundtrip checks; using the same numerical window on both builds is prohibited.
- eQTL effective sample size: the prepared extractor derives per-row `effective_n=an/2` only after verifying that every retained `an` value is a positive even integer. The source `an` and derivation rule remain in the output and audit JSON.

## Scope of externality

These are resource- and tissue-level external technical validations, not guaranteed sample-independent replications. GTEx tissues can share donors. GCST005195 may overlap cohorts later represented in Aragam 2022. HERMES releases can share contributing cohorts. Manuscript wording should therefore avoid `independent replication` unless cohort overlap is quantified and excluded.

## Access and publication conditions

No controlled-access application or individual-level data application is needed for these public aggregate files. Nevertheless, study-specific reuse terms must be checked. The HERMES README asks journal users to follow consortium authorship/publication conditions; this requires author confirmation before submission. Raw third-party files should not be redistributed unless their terms permit it.

## Storage blocker

The server filesystem has global free space, but this user is at the 1 TiB hard quota. The newly downloaded HERMES and GCST005195 files are therefore physically cached under `/tmp` and exposed under the project by symlinks. They are usable now but are not durable across cache cleanup or reboot. A quota increase of at least 2 GiB is needed to retain only the outer source archives; 10 GiB of practical headroom is recommended for extracted outcomes, regional intermediates, logs, and reruns.

## Machine-readable evidence

- `audit/data_manifest.tsv` and `audit/data_manifest.json`
- `audit/verification_summary.tsv`
- `audit/sha256_manifest.tsv`
- `audit/hermes2_member_md5.tsv`
- `audit/eqtl_tabix_index_audit.tsv` and `audit/eqtl_tabix_query_audit.tsv`
- `audit/1000g_remote_reachability.tsv` and `audit/1000g_superpopulation_counts.tsv`
- `metadata/` contains headers, official README material, phenotype definitions, source metadata, and archive listings.
