# Raw-to-prepared pipeline blocking audit

Audit scope: `scripts/external_validation_20260901/prepare_external_validation_from_raw.py` and its unit tests, before the formal 24-region run.

Audit standard: only P0/P1 defects that can stop the run, change a scientific result, depart from the frozen contract, or make an input/output falsely reproducible are listed. The file was changing during review; the first blocking snapshot had SHA-256 `94afd960155ecf787a230bfe1502cafafe95e9035514d454ff501fef1f3123bc` and 869 lines. A second audit is required after the fixes are complete.

## P0 — annotation-only eQTL aliases are treated as discordant scientific records

The eQTL registry key is gene, chromosome, position, REF, and ALT, but duplicate equality is currently evaluated over the entire source row, including `rsid`. eQTL Catalogue legitimately emits repeated analytic rows with different rsID annotations. In the real first frozen window (`CAD_L01`), direct source queries found 962, 986, 843, and 760 duplicate keys in QTD000131, QTD000136, QTD000251, and QTD000256, respectively. Every examined duplicate differed only in `rsid`; all analytic fields were identical. The current main loop calls the full-row duplicate function, so the formal run stops in its first region with `E_DISCORDANT_DUPLICATE_KEY`.

Minimum fix and acceptance:

- Define an explicit analytic equality payload. `rsid` alone may differ and its complete alias set must be written to an audit table. Gene/molecular-trait identity, build coordinate, REF, ALT, beta, SE, P value, allele count, MAF, allele count/frequency support fields, test type, and all other analysis-relevant fields must remain identical.
- If any analysis-relevant field differs for the same analytic key, retain fail-closed `E_DISCORDANT_DUPLICATE_KEY`; never pick the first or last row.
- Use the analytic duplicate function in both candidate counting and the post-liftover map. Do not leave it as an unused helper.
- Add tests showing (i) two rows differing only in `rsid` collapse once and preserve both aliases, (ii) a beta/SE/AN/MAF or trait-identity difference stops, and (iii) the real first-region smoke reaches all eight units and emits nonzero alias-audit counts.

## P1 — the run can consume different files while recording the old manifest provenance

The original implementation read metadata from `data_manifest.json` but did not bind the files actually opened to the manifest byte counts and checksums. This affects the four eQTL files, their tabix indexes, CAD, the extracted HERMES phenotype and index, the HERMES outer archive/member relation, the EUR sample panel, and the two chain files. In particular, changing the sample panel to any other 503 unique EUR-labelled samples would change LD while passing the existing sample-count check; the resulting unit provenance did not identify the panel.

Minimum fix and acceptance:

- Before any region or result is evaluated, verify the exact path, byte count, and frozen checksum for every locally consumed asset and index. The four large eQTL files may be hashed once per full run; no whole-chromosome 1000 Genomes VCF download is required.
- For the extracted HERMES phenotype, verify the provider-supplied member MD5 (or an equivalently frozen member SHA-256) and record the computed SHA-256, member path, and parent-archive SHA-256. A placeholder string is not an input checksum.
- Verify the EUR panel and both liftOver chains against their manifest SHA-256 values and include their hashes in the run-level provenance.
- For remote 1000 Genomes data, retain range requests. Record the stable URL plus the deterministic decoded regional fingerprint already produced by the code; verify the cached `.tbi` checksum for every chromosome used. This is sufficient and does not require downloading whole VCFs.
- Make the range reader actually use the verified local cached index (for example, the `index_filename` argument supported by the VCF reader). Merely hashing a local index while the reader silently uses the remote sibling `.tbi` does not bind the retrieval to that verified index.
- Add a mutation test for each asset class: a one-byte or size change must stop before `prepared/` is populated. Include the verified manifest SHA-256 and preparation-code SHA-256 in the run summary.

## P1 — frozen protocol, metadata, and region selection are not fully bound at preparation time

The original preparation entry point checked only protocol ID plus a Boolean `frozen` field, checked only the metadata-lock state string, and accepted any table with 24 unique region IDs. Thus a modified protocol, a lock belonging to another protocol, or a different set of 24 regions could be prepared and reported as a passing v2 run.

Minimum fix and acceptance:

- Verify `protocol.v2.freeze.json`: file bytes, file SHA-256, canonical JSON SHA-256, protocol ID, frozen timestamp, and read-only state.
- Verify that `metadata.locked.v2.json` names the same protocol canonical hash and that cases + controls = N and the recorded case fraction equals cases/N for CAD and HF.
- Bind the region table to a pre-result audit that contains its own SHA-256, protocol canonical SHA-256, and both GWAS selection-input SHA-256 values.
- Independently validate 12 CAD plus 12 HF regions, ranks 1–12, the expected IDs, significant lead threshold, lead within its own declared chromosome/window, and exact `max(1, lead-500000)`/`lead+500000` boundaries. Reject unexpected columns/empty fields that control selection.
- Add tests in which one protocol byte, one metadata count, one region coordinate, and one selection-source hash are altered; each must stop before data preparation.

## P1 — CAD effect direction depends on an audit that was not executable input

The pipeline interprets CAD `a1` as the effect allele. The signed-sentinel audit exists, but the original run did not require or verify it. A missing/stale audit or wrong A1 interpretation could reverse CAD effects without a gate failure.

Minimum fix and acceptance:

- Freeze the expected sentinel identifiers, raw A1 alleles, and raw beta signs/values in a versioned lock or audit whose SHA-256 is bound to the run.
- During the existing single CAD scan, recover all locked sentinel records directly from the exact CAD file and compare identifier, A1, beta sign/value within a stated tolerance, and source SHA-256. All must pass before any prepared unit is written.
- Record the sentinel-audit SHA-256 and observed values in `cad_single_scan_audit.json`.
- Add pass, missing-sentinel, flipped-A1, and changed-beta tests.

## P1 — a mismatched tabix index can silently return the wrong analytic rows

The original tabix reader checked row width but trusted the index to return the requested interval. Because index hashes were not bound and returned coordinates were not checked, a wrong `.tbi` could select the wrong gene or variants while producing syntactically valid rows.

Minimum fix and acceptance:

- Verify each `.tbi` against the frozen manifest (including a frozen checksum for the HERMES phenotype index).
- For every fetched row, require normalized chromosome equality and a one-based position inside the requested inclusive interval. Fail on absent/ambiguous coordinate columns.
- Test a mocked fetch containing one outside-interval row and require `E_TABIX_INDEX_ROW_OUTSIDE_QUERY`.
- Also test an index that returns an in-range subset. Prevent this for the locally generated HERMES index by freezing its checksum before the full run; interval checks alone detect wrong rows but cannot detect silently omitted rows.

## P1 — LD cache validation is incomplete

On an ordered-variant hash cache hit, the code validated only the primary mean-imputed matrix. It did not require or validate the complete-case archive. A stale, missing, or corrupted complete-case sensitivity matrix could therefore be materialized under a valid primary cache key.

Minimum fix and acceptance:

- On every cache hit, require both archives; check both variant orders, dimensions, and matrix equality (including NaNs) against the freshly computed primary and complete-case matrices.
- Bind the cache to the EUR-panel hash and decoded regional 1000 Genomes fingerprint, or retain the current equality recomputation plus verify both matrices.
- Add tests for missing, wrong-order, and altered complete-case cache files.

## P1 — output reuse can mix stale and current runs

The output directory, LD cache, and unit directories are opened with `exist_ok=True`. An interrupted run or reuse of the same path can leave stale regions/units and a misleading `prepared_units` count; there is no resume binding to protocol, code, manifest, region, and input hashes.

Minimum fix and acceptance:

- Default to a new empty output directory and refuse a nonempty output path.
- If resume is supported, require a run-binding record containing protocol, metadata lock, code, manifest, region registry, panel, chain, and source hashes; validate every existing unit before reuse.
- On successful full completion, require exactly 24 region audit files, 24 gene rows, 192 unique attempt IDs, and 192 unit directories matching the Cartesian product. Write `run_summary.json` atomically only after these checks.

## Indel decision

Indels do not need to be removed solely because they are indels. They may remain when the source position has a unique positive-strand exact point round trip and the complete, non-complemented REF/ALT strings match the biallelic GRCh37 1000 Genomes record as an exact allele pair/set, with the effect allele explicitly resolvable to the target REF or ALT. The first-region HF units contain approximately 148–152 such indels each, so the rule must be tested and their counts reported. The manuscript must describe this as point-coordinate round-trip plus exact target-reference allele reconciliation, not as full-interval liftOver. A full-interval start/end check is an optional stronger sensitivity, not a prerequisite for the current primary run.

Required tests: one valid insertion/deletion, one coordinate-only match with different allele strings, one reverse-complement-only match, one multiallelic target, and one ambiguous/negative-strand mapping. Only the first should be retained.

## Tests present at the blocking snapshot

Only three helper tests were present and all passed. They covered full-row duplicate collapse, LD correlation ordering/complete-case shape, and zero-overlap dimensions. They did not exercise the raw schema, real alias duplicates, protocol/manifest binding, tabix range validation, sentinel direction, HF/CAD denominator handling, allele orientation, LD-cache reuse, no-gene 192-attempt registration, provenance, or an end-to-end first-region smoke. Passing those three tests was therefore not sufficient authorization for the formal full run.

## Final second-round disposition: GREEN for raw full preparation

The blocking findings above were re-audited against preparation-script SHA-256 `129918a8adbfebc6da08e09022ba8331dd19f504863a20d3381d29d5a9d9b7fc`. No remaining P0/P1 defect was found. The fixes are executable from `main`, and the preparation-specific test suite reports 10 passing tests.

The final frozen one-region smoke (`raw_prep_smoke_v6`) completed with status `PASS`: one gene, eight registered attempts, and eight prepared tissue-by-outcome units. Each unit contains exactly the seven required artifacts. Ordered-overlap counts were 2,538--2,580 for CAD and 2,838--2,891 for HF, matching the pre-fix smoke fingerprints. Summary hashes, ordered-variant identities, LD variant order and dimensions, finite/symmetric/unit-diagonal LD matrices, protocol/region/runtime-index/source hashes, and the actual preparation-code hash were verified read-only. All 27 consumed source/index assets were present in the verification registry with positive byte counts, SHA-256 values, and passing status.

The remote 1000 Genomes reader uses the verified cached index. HTSlib records whose intervals overlap a query bin but whose POS anchor lies outside the POS-defined window are audited and skipped; contig mismatches still stop. One reciprocal indel representation produced a genuinely discordant unordered reference key in this region. The complete key was quarantined before overlap rather than selected by record order; it matched none of the eight analytic units, and the final analysis view reports zero discordant duplicate keys.

Annotation-only eQTL duplicates were collapsed only when their analytic payloads were exact, with 3,551 alias rows retained in the audit. QTL `an` is required to be finite, positive, and even before `an/2`; CAD and HF effect alleles are explicitly resolved to the target REF/ALT and both beta directions are transformed to target ALT. CAD uses the source per-row total N and locked 34,541 cases; HF preserves source per-row `N_total` and `N_case`, while the separate locked study-level case fraction remains an audit field rather than replacing row-specific denominators.

Two independent gate materializations each returned 8/8 `READY` with code `OK`, and the corresponding per-unit gate hashes were identical. The rerun comparison reported 8/8 byte and semantic agreement.

This GREEN authorizes starting the 24-region raw full preparation in a new empty directory on `/tmp`, where the recorded space preflight passed. It is not a downstream-analysis GREEN. The completed full output must still satisfy the exact 24-region/24-gene/192-attempt/192-unit Cartesian checks, then pass two independent gates and their repeat/semantic comparison before downstream work is released.
