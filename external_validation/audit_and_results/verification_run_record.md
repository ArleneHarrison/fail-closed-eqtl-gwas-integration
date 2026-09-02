# External technical-evaluation verification record

Run date: 2026-09-01 UTC. Scope: implementation and workflow verification;
no posterior, gene, locus, tissue, or clinical interpretation.

## Local Python collection

- Environment: Python 3.12.10; NumPy 2.3.5; pandas 2.3.3; SciPy 1.17.1;
  pytest 9.1.1.
- Command scope: `scripts/external_validation_20260901/tests`.
- Result: **25 passed in 0.32 s**.
- This collection is distinct from the 93-test core implementation suite and
  from the two R integration tests below.

## Version-pinned server checks

- `test_coloc_api_ceiling.R`: **PASS** for the explicit coloc 5.2.3 / susieR
  0.14.2 API route and iteration ceiling.
- `test_runner_minimal.R`: **PASS** for the minimal locked end-to-end runner.
- Implementation-freeze check: **PASS**; protocol canonical SHA-256
  `ec619b392221037bf63628919b22bba0b45d9f1c0ccd720b9ec52531e0bf0ba6`.
- The eight-unit serial/parallel, one-BLAS-thread gate smoke was recorded
  separately and agreed byte-for-byte in 8/8 units.
- The final bounded probe recorded CAD `COMPLETED/OK`, HF
  `FAILED/E_SUSIE_NONCONVERGENCE` at exactly 1,000 iterations, safe resume,
  and non-READY refusal without model invocation. These probe outputs are
  engineering evidence and are excluded from the formal primary result set.

## Formal pre-model workflow evidence

- Raw preparation: **PASS**, 24 regions, 24 region-level gene selections
  representing 22 unique gene identifiers, and 192/192 prepared units.
- The 389,322 summary rows (CAD 181,656; HF 207,666) are unit-level rows in
  the 192 factorial summaries and include repeated source content across
  attempts; they are not independent variants or participants.
- Independent integrity audit: **GREEN**, no failure; 18 ambiguous reference
  keys quarantined and none intersected an analysis request.
- Formal gates v1/v2: **PASS**, 192/192 attempts each, eight workers per pass,
  one BLAS thread per worker, identical script/protocol/metadata/attempt hashes.
- Exact rerun comparison: **192/192 PASS**, byte-identical and
  semantic-identical.
- Complete-case LD sensitivity: **PASS**, 192/192 units; descriptive only and
  unable to change primary gate state or dispatch.

## Deployed implementation hashes

- Gate/summary runner: `7c49fa48fe8c80c4510497b7f63bbe0ee406722b8398bcaa6f8fe632959432fb`
- Model dispatcher: `e53975ce5ef4ac088aa92f70b04f7a467466e0eb1ea92434af3bcfd1f8d16d79`
- R model runner: `f46fe8e94c9371b98dadac81ce025769f7c8a2888a3470f1d957c31d5e84680e`
- Implementation lock: `a771251b9e9db3e56210cfd8229393b6bba1110d696431eef82adafd0f82f345`
- Implementation freeze: `894930f9434689a8b745bd89d5d628a2b1b8d48b6d8e5cbaee4df11d9ca979ad`
- Complete-case sensitivity runner: `e62669bb1ec1ff19852f1c274423cb342d0e0eeee8ac54cd4a863f9015666955`

Formal model-dispatch terminal counts are intentionally not entered here until
the complete 192-unit dispatch summary and its hashes exist.
