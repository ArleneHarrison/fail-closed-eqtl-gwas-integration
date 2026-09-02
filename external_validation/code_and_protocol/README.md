# Prospectively frozen cardiovascular external technical validation

This directory implements amended frozen protocol v2 (`protocol.v2.json`). The
original v1 remains read-only in the same directory, and
`amendment_ledger.tsv` proves that only pre-run metadata locks changed. It is an
analysis component, not a result package. No biological claim or posterior is
created by the planning and gate stages.
`implementation_correction_ledger.tsv` separately records result-before-run
implementation conformance corrections and explicitly excludes all earlier
engineering probe outputs from the formal primary analysis.

## Scope and honest externality label

Use the exact phrase **prospectively frozen multi-tissue external technical
validation**. Do not call it a fully independent replication: the four GTEx
tissues appeared in earlier exploratory scripts, and the older and newer CAD
GWAS generations may contain overlapping participants. The primary factorial
uses external CAD `GCST005195` and `HERMES2_EUR_2025`; development CAD
`GCST90132314` is not a primary outcome and may be shown only as a labelled
development/benchmark comparator.

## Frozen design

- Disease regions: 12 from external CAD and 12 from HERMES2 EUR. Within each
  GWAS, variants with P <= 5e-8 are sorted by P, chromosome, position and ID;
  a 1-Mb distance-pruned greedy rule selects 12 anchors, each expanded by
  +/-500 kb. This stage can read only GWAS chromosome, position, P and ID.
- Gene: for each frozen region, candidate genes must have at least 20 rows with
  the required fields present in every one of the four tissues. The shared
  candidates are sorted lexicographically by `gene_id`; the first is selected.
  Effect sizes and standard errors are checked only for presence, never value,
  direction or rank. No common candidate gives `E_NO_ELIGIBLE_GENE`.
- Attempts: 24 regions x 4 tissues x 2 outcomes = 192. Regions remain anchored
  to their source outcome but are evaluated against both outcomes.
- Primary endpoints: gate state and error codes, ordered variant overlap, LD
  dimension/rank/condition number, and deterministic rerun agreement.
- Only `READY` attempts can reach SuSiE-RSS/`coloc.susie`. A stopped attempt
  cannot contain an RDS/posterior artifact.
- Gate JSON uses `gate_ready_for_dispatch` solely as a boolean derived from
  `status == READY`; the distinct scientific field `model_eligibility` remains
  `null` and must not be inferred from this dispatch marker.
- LD rank is diagnostic, not a universal gate failure: with 503 EUR reference
  samples, a region containing more variants can be structurally rank deficient.
  Finite values, symmetry, unit diagonal, PSD and exact variant order remain hard
  contracts. The pinned R runner allows rank-deficient PSD matrices without
  conditioning or jitter and records rank and infinite condition number.

This gene rule is executable and protects against outcome-driven selection. Its
cost is potentially high, informative attrition because requiring a common gene
across four tissues is deliberately strict. The attrition must be reported as a
technical endpoint, not repaired by selecting a different gene after seeing the
gate or model results.

## Input handoff contract

Each gate-ready-to-check unit is prepared at:

```text
prepared/<unit_id>/
  summary.csv
  ld.npz
  provenance.json
  overlap.json
```

`summary.csv` must contain the existing fail-closed gate fields, including an
ordered `target_variant`; `ld.npz` must contain `ld` and `variant_ids`; and each
provenance record must contain URL, build/version, access terms, retrieval date,
SHA-256, byte count and structural check. Inputs are never imputed or repaired.
`overlap.json` must record pre-intersection eQTL/GWAS/LD counts, ordered overlap
count and the SHA-256 of newline-delimited ordered variant IDs; disagreement with
the summary or LD archive gives `E_OVERLAP_AUDIT_INVALID`.
`external_gate.py` is an exact pinned snapshot of the project gate for servers
whose older code deployment predates that module; the environment audit records
whether the project module or snapshot was loaded.

"No imputation" applies to summary fields, alleles, metadata, identifiers and
variants. It does not hide LD genotype handling: the frozen primary LD policy
permits per-variant mean-dosage filling only after a <=2% missingness check and
requires complete-case-sample LD as a sensitivity analysis. Both matrices and
missingness traces must be retained.

## Commands

All paths below are examples and remain relative/configurable.

```bash
python run_external_validation.py validate-protocol \
  --protocol config/protocol.v2.json

python run_external_validation.py select-loci \
  --protocol config/protocol.v2.json --metadata-lock config/metadata.locked.v2.json \
  --cad data/GCST005195.gz \
  --hf data/HERMES2_EUR.gz --cad-separator '\s+' --hf-separator $'\t' \
  --cad-columns '{"chromosome":"chr","position":"bp","p_value":"pval","variant_id":"uniqid"}' \
  --hf-columns '{"chromosome":"CHR","position":"POS","p_value":"P","variant_id":"SNP"}' \
  --out work/regions.tsv

python run_external_validation.py select-genes \
  --protocol config/protocol.v2.json --regions work/regions.tsv \
  --tissue-tables '{"QTD000131":"work/qtd131_candidates.tsv","QTD000136":"work/qtd136_candidates.tsv","QTD000251":"work/qtd251_candidates.tsv","QTD000256":"work/qtd256_candidates.tsv"}' \
  --out work/genes.tsv

python run_external_validation.py plan --protocol config/protocol.v2.json \
  --regions work/regions.tsv --genes work/genes.tsv --out work/attempts.tsv

python run_external_validation.py gate --protocol config/protocol.v2.json \
  --metadata-lock config/metadata.locked.v2.json \
  --attempts work/attempts.tsv --prepared work/prepared \
  --out outputs/gates

python run_external_validation.py dispatch --protocol config/protocol.v2.json \
  --metadata-lock config/metadata.locked.v2.json \
  --attempts work/attempts.tsv --prepared work/prepared --gates outputs/gates \
  --out outputs/models

# Formal READY-only model execution uses the implementation-locked parallel
# dispatcher. It refuses non-READY units, publishes unit status atomically, and
# safely resumes only when the complete binding and output hashes still agree.
# The lock makes all result-affecting coloc 5.2.3/susieR 0.14.2 defaults on the
# active N-supplied path explicit: L, iteration ceiling, convergence policy,
# residual-variance policy, tol, coverage, minimum correlation, scaled prior
# variance, prior checking, z/LD blending, priors, overlap threshold and
# posterior trimming. Version-specific inactive API formals remain audit fields.
python freeze_model_implementation.py check \
  --protocol config/protocol.v2.json \
  --runner run_coloc_susie_rank_diagnostic.R \
  --lock config/model_implementation.lock.json \
  --freeze config/model_implementation.freeze.json

python dispatch_ready_models.py \
  --attempts work/attempts.tsv --prepared work/prepared --gates outputs/gates \
  --metadata-lock config/metadata.locked.v2.json \
  --protocol config/protocol.v2.json \
  --runner run_coloc_susie_rank_diagnostic.R \
  --implementation-lock config/model_implementation.lock.json \
  --implementation-freeze config/model_implementation.freeze.json \
  --out outputs/models_locked --jobs 12 --timeout-seconds 3600

# Run the gate a second time into outputs/gates_rerun, then compare all units.
python run_external_validation.py compare-reruns --attempts work/attempts.tsv \
  --left outputs/gates --right outputs/gates_rerun \
  --out outputs/rerun_agreement.tsv

# Read-only secondary sensitivity after the frozen primary gates_v1 exists.
# This requires exactly 192 attempts by default and never changes a primary
# gate state or the READY-only model-dispatch set.
python analyze_complete_case_ld_sensitivity.py \
  --attempts work/attempts.tsv --prepared work/prepared \
  --gates outputs/gates_v1 --protocol config/protocol.v2.json \
  --out-dir outputs/complete_case_ld_sensitivity --jobs 8
```

The complete-case report verifies identical primary/complete-case variant IDs
and order, then records complete-case sample count and LD integrity, rank and
condition diagnostics under the already frozen diagnostic-rank policy. Relative
Frobenius distance is `||R_complete-R_primary||_F / ||R_primary||_F`; the max,
median and 95th percentile use absolute upper-triangle off-diagonal changes in
correlation. No difference threshold is selected. Zero-overlap attempts are
retained as explicit `N_A_ZERO_OVERLAP`, and a non-computable complete-case
matrix is reported as a failure without repair. These secondary descriptions
cannot modify primary gate eligibility or downstream dispatch.
Independent unit workers can be selected with `--jobs`; executor mapping keeps
the frozen attempt-registry order, and production runs pin each worker's BLAS
thread count to one.

Pin exact source hashes after retrieval with `verify-sources`. Do not start the
primary run while source metadata, genome build or authoritative study-level
case fraction is unresolved. Protocol v2 records the prospectively resolved
GCST005195 ancestry, build and study-level case fraction in
`metadata.locked.v2.json`; protocol v1 is retained only as the immutable
pre-amendment record.
The CAD example uses the verified native `CAD_UKBIOBANK.gz` whitespace-delimited
header (`uniqid chr bp ... pval ...`). HERMES2 columns must be inspected and
locked rather than assumed from this example; if its delimiter differs, perform
lossless normalization to the four declared selection columns and hash both the
source and normalized file.

## Tests

```bash
python -m pytest -q tests
```
