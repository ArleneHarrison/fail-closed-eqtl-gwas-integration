# Downstream SuSiE-RSS / coloc.susie technical benchmark status

Date: 2026-09-01

Scope: technical audit only. No posterior, locus, gene, tissue, CAD, or HF result
from this smoke benchmark is used for biological interpretation or manuscript
claims.

## Input conformance audit

- Smoke input: frozen region `CAD_L01`, four GTEx heart/vascular tissues, CAD
  and HF outcomes (eight units).
- Gate status: 8/8 `READY`; gate rerun had already agreed.
- Ordered overlap: 2,538--2,891 variants per unit. All eight LD matrices had
  numerical rank 502 and the exact summary/LD identifier order agreed.
- Effect orientation: raw-source effects had been converted independently to
  the 1000 Genomes LD ALT allele. `alignment_status` retained whether source
  alleles required a flip; only `aligned` and `flipped` occurred.
- MAF: `eqtl_maf` and `gwas_maf` were present as separate finite columns in all
  units and were not substituted by the LD-panel MAF.
- Case fraction: the runner received only the frozen study-level value from
  `metadata.locked.v2.json`: CAD 0.116485962397774 and HF
  0.0816774392949421. HERMES row-level `N_case/N_total` varied by variant, as
  expected, and was not used to redefine the frozen study-level case fraction.
  The scalar `N` supplied to SuSiE-RSS was the runner's predeclared median rule:
  CAD 296,525 and HF 1,631,770--1,631,855 across these smoke units.

## Benchmark of the pre-correction runner

Environment: R 4.5.1, coloc 5.2.3, susieR 0.14.2. Runner SHA-256:
`c6474f17039d645d3cf243b55c66cf56a773cefcde8ff6a55a0240ee48d9ce31`.
Each unit had a 1,800-second external boundary.

| Unit class | Units | Outcome | Result | Elapsed per completed run | Peak RSS |
|---|---:|---|---|---|---|
| Four tissues | 4 | CAD | 4/4 completed; one returned no comparable credible-set pair | 43.7--50.0 s | 882,776--902,316 KB |
| Four tissues | 4 | HF | 4/4 reached the 1,800-s boundary; no completed result | 1,800.1 s | exact GNU-time peak unavailable after external termination; live observed lower bounds were approximately 854,108--1,002,064 KB |

The two CAD repeats were byte-identical for every retained model result in all
four units. No second completed HF repeat existed, so HF byte reproducibility
was not claimed. Stable benchmark codes were `OK`,
`OK_NO_COMPARABLE_CREDIBLE_SETS`, and `E_MODEL_TIMEOUT`.

## Root cause and pre-primary implementation correction

The pre-correction runner passed `max_iter=1000` but did not override coloc
5.2.3's `runsusie(..., repeat_until_convergence=TRUE)` default. When a SuSiE
fit did not converge, coloc multiplied the iteration count by 100 and reran it.
Thus the written 1,000-iteration contract was not actually enforced, explaining
the HF runtime. This was an implementation-conformance defect, not a scientific
result.

Before any formal primary model dispatch, the runner was corrected to assert
coloc 5.2.3 and susieR 0.14.2 and to pass every result-affecting default on the
active, `N`-supplied API path explicitly, without selecting values from model
outputs: `L=10`, `seed=20260901`, `max_iter=1000`,
`repeat_until_convergence=FALSE`, `estimate_residual_variance=FALSE`,
`tol=0.001`, `coverage=0.95`, `min_abs_corr=0.5`,
`scaled_prior_variance=0.2`, `check_prior=TRUE`, and `z_ld_weight=0`.
The coloc settings were explicitly fixed at `p1=1e-4`, `p2=1e-4`,
`p12=5e-6`, `overlap.min=0.5`, and `trim_by_posterior=TRUE`.
`prior_variance=50` is retained only as an inactive version-specific API audit
field because susieR 0.14.2 does not consume it as the model prior when `N` is
supplied; the effective parameter is the separately locked
`scaled_prior_variance=0.2`. `back_calculate_lbf=FALSE` is passed explicitly
but marked as accepted and unreferenced by the coloc 5.2.3 function body. The
scientific design and endpoints were unchanged. Dynamic server tests proved
the API route, effective scaled-prior transfer, and 1,000-iteration ceiling.

Final corrected runner SHA-256:
`f46fe8e94c9371b98dadac81ce025769f7c8a2888a3470f1d957c31d5e84680e`.
The machine-readable implementation lock SHA-256 is
`a771251b9e9db3e56210cfd8229393b6bba1110d696431eef82adafd0f82f345`,
bound to protocol canonical SHA-256
`ec619b392221037bf63628919b22bba0b45d9f1c0ccd720b9ec52531e0bf0ba6`.
The lock states explicitly that prior engineering smoke model runs had been
seen, that no formal primary model result had been seen before locking, and
that smoke outputs are excluded from the formal primary result set. Parameter
choices were made for implementation conformance, not from posterior values.

## Audit artifacts

- `downstream_benchmark_smoke_v3/benchmark_summary.tsv`
- `downstream_benchmark_smoke_v3/benchmark_summary.json`
- `downstream_benchmark_smoke_v3/units/*/input_audit.json`
- `downstream_benchmark_smoke_v3/units/*/unit_benchmark.json`
- `scripts/external_validation_20260901/benchmark_coloc_susie_smoke.py`
- `scripts/external_validation_20260901/run_coloc_susie_rank_diagnostic.R`
- `scripts/external_validation_20260901/config/model_implementation.lock.json`
- `scripts/external_validation_20260901/config/model_implementation.freeze.json`
- `scripts/external_validation_20260901/dispatch_ready_models.py`

The smoke artifacts remain a technical benchmark only and must not be merged
with formal primary model outputs.

## Superseded v6 engineering probe

The earlier v6 probe used the preceding, incomplete implementation lock. It is
retained only as engineering evidence and excluded from the formal primary
analysis. It observed:

- CAD: `COMPLETED / OK`, 48.06 s, peak RSS 900,860 KB.
- HF: `FAILED / E_SUSIE_NONCONVERGENCE`, 893.89 s, peak RSS 1,081,396 KB. It
  stopped naturally after exactly 1,000 iterations; no x100 retry occurred.
- Resume check: 1.83 s. The completed CAD result was skipped only after binding
  and output-hash verification; the prior HF failure was retained without an
  implicit retry.
- Refusal check: a synthetic non-READY gate produced `E_GATE_NOT_READY` without
  invoking R or producing an RDS file.

It must not be described as the final locked probe or merged with formal model
outputs. A new CAD/HF bounded probe under the final runner and implementation
lock is reported separately below.

## Final result-before-run conformance probe and dispatcher decision

The final runner, lock, freeze and dispatcher were deployed read-only before
this probe. Gate JSON was regenerated with the unambiguous
`gate_ready_for_dispatch` marker while the scientific `model_eligibility`
field remained `null`. One CAD and one HF unit were run with two workers and a
3,600-second external boundary:

- CAD: `COMPLETED / OK`, 45.11 s, peak RSS 901,292 KB.
- HF: `FAILED / E_SUSIE_NONCONVERGENCE`, 885.14 s, peak RSS 1,071,120 KB. The
  stable reason records non-convergence after exactly 1,000 iterations; no x100
  retry occurred.
- Safe resume: 0.61 s. CAD was `SKIPPED_VERIFIED_COMPLETED` only after binding
  and result-hash verification; HF was
  `SKIPPED_PREVIOUS_FAILED_USE_RETRY_FAILED_TO_RERUN` and was not implicitly
  rerun.
- Final non-READY semantic probe: `REJECTED / E_GATE_NOT_READY`,
  `model_invoked=false`, and zero RDS files.

Verification completed before the decision: 25/25 local Python tests, server
implementation-freeze check, dynamic coloc/susieR API route and iteration
ceiling test, and minimal end-to-end R runner test all passed. No formal primary
model dispatch had started when these corrections, locks and probes were
completed.

Dispatcher decision: **GREEN** for formal execution. Twelve workers are a
conservative default on the observed server (112 logical CPUs and roughly
454 GB available memory): the measured per-job peak was about 1.1 GB. Disk is
the tighter resource because each READY unit materialises a CSV LD binding.
With approximately 105 GB free on `/tmp` at the decision point, monitor free
space during the run and pause new dispatch if the projected locked-input CSVs
would exhaust the volume. Do not increase `L`, the iteration ceiling, or the
coloc priors in response to `E_SUSIE_NONCONVERGENCE`; it is a declared model
failure endpoint.

The non-posterior final probe audit is retained under
`formal_v6_bounded_probe_final_defaults_v2/`; it contains dispatch status,
bindings, gate JSON, policy/status CSVs, timing, memory and logs, but no copied
posterior RDS or biological interpretation.
