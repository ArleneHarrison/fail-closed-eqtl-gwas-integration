# Pre-result methods recheck

> Historical checkpoint: this document reviewed the earlier 20-test external
> Python collection. Five complete-case/dispatch regressions were added before
> final packaging; the current authoritative record is 25/25 in
> `results/verification_run_record.md`. The findings below are retained as an
> audit trail and were resolved before formal model dispatch.

## Scope

- Manuscript reviewed: `paper/submission/interdisciplinary_sciences_20260901/latex_source/main_external_validation_draft.tex`
- Initial reviewed snapshot SHA-256: `9823ff8087a354975f84f7da90af3a491763f8432082e68018a9421fcfd87720`
- Current manuscript SHA-256 at final read-only verification: `f95e09fdc5bb1bcf784d69911eac7140be14dbbf698c62e47ebd7c1153225ee7`. The manuscript changed concurrently during this review; the cited methods passages were rechecked against the current file and all findings below still apply.
- Review type: second-round, read-only methods/code/audit consistency check.
- Result placeholder: retained. No primary result, posterior or biological interpretation was inserted.
- Focus: external preprocessing, HERMES denominators, BLAS threading, SuSiE/coloc implementation, and the reported 20 external tests.

## Overall disposition

No remaining **non-result P0** was identified in the reviewed snapshot. The new methods text is substantially concordant with the frozen files and deployed code. Seven **P1** items should nevertheless be resolved before the real result block is filled, mainly because one promised sensitivity analysis is not yet executed, some result-determining SuSiE defaults remain implicit, and the current supplementary staging directory does not contain the final implementation-lock/dispatcher materials.

## Facts independently reconciled

### External preprocessing

The manuscript lines 43--49 agree with `prepare_external_validation_from_raw.py` SHA-256 `129918a8adbfebc6da08e09022ba8331dd19f504863a20d3381d29d5a9d9b7fc`, which is identical locally and in the deployed server formal directory.

- Server environment confirms Python 3.13.12, NumPy 2.4.6, pandas 2.3.3, pysam 0.24.0 and pyliftover 0.4.1.
- `an/2` is used only after finite, positive, even-allele-number validation.
- `eqtl_maf`, `gwas_maf` and `ld_maf` are constructed separately.
- HERMES `N_case` and `N_total` are retained per variant; `gwas_n` is row-level `N_total`, `gwas_cases` is row-level `N_case`, and the locked study-level case fraction is a separate field.
- CAD uses the exact file denominator and fixed 34,541 case count.
- Primary LD excludes variants with more than 2% missing genotype dosage, mean-fills only retained variants, and writes a separate complete-case matrix.
- Overlap counts and eQTL/GWAS-denominated fractions are explicitly defined in manuscript line 47.

### HERMES denominator statement

Manuscript line 32 is consistent with `metadata.locked.v2.json` and the raw-preparation mapping:

- maximum aggregate cases: 139,533;
- maximum aggregate controls: 1,568,809;
- maximum aggregate total: 1,708,342;
- locked study-level case fraction: 0.0816774392949421;
- downstream scalar `N`: median of retained row-level `N_total`, not `N_max`.

The manuscript correctly distinguishes the model case fraction from per-variant denominator fields.

### BLAS thread control

The deployed and local `dispatch_parallel_gates.py` files share SHA-256 `be8dd6d156bdd566d72f8706d56558f5c4e5340ddb890bb0a1afb278ee405464`. The wrapper sets `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMEXPR_NUM_THREADS` to one before each gate subprocess.

The server smoke audit records `per_worker_blas_threads: 1`; serial versus parallel, single-thread comparisons passed byte-for-byte in 8/8 units. Thus manuscript line 49 is supported for the smoke. The final past-tense claim about “both formal passes” still needs the two formal dispatch audits; see P1-5.

### Locked SuSiE/coloc implementation

The manuscript line 51 agrees with both the runner and lock:

- R 4.5.1;
- coloc 5.2.3;
- susieR 0.14.2;
- $L=10$;
- seed 20260901;
- `max_iter=1000`;
- `repeat_until_convergence=FALSE`;
- `estimate_residual_variance=FALSE`;
- $p_1=10^{-4}$, $p_2=10^{-4}$ and $p_{12}=5\times10^{-6}$;
- median retained-row eQTL and GWAS sample sizes;
- nonconvergence mapped to `E_SUSIE_NONCONVERGENCE` without a larger-iteration retry.

Runner SHA-256 `c874950dfdacded856c598593bd4d5b7e1119ac4cdd9361cffce00ed6efab2c3` and implementation-lock SHA-256 `e58b2b8532edb4adfa62ce68f25bd5f34b71fc90dd1110f9ea22be7346a86904` match manuscript line 51, the local freeze, and the deployed server files.

Two dynamic R checks were rerun during this audit in the pinned server environment: the coloc `repeat_until_convergence` ceiling route passed, and the minimal locked runner passed.

### Twenty external tests

The current Python collection contains exactly 20 tests and reran as `20 passed`. It comprises ten protocol/selection/gate/runner-source tests and ten raw-preparation/integrity tests. This supports a “20 Python tests passed” claim, but it does not include the two dynamic R tests or the serial/parallel integration smoke. Those three evidence classes must be named separately; see P1-3.

## Remaining non-result P1 issues

### P1-1 — “Complete-case sensitivity analysis” is only materialized, not analyzed

**Manuscript locations:** lines 38 and 47.

**Code evidence:** `prepare_external_validation_from_raw.py` writes `ld_complete_case.npz` and a complete-case sample count. No gate, dispatcher, model runner, summary builder or figure script consumes `ld_complete_case.npz`; a repository-wide code search finds it only in preparation and preparation tests.

**Why this matters:** the text and protocol call this a “required sensitivity analysis,” but the implementation currently provides only a sensitivity artifact. No comparison metric, pass/fail rule, alternative gate diagnostic or downstream rerun is defined. Calling an unused matrix an analysis is inaccurate.

**Required resolution before result fill:** choose and lock one of the following without consulting primary outcomes.

1. **Diagnostic sensitivity:** run the same LD integrity/rank diagnostics on the complete-case matrix and prespecify comparisons with primary LD (complete-case sample count, Frobenius distance, maximum/median absolute $\Delta r$, rank and gate-diagnostic concordance). State explicitly that it cannot alter primary gate state or dispatch.
2. **Model sensitivity:** in addition to the above, rerun the downstream model on valid complete-case LD and prespecify comparison outputs and nonconvergence handling.
3. **Artifact only:** if neither analysis will be run, replace “requires ... as a sensitivity analysis” with “materializes and hashes a complete-case LD matrix as a prespecified sensitivity artifact; it is not used for primary gating or model dispatch.” This is the least expansive correction but must also be reconciled with protocol wording.

Do not decide among these options after inspecting formal primary model outputs.

### P1-2 — Result-determining SuSiE defaults remain implicit

**Manuscript location:** line 51.

**Code evidence:** the runner explicitly locks $L$, iteration ceiling, repeat behavior and residual-variance policy, but it does not pass convergence tolerance, credible-set coverage or credible-set purity. In susieR 0.14.2 the installed defaults are `tol=0.001`, `coverage=0.95`, `min_abs_corr=0.5`, `prior_variance=50`, `check_prior=TRUE` and `z_ld_weight=0`. `coloc.susie` also retains `back_calculate_lbf=FALSE` by default.

**Why this matters:** tolerance determines `E_SUSIE_NONCONVERGENCE`; coverage/purity determine whether comparable credible-set pairs exist. These can therefore change stable model terminal states, not merely posterior formatting.

**Required resolution:** before formal primary dispatch, either add these values explicitly to the runner and implementation lock (with new hashes and a documented pre-result implementation correction) or state in line 51 that all unlisted arguments used the version-pinned coloc 5.2.3/susieR 0.14.2 defaults and list at least tolerance, coverage and `min_abs_corr`. The explicit-lock option is stronger. Add corresponding assertions to the R minimal test.

Also add one clarifying sentence that the $L=10$ and coloc prior values were installed defaults made explicit in the post-smoke, pre-primary implementation lock; they were not fields in protocol v2 and were not selected from posterior values.

### P1-3 — “20/20 external preparation and dispatch tests” conflates three evidence classes

**Manuscript locations:** line 95, line 155 and Table 2 line 175.

**Observed evidence:**

- 20 Python tests pass;
- two separate dynamic R tests pass in the pinned server environment;
- an 8-unit single-thread serial/parallel gate smoke passes byte-for-byte;
- a separately documented two-unit bounded model probe and non-READY refusal probe exist but are not members of the 20-test pytest collection.

**Problem:** “external-preparation and dispatch suite passed 20/20” can be read as 20 dynamic end-to-end dispatch tests. Several of the 20 are source/lock assertions, and neither R test is counted in 20.

**Exact non-result methods replacement for the last part of line 95:**

```tex
The external Python suite comprised 20 protocol, selection, preparation, integrity, gate, summary, and implementation-lock tests. Two additional version-pinned R tests exercised the coloc iteration-ceiling route and a minimal end-to-end runner input. A separate eight-unit integration smoke compared serial and parallel single-thread gate execution. Test outputs, environments, and script hashes were retained separately.
```

After the final packaging rerun, report pass counts for each class separately in Results/Table 2 and include the exact environment of the 20-test run. Do not imply that the old 93-test environment also applied to the external 20 unless they were actually rerun together in that environment.

### P1-4 — `model_eligibility` and `model_eligible` have conflicting semantics

**Manuscript locations:** lines 38, 49 and 54.

**Code evidence:** `external_gate.py` emits `model_eligibility: null`, correctly matching line 54. `run_external_validation.py` then adds `model_eligible = (status == READY)`, and `dispatch_ready_models.py` refuses a unit unless both `status == READY` and `model_eligible is True`.

**Why this matters:** the two near-identical fields imply opposite things. In code, `model_eligible` is only a READY-derived dispatch marker, not a separately adjudicated scientific/model-eligibility result. A reviewer or future caller could treat it as evidence that all model assumptions passed.

**Required resolution:** rename the derived field to `ready_for_dispatch_check` or `gate_ready_for_dispatch`, update the dispatcher/tests, and retain `model_eligibility=null`. If a code rename is no longer feasible before the formal run, add an explicit sentence at line 49 and to the machine-readable schema:

```tex
The legacy Boolean field \texttt{model\_eligible} was a dispatch-control alias for \texttt{status==READY}; it did not populate the separate scientific field \texttt{model\_eligibility}, which remained null.
```

The rename is preferable because the manuscript should not need to explain an avoidable semantic collision.

### P1-5 — Past-tense formal-pass and registry claims must be bound to final audit files

**Manuscript location:** line 49.

**Problem:** code and the eight-unit smoke support single-thread dispatch, but at this pre-result checkpoint the two full formal gate-pass audit files were not present in the locally archived result package. “Gate workers used one BLAS thread each in both formal passes” must not be accepted solely because the wrapper is configured that way.

**Required release check:** for each formal pass require:

- 192 attempt coverage;
- `per_worker_blas_threads=1`;
- the same deployed gate-script hash;
- protocol, metadata-lock and attempts hashes;
- byte-copy merge policy;
- PASS status;
- exact rerun agreement after canonical comparison.

Keep the manuscript past tense only after these records exist and have been copied to durable storage. Otherwise use planned/future tense until execution.

### P1-6 — The HERMES row-audit wording omits invalid-row exclusion and the median-$N$ limitation

**Manuscript locations:** lines 32 and 51.

**Code evidence:** finite HERMES rows are retained; rows with missing/non-finite numeric fields or invalid allele pairs are counted and excluded before the analytic map. Gate checks later require positive `gwas_n` and `gwas_cases <= gwas_n`. The scalar downstream $N$ is the median of retained row-level `N_total` values.

**Problem:** “retained and audited row by row” could be read as retaining every deposited row. The median-$N$ choice is transparent but still an approximation when HERMES denominators vary by variant; it affects exploratory posterior calculations.

**Suggested replacement fragment at line 32:**

```tex
Finite, allele-resolved per-variant case and total sample sizes were retained; rows failing the declared numeric or allele checks were counted and excluded before harmonization.
```

**Suggested addition at line 51:**

```tex
The median-$N$ reduction was a pre-primary implementation rule required by the scalar-$N$ interface; it does not make the contributing per-variant denominators constant, and resulting posteriors are interpreted only as exploratory workflow outputs.
```

The final Results should report the number of invalid/excluded GWAS rows and the range of retained per-variant $N$ by outcome, without treating these as biological findings.

### P1-7 — The current verified ESM staging directory does not contain the final formal dispatcher/lock evidence

**Affected package:** `paper/submission/interdisciplinary_sciences_20260901_external_validation_stage_verified/ESM_1_Code_and_Reproducibility_Package/external_validation/code_and_protocol`.

**Missing or stale items at this checkpoint:**

- `dispatch_parallel_gates.py`;
- `dispatch_ready_models.py`;
- `freeze_model_implementation.py`;
- `model_implementation.lock.json` and `.freeze.json`;
- `test_runner_minimal.R` and `test_coloc_api_ceiling.R`;
- `DOCUMENTATION_CORRECTIONS.tsv`;
- `integrity_freeze_deviation_20260902.tsv`;
- final BLAS serial/parallel comparison and formal dispatch audit files;
- a recorded 20-test output/environment manifest.

**Why this matters:** manuscript lines 49--51 and 95 make reproducibility claims that the currently staged anonymous ESM cannot yet audit. The stage also contains cache/compiled artifacts (`__pycache__`, `.pytest_cache`) that should not substitute for source and test logs.

**Required resolution:** rebuild the anonymous ESM from the locked final source tree after the formal passes but before result-text generation. Include hashes for every deployed script and lock, the two dynamic R test outputs, the external 20-test environment/output, and the amendment/correction ledgers. Exclude caches, local absolute paths and identifying metadata.

## Additional consistency notes

1. The external formal model uses $p_{12}=5\times10^{-6}$, whereas the older synthetic stress experiment at manuscript line 89 uses $p_{12}=10^{-5}$, $L=3$ and 200 iterations. This is not a code contradiction because they are distinct experiments, but one sentence should explicitly say that the locked external settings do not retroactively replace the synthetic-stress settings.
2. The integrity manifests for regions/runtime indices were created after early preparation diagnostics but bind unchanged pre-existing bytes; `integrity_freeze_deviation_20260902.tsv` documents this. Cite that deviation ledger in the Online Resource so “prospectively frozen” is understood as a design/result-naive claim, not a claim that every supplementary manifest existed before any engineering smoke.
3. `DOCUMENTATION_CORRECTIONS.tsv` records a result-discovery glob correction after formal raw preparation but before formal gate/model results. The correction should remain in the anonymous audit trail because it could otherwise look like post-result outcome recovery.
4. The manuscript names pyliftover and UCSC chain files but does not yet cite a stable software/source record for them. Add access URL/date/checksum metadata in the Online Resource and a software citation or repository reference where appropriate.

## Pre-result release gate

The result placeholder should remain until all of the following non-result conditions are satisfied:

1. Complete-case sensitivity is operationally defined or correctly relabelled as an artifact.
2. SuSiE convergence/credible-set defaults are explicitly locked or fully reported.
3. `model_eligible` semantic collision is removed or explicitly documented.
4. Twenty Python tests, two dynamic R tests and the gate integration smoke are reported as separate evidence classes with environments and hashes.
5. Two formal single-thread gate audits and their rerun comparison are copied to durable storage and validated.
6. The final anonymous ESM contains the deployed dispatchers, implementation locks, correction/deviation ledgers, R tests and audit outputs.
7. Only after steps 1--6 should the locked result summary be used to replace `RESULTS_PENDING_LOCKED_RUN`.
