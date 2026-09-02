# Final strict peer-review self-check

## Review scope and limits

This audit follows the supplied AI self-check instruction images and evaluates the final anonymous manuscript, title page, cover letter, figures, tables, lightweight result registries, code, tests, and package structure. It is a computational method-development and software-verification study, not a clinical trial, observational patient-level study, diagnostic-accuracy study, prediction model, meta-analysis, or animal experiment. Consequently, patient flow, exposure allocation, clinical follow-up, sample-size power, CONSORT/STROBE/STARD/TRIPOD/PRISMA/ARRIVE items, clinical effect estimates, and clinical confidence intervals are not applicable.

I could not independently verify original patient-level records, source-study ethics approvals or consents, undistributed third-party raw archives, author identities/contributions, funding, competing interests, submission authorization, HERMES redistribution terms, or a repository DOI/license. These are explicitly retained as author-owned confirmation items; no inference of misconduct is made.

## Study identification

- Type: computational method-development, software verification, and technical transportability evaluation.
- Design: deterministic contract tests, seeded fault injection and simulations, stored-artifact re-audit, LD-reference transfer, retrospective status audit, and an internally frozen 24-region × 4-tissue × 2-outcome factorial technical evaluation.
- Objects: aggregate eQTL/GWAS summary statistics and 1000 Genomes reference genotypes; no newly recruited participants.
- Primary technical endpoints: gate and terminal-code distributions, ordered overlap, LD dimension/rank/condition diagnostics, exact rerun agreement, and version-pinned model terminal states.
- Core conclusion: the implementation makes tested regional-integration rules deterministic and auditable while retaining explicit failure boundaries; it does not establish cardiovascular association, mechanism, causality, tissue specificity, clinical utility, or superiority.

## Critical issues (P0)

No unresolved P0 issue remains in the computational manuscript/package after correction.

Resolved P0 items:

1. Canonical LaTeX entry point, bibliography style, and seven-figure bindings were made unique and self-contained.
2. `regions.tsv` was restored byte-for-byte to the frozen server artifact (SHA-256 `9c6aff019ffcd0d081a78f461760c8d014281435494b4052f2b10e5a7f9d87e9`).
3. Final-mode packaging now requires exactly 192 terminal model rows and validates runner, implementation-lock, protocol, posterior flags, and RDS hashes.
4. ESM anonymization now creates a path-neutral copy and an explicit origin-to-anonymized hash map instead of silently invalidating source hashes.
5. Figure 3's “independent numerical oracle” wording was corrected to “cross-decomposition check”; Figure 7's “prospectively frozen” wording was corrected to “internally frozen”.
6. The packaged prior-verification code/tests were replaced by the complete current audited tree, eliminating five stale tests that depended on absent historical Markdown files.

## Major issues (P1)

The following are not scientific-result defects, but they block actual journal upload until the corresponding author resolves them:

1. Written confirmation from every author of authorship, order, CRediT roles, accountability, final approval, and authorization to submit.
2. Final institutional confirmation of funding and competing-interest statements.
3. Confirmation of HERMES/resource reuse and redistribution terms for the exact intended public release.
4. Public repository URL, immutable archival DOI, and author-approved software license.
5. Confirmation that the manuscript is not under consideration elsewhere and that journal account/metadata are correct.

These items are marked in `AUTHOR_CONFIRMATION_REQUIRED.md`; the release remains `submittable: false` until they are completed.

## Minor issues (P2)

No unresolved formatting or internal-consistency P2 issue was observed in the final PDF. The anonymous manuscript is long (8,109 main-text words; 42 referee-spaced pages including references and display pages), but the length reflects detailed executable-contract methods and audit boundaries. Editorial shortening may be requested, but deleting the denominators, locks, negative results, or limitations would reduce reproducibility.

## Design, novelty, and validity review

- The research question is explicit and matches the methods, results, and conclusions.
- Novelty is framed as an integration-layer executable contract, not as replacement or empirical superiority over GWAS-SSF, MungeSumstats, DENTIST, or SuSiE-RSS.
- The internal freeze is accurately described as occurring before formal primary execution but after disclosed engineering smoke runs; it is not called public preregistration.
- Region and gene selection rules are deterministic and outcome/eQTL eligibility based; posterior results were not used for selection.
- CAD and HF resources, tissue identifiers, coordinate builds, allele orientation, case fractions, median-N reduction, LD construction, gate precedence, model versions, and iteration ceiling are specified.
- The two formal gates are separately executed but use the same preparation implementation; the manuscript does not treat them as independent implementations.
- All 192 real matrices are rank deficient, so the finite `10^12` condition threshold was not exercised; this negative design fact is stated.
- Complete-case LD retained all 503 samples and produced zero differences; the manuscript explicitly rejects a missingness-robustness interpretation.
- Thirty-two HF attempts stopped at the locked 1,000-iteration ceiling. The text conditions this result on susieR 0.14.2, the median-N rule, and the frozen ceiling; it does not generalize to HF biology or all HF analyses.
- No posterior value is interpreted. Completed/no-comparable-set/failure codes are workflow states only.

## Result and numeric consistency review

- Frozen design: 24 regions; 4 tissues; 2 outcomes; 192 attempts.
- Selection: 24 region-level gene selections; 22 unique genes.
- Prepared audit: 389,322 unit-level rows, explicitly denominated as factorially repeated.
- Preparation exclusions: 18 ambiguous reference keys quarantined; none intersected analytic requests; zero discordant analytic duplicates remained.
- Gates: 192/192 READY in each pass; 192/192 exact byte and semantic rerun PASS.
- Overlap: 21–3,807 variants, median 2,245.5.
- LD rank: 19–502; all 192 matrices rank deficient; 188 unique primary and 188 unique complete-case content hashes across 192 bindings.
- Complete case: 503 samples in all attempts; order/rank concordant; all recorded distance/delta metrics exactly zero.
- Formal models: 160 completed with posterior artifacts and RDS hashes; 32 failed without posterior artifacts.
- Terminal codes: 15 `OK`, 145 `OK_NO_COMPARABLE_CREDIBLE_SETS`, 32 `E_SUSIE_NONCONVERGENCE`.
- Outcome split: CAD 96 completed (11/85); HF 64 completed (4/60) plus 32 nonconvergence failures.
- Tissue split: each tissue contributes 48 units and eight HF nonconvergence failures.
- Synthetic consumer experiment: 100 fixed seeded simulations × 5 conditions = 500 condition evaluations; exact-order protection is direct, ancestry protection is conditional on truthful provenance and caller policy, and sign corruption remains unprotected.

The abstract, main text, tables, figure captions, cover letter, terminal registry, and derived count table were checked against these denominators. No conflicting count was retained.

## Reporting, language, and figures

- Title and abstract accurately describe a technical evaluation rather than a cardiovascular discovery study.
- Abbreviations are defined; no P=0.000, “trend” for non-significance, causal language, clinical-benefit language, or unsupported biological claim is used.
- Tables use three-line styling without vertical rules or colored backgrounds.
- Seven main figures have 600-dpi PNG plus editable SVG sources and origin/binding manifests; Figure 7 also has PDF.
- The final anonymous PDF, title page, and cover letter compiled without undefined citations, undefined references, overfull boxes, clipping, overlap, or broken captions.
- All 42 anonymous-manuscript pages, both title-page pages, and the cover-letter page were visually inspected. A remaining Figure 7 label collision and overstatement were corrected and rechecked.

## Reproducibility and release assessment

- Core current source tree: 93/93 tests passed.
- External Python suite: 25/25 tests passed.
- Formal server registry: 192 unique terminal units; 192 model invocations; 160 hashed posterior artifacts; 32 failures with none.
- Raw restricted GWAS/eQTL/genotype files and posterior RDS objects are excluded from ESM; retrieval metadata and hashes are supplied instead.
- ZIP archives contain one canonical `main.tex`, `sn-jnl.cls`, `sn-basic.bst`, bibliography, seven figure bindings, current tests/code, anonymization map, and SHA-256 manifest.

## Reviewer recommendation

The computational/scientific package is suitable for journal submission after the five author-owned P1 items are completed. The principal review risk is novelty/priority rather than an unresolved internal-validity defect: editors may regard the work as a rigorous QC/contract implementation unless they value the locked cross-resource technical evaluation and explicit downstream failure quantification. The manuscript appropriately avoids upgrading technical reproducibility into biological validation.
