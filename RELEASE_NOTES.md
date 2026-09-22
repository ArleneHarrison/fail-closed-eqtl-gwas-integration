# v0.3.1 — corrected denominators and cross-region validation

**Archived release DOI:** [10.5281/zenodo.22893851](https://doi.org/10.5281/zenodo.22893851)  
**All-versions concept DOI:** [10.5281/zenodo.22241026](https://doi.org/10.5281/zenodo.22241026)

Version 0.3.1 accompanies the manuscript, “An executable fail-closed integrity contract for regional eQTL–GWAS integration: LD-graph triage and cardiovascular workflow evaluation.” It corrects one stale archive-status sentence retained in v0.3.0; software, data, model checkpoints, figures, and numerical results are unchanged.

## New or corrected in this release

- Corrected the multi-signal outcome accounting to 600 base datasets evaluated under three paired LD conditions: 1,800 condition evaluations, 48 nonconvergences, 66 completed evaluations without finite comparable PP.H4, and 1,686 finite posteriors.
- Recomputed ABF results independently for all 600 base datasets and retained a condition-level attempt ledger.
- Recalculated uncertainty with outer-replicate clustering and paired resampling.
- Added an eleven-region synthetic challenge comprising 3,520 content inputs plus 440 explicit identifier-order faults.
- Added unsigned-graph, pooled-MLP, logistic-regression, and SuSiE-RSS conditional-residual comparators, together with graph ablations and validation-selected thresholds.
- Added numerical-interface diagnostics separating projected float64 simulation LD from float32 casting failures.
- Added the code-verified LD-GuardNet architecture schematic in PNG, SVG, and vector PDF formats.
- Added Liqun Zhang to the current manuscript authorship metadata and synchronized the eight-author order.
- Replaced broad reliability and external-validation wording with the narrower integrity-contract, synthetic-warning, and internally frozen workflow-evaluation terminology used by the manuscript.

## Current quantitative boundaries

- LD-GuardNet binary anomaly AUROC in the eleven-region challenge: 0.916.
- Valid-input false alerts at the validation-selected threshold: 57/440.
- Sensitivity for 1% sign errors: 155/440.
- Cardiovascular workflow attempts: 192; 15 with comparable credible-set pairs, 145 completed without comparable pairs, and 32 nonconvergences.

## Interpretation boundary

LD-GuardNet is a synthetic warning prototype, not a deterministic eligibility gate, ancestry authenticator, repair method, or deployment-validated detector. The cardiovascular analysis demonstrates implementation consistency across selected resources; it is not an independent replication or a biological/clinical study. No new disease association, causal mechanism, therapeutic claim, clinical utility, or general superiority over existing tools is asserted.

## Data and licensing

Raw third-party GWAS, eQTL, reference-genotype data, posterior objects, and controlled-access material are excluded. Software is Apache-2.0; project-authored documentation, figures, tables, and generated audit data are CC BY 4.0. See `LICENSE_SCOPE.md`.
