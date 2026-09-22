# Review-triggered revision analysis plan

Created on 2026-09-20, after inspection of all v2 results. This is a retrospective
amendment, not a preregistration or a reconstruction of an earlier freeze.

1. Retain original v2 inputs, predictions and 1,800 condition outcomes. Correct
   denominators using an explicit attempted/failed/no-posterior/finite-posterior
   ledger. Run ABF independently of SuSiE on the same seeded data.
2. Recalculate uncertainty from existing predictions with 2,000 bootstrap samples
   of outer replicate IDs (100 clusters, each containing all regimes and states).
   Report the paired AUROC difference against logistic regression. All intervals
   condition on the one training locus and fitted models.
3. Compare the frozen original three-seed network with two information-controlled
   ablations: identical parameter count and training settings, but identity-only
   propagation (non-graph pooled MLP) or identical absolute-LD channels (unsigned
   graph). Use original training/validation seeds only for fitting/early stopping.
   Do not select or suppress a model based on the new challenge outcome.
4. Freeze a new challenge using every available R01--R12 archived regional LD
   matrix (R06 has no matrix). Confirm genomic intervals do not overlap the
   original approximately 2-Mb chromosome-11 training locus. Use 100 evenly spaced
   stored variants per region, with source and selected-matrix hashes. These are
   new regions to the network, but use the same 1000 Genomes EUR donors; this is
   not independent-cohort validation. All z vectors remain simulated.
5. Per challenge region use 40 seeded base draws, cyclically covering four causal
   architectures and N=2,000/20,000, and eight paired states: valid, random z-value
   permutation, local block permutation, random sign flips at 1/5/10/20 percent,
   and permutation of the supplied LD rows/columns with labels held fixed. The
   last is covariance-content corruption, not an ancestry mismatch experiment.
6. Select each anomaly threshold as the 95th percentile of valid-input validation
   scores; freeze it before scoring the new challenge. Report all test strata,
   valid-input false-positive rate and per-fault sensitivity. This is an empirical
   operating rule, not a statistical false-positive guarantee under shift.
7. Run the actual deterministic gate on challenge inputs, separately including an
   explicit identifier-order mismatch. Do not claim z-only permutations are caught
   by the identifier check. No claim of avoided biological conclusions will be made.
8. Add the installed susieR 0.14.2 kriging_rss diagnostic as a binary anomaly
   comparator, using the maximum absolute standardized conditional residual across
   both traits and variants. Select its threshold using the same validation-valid
   quantile. Archive package version, errors and raw scores; do not equate it with
   DENTIST or a four-class predictor. No new test-set tuning is allowed.
9. Retain unfavourable results. If the new network does not outperform appropriate
   controls, describe it as exploratory and do not make architecture-superiority
   or deployment claims. Historical provenance uncertainty and outstanding data
   terms remain explicitly unresolved unless source evidence closes them.
