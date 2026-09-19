# v0.2.0 — multi-signal and signed-LD graph extension

**Archived release DOI:** [10.5281/zenodo.22848192](https://doi.org/10.5281/zenodo.22848192)  
**All-versions concept DOI:** [10.5281/zenodo.22241026](https://doi.org/10.5281/zenodo.22241026)

Version 0.2.0 adds the computational extensions used in the revised manuscript, ``Fail-closed validation with LD-graph anomaly triage for regional eQTL–GWAS integration: multi-tissue benchmarking in coronary artery disease and heart failure.''

## New in this release

- A 1,800-replicate, four-architecture, three-sample-size and three-LD-policy multi-signal robustness benchmark.
- LD-GuardNet, a 12,772-parameter signed-LD graph neural network for optional anomaly triage.
- A leakage-aware three-seed evaluation with 3,200 out-of-distribution graphs from held-out ancestry panels, sample sizes, causal architectures, and random seeds.
- A logistic-regression comparator, bootstrap intervals, calibration output, condition-level predictions, three compact model checkpoints, and a model card.
- Editable PNG/SVG figures, complete environment records, and a separate SHA-256 manifest for the extension.
- Complete manuscript author order: Zhiyong Tao, Yiwei Zhang, Chulin Zhong, Jin Zhou, Wenkai Gao, Yujiao Ran, and Guohua Fan.

## Locked formal results

- Multi-signal attempts: 1,800; completed: 1,752.
- In a no-shared-causal setting at N=20,000, high `coloc.susie` PP.H4 counts changed from 0/50 with matched EUR LD to 19/39 and 26/50 with substituted EAS and AFR LD.
- LD-GuardNet out-of-distribution macro AUROC: 0.943, versus 0.863 for the aggregate logistic baseline.
- LD-mismatch argmax recall: 0.201; anomaly Brier score: 0.164 versus 0.163 for the logistic baseline.

## Interpretation boundary

LD-GuardNet is an optional warning layer, not a deterministic gate or ancestry authenticator. The extension is synthetic and uses one 100-variant region. It does not establish association, empirical colocalization, causal biology, clinical utility, deployment readiness, or general superiority over existing tools.

## Data and licensing

Raw third-party GWAS, eQTL, reference-genotype data, posterior objects, and controlled-access material are excluded. Software is Apache-2.0; project-authored documentation, figures, tables, and generated audit data are CC BY 4.0. See `LICENSE_SCOPE.md`.
