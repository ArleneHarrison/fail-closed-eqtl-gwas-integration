# Additional file 3: multi-signal and LD-graph extension

This archive contains the fully generated inputs, outputs, source code, model checkpoints, figures, and environment records for two synthetic extensions reported in the manuscript.

## Contents

- `code/run_multisignal_coloc_robustness.R`: runs the four-architecture, three-sample-size, three-LD-policy benchmark.
- `code/render_multisignal_robustness_figure.R`: regenerates the multi-signal figure from the summary table.
- `code/run_ld_guardnet_benchmark.py`: generates the leakage-aware train/validation/out-of-distribution splits, trains LD-GuardNet with three seeds, fits the logistic baseline, evaluates both models, and renders the graph-model figure.
- `multisignal_coloc/`: 1,800 attempted condition rows, aggregate results, session information, manifest, and PNG/SVG figure.
- `ld_guardnet/`: formal benchmark report, 3,200 out-of-distribution predictions, seed-level metrics, three compact model checkpoints, training history, and PNG/SVG figure.
- `SHA256SUMS.tsv`: byte-level integrity manifest for every archived file except the manifest itself.

## Scope boundary

All examples are synthetic and use the same 100-variant genomic region. The experiments quantify software-consumer sensitivity and anomaly-triage behavior. They do not provide evidence of cardiovascular association, empirical colocalization, biological mechanism, causal effect, clinical validity, or deployment readiness.

LD-GuardNet is an optional probabilistic warning layer. It cannot override a deterministic contract failure, authenticate ancestry, or certify inputs that it classifies as valid. The formal out-of-distribution run met the prospectively locked inclusion rule, but default-threshold LD-mismatch recall remained low and calibration did not improve over the logistic baseline.

## Reproduction

Run the scripts from the project root and pass the directory containing `common_ld_EUR.tsv`, `common_ld_EAS.tsv`, `common_ld_AFR.tsv`, `common_ld_SAS.tsv`, and `common_ld_AMR.tsv`. The scripts accept only explicit input and output paths and record package versions, seeds, settings, and interpretation boundaries in their generated manifests.

