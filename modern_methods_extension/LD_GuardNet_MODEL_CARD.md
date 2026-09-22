# LD-GuardNet model card — review revision v3

## Intended use and status

Exploratory ranking of selected synthetic content anomalies after structural
validation. Not an eligibility gate, ancestry authenticator, biological model
or clinically deployable predictor. No guarantee of input correctness.

## Architecture and original experiment

Three independently seeded 12,772-parameter networks; two positive/negative LD
message blocks, 48 hidden units, mean/max pooling, four-state softmax. Training:
5280 graphs; validation: 1440; original shifted-regime test: 3200 graphs from 800
base data pairs at one locus. Original faults use a fixed value permutation and
ten regularly spaced sign masks. Population panels, sample sizes and causal
architectures shift at test time; genomic region does not.

Original macro AUROC 0.943 (corrected outer-replicate bootstrap 95% CI
0.938–0.948), F1 0.741 (0.725–0.756), accuracy 0.774 (0.763–0.786).
Bootstrap: 2000 resamples of 100 replicate clusters, each with 32 regime/state
graphs. Paired AUROC gain over logistic: 0.080 (0.071–0.089).
LD-mismatch argmax recall 0.201; anomaly Brier 0.164 vs logistic 0.163.
Brier is an overall probability loss, not a pure calibration statistic.

## Review-triggered challenge

All 11 available non-overlapping chromosome-11 reference regions; 100 variants
per region; 40 base draws and 8 paired states per draw = 3520 graphs.
Same EUR reference donors as development; synthetic summary statistics only.
New random/block value permutations, random 1/5/10/20-percent sign errors and
label-preserving LD covariance scrambling. The latter is not ancestry mismatch.
Thresholds use the original validation-valid 95th percentile, not test labels.

| Model | Binary anomaly AUROC | Valid-input false alerts / 440 |
|---|---:|---:|
| Signed graph | 0.916 | 57 |
| Unsigned graph | 0.836 | 110 |
| Pooled non-graph MLP | 0.683 | 245 |
| Logistic | 0.802 | 165 |
| SuSiE kriging regional residual score | 0.761 | 114 |

Signed-graph sign-error sensitivities for 1/5/10/20 percent are
155/440, 343/440, 414/440 and 437/440. The kriging diagnostic flags 311/440
one-percent errors, at a higher false-alert rate. Signed propagation is not
uniformly superior: original four-class macro AUROC was 0.945 for unsigned
propagation versus 0.943 for signed propagation.

## Numerical-interface warning

Structural input validation uses source float64 LD before tensor conversion.
Float32 casting caused numerical gate failure in 9/11 regions. Such failures
are not evidence of detecting false statistics. The precision amendment and
initial metadata-harness errors are retained as engineering evidence, not
primary accuracy outcomes. With declared source64 inputs the gate rejected
440 explicit ID-order faults but passed all 3520 content/value inputs.

## Reproducibility and reporting boundaries

Bundled raw simulated z inputs, selected LD, original checkpoints, ablation
checkpoints, thresholds and replay script reproduce the three network ensembles
to absolute score difference below 2e-6. Historical v2 reports remain explicitly
labelled and must not replace the corrected tables. The revision is an
exploratory, review-triggered amendment, not retrospective preregistration.
External observed-error datasets, independent reference cohorts, calibration at
real anomaly prevalence and downstream utility evaluation remain untested.
