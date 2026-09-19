# LD-GuardNet model card

## Intended use

Rank selected regional summary-statistic/LD inputs for manual review after deterministic validation. The modeled states are valid input, summary-order mismatch, 10% sign corruption, and LD mismatch.

## Architecture

The 12,772-parameter network applies two signed-LD message-passing blocks, treating positive and negative LD as separate normalized operators. It pools node embeddings by their mean and maximum and returns four softmax scores. Three independently seeded models are averaged for the formal result.

## Formal split

- Training: 5,280 graphs; EUR/EAS/AFR; N=5,000; one- and two-shared-signal architectures.
- Validation: 1,440 new-seed graphs under the training regimes.
- Out-of-distribution test: 3,200 graphs; held-out SAS/AMR panels; N=2,000 or 20,000; shared-plus-specific or distinct-only architectures.

Valid graphs used ancestry-matched LD from every data-generating panel, preventing panel identity alone from defining anomaly status. The out-of-distribution ancestry panels were not used for training or model selection.

## Formal result

- Macro one-vs-rest AUROC: 0.943 (bootstrap 95% CI 0.938-0.948).
- Macro F1: 0.741 (0.730-0.754).
- Accuracy: 0.774 (0.760-0.789).
- Logistic baseline macro AUROC: 0.863.
- Class AUROC: valid 0.902; order mismatch 1.000; sign corruption 0.999; LD mismatch 0.872.
- Argmax recall: valid 0.935; order mismatch 1.000; sign corruption 0.960; LD mismatch 0.201.
- Anomaly Brier score: 0.164 versus 0.163 for the logistic baseline.

## Limitations and prohibited uses

The formal test is synthetic and uses one genomic region. Low thresholded LD-mismatch recall and lack of calibration improvement preclude use as a deterministic gate or ancestry authenticator. Do not use the scores as evidence of association, colocalization, causal biology, clinical utility, or input correctness. External cross-locus and real-error validation is required before operational deployment.

