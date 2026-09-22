# Supplementary File 3: corrected simulations and review-triggered graph validation

## v4 presentation update

`architecture_v4/` contains the new Figure 9 architecture, source-code/tensor
verification and PDF/SVG/PNG exports. The original v3 Figures 9 and 10 become
Figures 10 and 11 in the current manuscript. Existing historical figure
filenames retain their original numbers. No model or numerical result changed.

## Current versus historical outputs

- `multisignal_coloc/` and `ld_guardnet/` preserve original v2 outputs, including
  the old conditional-ABF missingness and graph-row bootstrap intervals. Their
  historical reports must not be used as the corrected manuscript estimates.
- `multisignal_revised/` is authoritative for the corrected multi-signal result:
  600 base datasets, 1800 paired conditions, 48 SuSiE failures, 66 completed
  outcomes without finite PP.H4, and 1686 finite posteriors. ABF was independently
  recalculated for all 600 base datasets and linked to all conditions.
- `revision_validation/bootstrap_intervals.json` supplies revised outer-replicate
  intervals and the paired difference from logistic regression (2000 resamples).
- `revision_validation/` also contains all 3520 unseen-region challenge outputs,
  input z vectors and LD, six ablation checkpoints, operating thresholds, regional
  scores, kriging comparisons and numerical-interface diagnostics.
- `revision_figures/` contains current figures, editable SVGs, calibration-bin
  tables and provenance. Other embedded figures in historical folders are not
  the current manuscript figures.
- `REVISION_ANALYSIS_PLAN.md` and `NUMERICAL_INTERFACE_AMENDMENT.md` disclose
  review-triggered and engineering amendments. No retrospective preregistration
  is claimed.

## Reproduction without raw third-party data

Use Python 3.12, numpy 2.3.5, pandas 2.3.3, torch 2.11.0 (CPU), scikit-learn 1.9.0
and matplotlib. From the unpacked archive:

```text
python code/replay_revision_predictions.py
python code/summarize_revision_validation.py
```

The replay verifies all three network ensembles against the bundled 3520 scores
(absolute tolerance 2e-6); it neither retrains nor accesses raw GWAS/genotypes.
The full training/challenge builder additionally takes `--workspace` pointing to
the original project tree and uses the exact archived-source paths/hashes in
`challenge_region_manifest.tsv`. Source selection is documented and selected LD
inputs are also bundled. The original five population LD matrices are supplied
in `original_ld_inputs/` for regeneration of the original simulations.

The named diagnostic can be recomputed with R, coloc 5.2.3 and susieR 0.14.2:

```text
Rscript code/revision_kriging_comparator.R revision_validation
Rscript code/revise_multisignal.R code original_ld_inputs multisignal_revised_replay
```

Run into a new output directory to preserve the supplied results. Source scripts
for fresh simulations now compute ABF independently and require explicit
bootstrap groups. The replay input tables carry enough numerical precision to
reproduce the stored network scores without the original regional genotype data.

## Gate and precision results

`source_float64_gate_challenge.tsv` is the corrected contract probe, with separate
source64 and float32 representations. Here source64 is projected float64
simulation LD (eigenvalue floor 1e-8 and diagonal renormalization outside the
gate), not unmodified archived LD. Source64 accepts all 3520 content inputs
and rejects 440 explicit identifier-order faults. Float32 casting introduces PSD
failures in nine regions. This is not evidence of content-error detection.
`actual_gate_challenge.tsv` and `engineering_*_quant_trait_gate.tsv` retain the
initial engineering harness with an unsupported quantitative GWAS metadata label;
do not count those metadata failures as successful anomaly detection. The proper
contract adapter uses declared synthetic case-control metadata and case fraction
0.2; it does not turn simulated Gaussian summary vectors into observed GWAS.

## Scientific boundaries

The challenge spans eleven regions unseen by the graph models, but all share
1000 Genomes EUR reference donors and all summary statistics are simulated.
The original one-region test remains a different experiment. No clinical,
biological association, real-error generalization or deployment claim is made.
The signed graph's new binary AUROC was 0.916, but its valid-input false-alert
rate was 57/440 and 1-percent sign-error sensitivity was 155/440. Unsigned graph
macro AUROC slightly exceeded signed graph macro AUROC in the original task.
All results are retained, without performance-based publication filtering.

## Integrity and licences

`SHA256SUMS.tsv` lists current bytes, excluding itself and transient interpreter
caches. This corrected extension is archived in Zenodo release v0.3.1
(https://doi.org/10.5281/zenodo.22893851).
Project code: Apache-2.0. Project-authored documentation/generated data: CC BY 4.0.
Third-party rights are unchanged. No raw participant genotype or GWAS/eQTL files
are redistributed.
