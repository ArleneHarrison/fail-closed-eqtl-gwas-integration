# LD-GuardNet architecture: presentation update v4

The diagram is a deterministic schematic of the existing executable network,
not a new model and not generated performance evidence. The builder imports
the actual SignedLDGuardNet class and verifies both concatenated block input
dimensions, hidden/output dimensions, the classifier input, all 12,772
parameters, and both signed normalization operators.

Figure panels:

- a: Two standardized regional z-score vectors produce five features per
  variant. LD off-diagonal signs define two nonnegative weight matrices;
  self-loops and symmetric degree normalization create the two operators.
- b: Two blocks concatenate self, positive and negative messages and apply
  Linear / LayerNorm / GELU / Dropout. Hidden width is 48.
- c: Node mean and maximum are concatenated to 96 features; a 48-unit head
  produces four logits, followed by softmax. Three independently trained
  probability vectors are averaged. Anomaly score is 1 minus mean valid score.

Four labels are simulated training states, not authenticated diagnoses of a
real-data error. The same architecture is instantiated for each of the three
seed offsets. Dropout is training-only. No attention, learned graph edges,
residual addition, automatic repair or model-eligibility override is depicted
because the source implementation does not contain these operations.

The diagram does not imply that gate dispatch and graph scoring are already
one end-to-end automatic decision pipeline. Structural validation and
probabilistic content warning remain separate evidence layers.

Exports: 170 × 176 mm, vector PDF with embedded fonts, editable-text SVG, and
600 dpi PNG. The manuscript embeds the PDF. No model weights or numerical
results were changed.

Rebuild: python code/render_guardnet_architecture.py

New figure numbering: new architecture is Figure 9; the old v3 Figures 9 and
10 become Figures 10 and 11. Historical filenames elsewhere in Additional
file 3 retain their original numbering and must not be used to infer the
current manuscript sequence.

Official figure guidance checked on 2026-09-20:
https://link.springer.com/journal/12864/submission-guidelines
