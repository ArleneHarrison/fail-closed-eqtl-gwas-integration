# External-validation figure fixture

`generate_external_validation_figure_fixture.py` creates a small structural
fixture only. Every table carries `_synthetic_fixture=TRUE`. The production
figure program refuses those inputs unless all of the following are true:

1. `--allow-synthetic-test-output` is supplied;
2. the output directory basename contains `synthetic_test`;
3. the output is visibly watermarked `SYNTHETIC TEST ONLY — NOT FOR SUBMISSION`.

Example smoke test from the project root:

```bash
python scripts/external_validation_20260901/tests/fixtures/generate_external_validation_figure_fixture.py \
  --out-dir tmp/external_validation_figure_fixture

python scripts/external_validation_20260901/make_external_validation_figure.py \
  --summary tmp/external_validation_figure_fixture/summary.tsv \
  --rerun tmp/external_validation_figure_fixture/rerun_agreement.tsv \
  --regions tmp/external_validation_figure_fixture/regions.tsv \
  --genes tmp/external_validation_figure_fixture/genes.tsv \
  --out-dir tmp/external_validation_figure_synthetic_test \
  --expected-regions 4 --expected-attempts 32 \
  --allow-synthetic-test-output
```

These outputs are not manuscript evidence and must never be copied into the
submission package.
