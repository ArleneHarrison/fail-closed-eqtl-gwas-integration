# ESM_1: Code and reproducibility package

This anonymous peer-review package accompanies *An executable fail-closed contract for reliable regional eQTL–GWAS integration*. It contains source code, automated tests, generated fixtures, benchmark outputs, tables, figures, environment records, and checksums. It does not redistribute raw third-party eQTL, GWAS, or reference-genotype files.

## Quick verification

1. Create a Python 3.12 virtual environment.
2. Install `code/requirements_submission.txt`.
3. Change directory to `code`.
4. Run `python -m pytest tests -q -p no:cacheprovider`.
5. Expected result: `93 passed`.

The submission-packaging rerun used Python 3.12.10, NumPy 2.5.2, pandas 2.3.3, SciPy 1.18.0, and pytest 9.1.1. The downstream consumer experiment additionally used the versions recorded in `generated_outputs/downstream_stress/R_sessionInfo.txt`.

## Scope and limitations

The generated outputs support software-integrity, numerical, LD-portability, and controlled downstream-sensitivity claims only. They do not establish biological association, causal inference, clinical utility, or superiority over other tools. The historical fixed-case LD matrix discussed in the manuscript is not available; that specific historical case cannot be numerically reproduced from this package.

## Current licensing (supersedes the historical review-only notice)

The former confidential-review-only restriction in this README is superseded
for project-authored material by the public release's Apache-2.0 software licence
and CC BY 4.0 documentation/generated-material licence. See the package-root
LICENSE_SCOPE.md and licence texts. This correction does not extend any licence
to third-party data or code. The original v0.2.0 archived bytes are unchanged;
this is a corrected local revision package, not a replacement of that DOI record.
