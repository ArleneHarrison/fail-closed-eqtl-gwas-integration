# Prior verification package

This directory contains the code and reproducibility materials used to verify the fail-closed regional eQTL–GWAS integration contract before the locked cardiovascular technical evaluation. It includes source code, automated tests, synthetic fixtures, benchmark outputs, tables, figures, environment records, and checksums. It does not redistribute raw third-party eQTL, GWAS, or reference-genotype files.

## Quick verification

1. Create a Python 3.12 virtual environment.
2. Install `code/requirements_submission.txt`.
3. Change directory to `code`.
4. Run `python -m pytest tests -q -p no:cacheprovider`.
5. Expected result: `93 passed`.

The submission-packaging rerun used Python 3.12.10, NumPy 2.5.2, pandas 2.3.3, SciPy 1.18.0, and pytest 9.1.1. The downstream consumer experiment additionally used the versions recorded in `generated_outputs/downstream_stress/R_sessionInfo.txt`.

## Scope and limitations

The generated outputs support software-integrity, numerical, LD-portability, and controlled downstream-sensitivity claims only. They do not establish biological association, causal inference, clinical utility, or superiority over other tools. The historical fixed-case LD matrix discussed in the manuscript is not available; that specific historical case cannot be numerically reproduced from this package.

## Licensing

Project-authored software in this directory is licensed under Apache-2.0. Project-authored documentation, figures, tables, and generated audit data are licensed under CC BY 4.0. Third-party materials remain subject to their original terms. See the repository-root `LICENSE_SCOPE.md` for details.
