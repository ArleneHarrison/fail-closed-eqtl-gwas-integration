# Fail-closed regional eQTL–GWAS integration

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22241026.svg)](https://doi.org/10.5281/zenodo.22241026)

Versioned code and reproducibility materials for:

> **An executable fail-closed contract for regional eQTL–GWAS integration: a multi-tissue technical evaluation in coronary artery disease and heart failure**

This repository implements an auditable, fail-closed contract for preparing and evaluating regional eQTL–GWAS analyses. It treats coordinate identity, allele orientation, ordered variant overlap, linkage-disequilibrium provenance, matrix integrity, metadata completeness, and deterministic reruns as explicit preconditions. Inputs that do not satisfy the contract terminate with recorded states instead of being silently repaired or passed downstream.

The cardiovascular evaluation uses coronary artery disease and heart-failure resources to study workflow behaviour across multiple tissues. It is a **technical validation**, not a report of new biological associations, causal effects, drug targets, or clinical utility.

## What is included

- `prior_verification/`: implementation, 93-test verification suite, fixtures, simulations, numerical checks, and generated audit outputs used before the locked cardiovascular evaluation.
- `external_validation/`: prospectively frozen multi-tissue cardiovascular protocol, locked implementation, 25-test validation suite, source/provenance audits, terminal registries, and reporting artifacts.
- `PEER_REVIEW_PACKAGE_SHA256_MANIFEST.tsv`: the immutable 323-record checksum manifest of the peer-review supplementary ZIP from which this public repository was prepared. It verifies that original ZIP after extraction; public-release metadata intentionally changes the root and prior-verification README files.
- `RELEASE_SHA256_MANIFEST.tsv`: repository-wide checksums for the public release, excluding Git metadata and the manifest itself.

Formal model results are represented by terminal registries and artifact hashes. Raw third-party eQTL, GWAS, and reference-genotype files, as well as posterior objects, are not redistributed.

## Quick verification

### Prior verification suite

Use Python 3.12, then:

```bash
cd prior_verification/code
python -m pip install -r requirements_submission.txt
python -m pytest tests -q -p no:cacheprovider
```

Expected result: `93 passed`.

### External validation suite

From `external_validation/code_and_protocol`:

```bash
python -m pytest -q -p no:cacheprovider
Rscript tests/test_coloc_api_ceiling.R
Rscript tests/test_runner_minimal.R run_coloc_susie_rank_diagnostic.R
```

Expected result: `25 passed` and both R checks pass. The locked R checks require the package versions declared by the protocol and report a closed failure when those requirements are not met.

The release was checked with Python 3.12.10. Exact Python and R package versions are recorded in the included environment and session-information files.

## Reproducing the cardiovascular evaluation

The execution contract and commands are documented in [`external_validation/code_and_protocol/README.md`](external_validation/code_and_protocol/README.md). Reproduction with the original resources requires users to obtain the cited third-party datasets under their own access terms and verify each source against the locked provenance metadata. The pipeline does not download or redistribute restricted source data automatically.

The frozen design contains 24 disease regions, four tissues, two outcomes, and 192 attempts. Only attempts satisfying all preconditions can enter the READY state and reach the locked model runner. A stopped attempt cannot contain a posterior artifact.

## Interpretation boundary

The included evidence supports claims about software integrity, deterministic gating, numerical behaviour, LD portability, and controlled downstream sensitivity. It does not establish biological association, causality, clinical validity, clinical utility, or general superiority over existing tools. See the protocol and audit reports for the prespecified failure states and known limitations.

## Data availability

This repository contains project-authored code, synthetic fixtures, configuration, provenance records, aggregate audit outputs, figures, tables, and hashes. It does **not** contain raw or row-level third-party GWAS, eQTL, or genotype-reference data. Dataset identifiers, retrieval metadata, and source-specific access information are retained so eligible users can reacquire the inputs from their original custodians.

## Licensing

- Software source, tests, scripts, and executable configuration are licensed under the [Apache License 2.0](LICENSE).
- Project-authored documentation, figures, tables, and generated audit data are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt).
- Third-party datasets and software remain under their original terms and are not relicensed here.

The exact file-level scope is described in [`LICENSE_SCOPE.md`](LICENSE_SCOPE.md) and `.reuse/dep5`.

## Citation

Please cite the archived release:

> Tao, Z., Zhang, Y., Zhong, C., Zhou, J., & Fan, G. (2026). *An executable fail-closed contract for regional eQTL–GWAS integration: a multi-tissue technical evaluation in coronary artery disease and heart failure* (Version 0.1.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22241026

The concept DOI for all versions is `10.5281/zenodo.22241026`. Version 0.1.2 restores the complete author order: Zhiyong Tao, Yiwei Zhang, Chulin Zhong, Jin Zhou, and Guohua Fan. The version-specific DOI is added here after Zenodo archives the release.

## Integrity

To verify the public release from the repository root in PowerShell:

```powershell
Import-Csv .\RELEASE_SHA256_MANIFEST.tsv -Delimiter "`t" | ForEach-Object {
  if ((Get-FileHash -Algorithm SHA256 -LiteralPath $_.path).Hash.ToLower() -ne $_.sha256) {
    throw "Checksum mismatch: $($_.path)"
  }
}
```

## Contact and contributions

Questions and reproducible bug reports are welcome through GitHub Issues. Please do not upload controlled-access or third-party row-level data to an issue or pull request.
