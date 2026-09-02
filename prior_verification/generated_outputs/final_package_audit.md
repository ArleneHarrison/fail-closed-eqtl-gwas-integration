# Final submission-package consistency audit

## Status: **TECHNICAL AUDIT PASSED — AUTHOR COMPLETION REQUIRED**

This is a file-identity and static-placeholder audit of the journal-neutral editable package. It does not select a journal, supply author metadata, or authorize submission.

## Verified identities

- Locked plot-data SHA-256: `3324e9a2a8d3f70011caff09a14e0cec52c74880800be4da9adbde1da8065641`
- Rendered figures: `1`, `2`, `3`, `4`, `5`, `6` (each with PNG, SVG, and immutable manifest)
- Rendered tables: `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8` (each from a locked TSV recorded in the DOCX manifest)
- Figure-manifest SHA-256 values: `{"1": "6caa2a5cdb1050b412545f80330aac11d81ba598f683d4948fbee2a5a23ae9ca", "2": "f55ae84547337619225afa2e05b723eecbfe1602a6d0d1d3ee6ccde5137a047a", "3": "8351584280dfd5f12f3e87799d09cacc02ae6d042deb8ac84e60edb5bf8ebc33", "4": "a8fc3ac5ccf6fe08beaa474c0f787edeb10abfb279f3ab84b3bce221d5aa9dfc", "5": "6ea7df16ea2a7cef68068841585c379862d5085bfeb8109d0c2b71830392b828", "6": "e11b2e1d1b32941b8f52f78f8ff83dc363441450c77359907a8ca977bb254093"}`
- DOCX-manifest SHA-256: `b73dcccf5d5bcfbe35aa91fab345e6bada75b92777d3dbc4d3b8690664394744`
- Editable DOCX SHA-256: `4981a132bc69f8ab70f118a1fbe13e135b3b8f9c4a49329615b0608d336b3ca6`

## Static checks

- PASS — `all_figure_source_hashes_match_their_manifests`
- PASS — `source_table_has_all_thirteen_expected_cases`
- PASS — `contiguous_expected_figure_manifests_present`
- PASS — `all_figure_png_and_svg_hashes_match_manifests`
- PASS — `all_figures_record_restricted_claim_contracts`
- PASS — `docx_embeds_all_manifested_figures`
- PASS — `docx_records_contiguous_verified_tables`
- PASS — `docx_hash_matches_manifest`
- PASS — `docx_names_all_rendered_figures`
- PASS — `docx_preserves_author_placeholders`
- PASS — `docx_labels_author_completion_required`
- PASS — `no_submission_metadata_asserted`

## Visual QA

The editable DOCX was rendered and visually inspected page-by-page. Figures 1-5 document technical contracts, row accounting, numerical diagnostics, policy dependence, and evidence boundaries. Figure 6 is a synthetic downstream-consumer stress test. None is presented as an empirical association, biological result, clinical implication, or comparative-performance claim.

## Remaining hard gates

- Author order, affiliations, corresponding author, and CRediT confirmation.
- Funding, competing interests, ethics/consent wording, and AI-use disclosure confirmation.
- Code/data licence, archive/repository, and derivative-data sharing decisions.
- Target-journal selection and target-specific proof, followed by corresponding-author approval of the complete package.

**Submission instruction:** Complete and verify every author-owned field before upload or submission.
