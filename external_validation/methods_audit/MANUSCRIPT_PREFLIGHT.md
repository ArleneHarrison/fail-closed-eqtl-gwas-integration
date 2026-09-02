# Manuscript preflight: Interdisciplinary Sciences: Computational Life Sciences

## Review scope and current disposition

- Reviewed file: `paper/submission/interdisciplinary_sciences_20260901/latex_source/main_external_validation_draft.tex`
- Snapshot SHA-256: `b5311957e9b3c48b46f740188b4207e5c177526e7c042090ecdf2d182b9f5e70`
- Review mode: strict, read-only peer-review-style audit. The manuscript text and result placeholder were not changed.
- Reliable checks completed: internal consistency of the supplied TeX/BibTeX and companion submission files; clean LaTeX compilation; comparison with the journal's current official submission instructions.
- Not independently verified here: raw-data truth, final locked-run outputs, ethics documents from source cohorts, authorship/contributor agreements, software-license rights, or the validity of every third-party record.

**Disposition:** the scientific framing is substantially safer than the earlier biological-discovery framing, but this is not yet a submission file. One P0 scientific blocker and one P0 package blocker remain. The placeholder must stay until the locked run is complete; no numerical value should be inferred or inserted manually.

## Must fix before submission

### P0-1 — The abstract, Results, Discussion and Conclusions still contain a result placeholder

**Locations:** lines 14, 91--94, 316 and 391.

**Problem:** `RESULTS_PENDING_LOCKED_RUN` and forward-looking editorial text are appropriate safeguards during analysis, but a manuscript containing them is incomplete. The abstract currently describes what the final abstract *will* report rather than reporting a result. The same unresolved state is repeated in Results, Discussion and Conclusions.

**Required action:** keep the existing block unchanged until the frozen input manifest, protocol hash, run manifest, exact-rerun comparison and result-table checks pass. Then update all four locations from the same locked summary table in one operation. Do not report selected loci or posterior examples unless the prespecified reporting rule allows them.

**Exact result-sentence template (tokens are deliberately not values):**

```tex
\textbf{Results:} Of 192 prespecified attempts, [N_READY] were READY and [N_INELIGIBLE] were INELIGIBLE; the leading terminal codes were [CODE: N] and [CODE: N]. Ordered overlap was [PRESPECIFIED SUMMARY], and LD rank was [PRESPECIFIED SUMMARY]. Canonical primary outputs agreed in [N]/[N] exact rerun comparisons. [N_DISPATCHED] READY attempts were dispatched and [N_COMPLETED] completed the exploratory downstream model.
```

Report all denominators, including zero-event categories. If no attempt is READY, state that no downstream posterior was produced. Retain the existing non-biological interpretation in the final conclusion.

### P0-2 — The upload directory still points to the old manuscript

**Locations:** `latex_source/README.txt` lines 1--4; `latex_source/main.tex`; existing `latex_source/main.pdf`; companion title page and cover letter.

**Problem:** the README instructs the editor to compile `main.tex`, but the revised work is in `main_external_validation_draft.tex`. The existing `main.tex` and `main.pdf` are the pre-validation version. Uploading the directory as-is could submit the wrong manuscript even if the draft is scientifically completed.

**Required action:** build a clean final submission directory. In that directory only, rename the locked revised source to `main.tex`, compile its PDF, include both `.bib` files, class/BST files and final figures, and exclude stale `.aux`, `.bbl`, `.blg`, `.fls`, `.fdb_latexmk`, old PDFs and the revision-stage draft. Update the README to name the final entry point. Verify the final PDF and source archive, not the current old `main.pdf`.

### P1-1 — READY is described both as non-model eligibility and as sufficient for automatic model dispatch

**Locations:** lines 38 and 43; also lines 14 and 94.

**Problem:** line 43 correctly says `READY` is policy conformance only and that model eligibility remains unset. Line 38 then says only READY attempts can be dispatched, which can be read as READY being sufficient for SuSiE-RSS/`coloc.susie`. The external model-dispatch checks, package defaults, non-convergence rule and posterior reporting rule are not stated. This is the main methods/claims inconsistency.

**Replace the first two sentences of line 38 with:**

```tex
The frozen factorial comprises 24 regions $\times$ four tissues $\times$ two outcomes, or 192 attempts. Prespecified primary endpoints are gate-state and terminal-code distributions, ordered variant overlap, LD dimension, numerical rank and condition number, and exact canonical-output rerun agreement. READY is a necessary but not sufficient condition for model dispatch. A separate version-pinned dispatch step checks the complete model input and records package version, model settings, convergence status and any model-level stop reason; only attempts passing both steps may reach SuSiE-RSS and \texttt{coloc.susie}. Non-READY attempts and model-dispatch failures produce no posterior.
```

Then add the exact `susieR` and `coloc` versions, settings, priors/defaults, convergence rule and failure codes actually used. If the frozen protocol did not prespecify a choice, label it as a recorded software default or exploratory rule rather than retrospectively calling it prespecified.

### P1-2 — The “prospectively frozen” account omits the v1-to-v2 amendment

**Locations:** lines 31--34 and 316.

**Problem:** line 32 says “Before running the new primary data, we froze...” The auditable record instead contains a v1 protocol and a metadata-resolution v2 amendment frozen before any external-evaluation result was viewed. Omitting that amendment makes the chronology less transparent and risks a reviewer construing metadata resolution as post-result adaptation.

**Replace the first sentence of line 32 with:**

```tex
Before viewing any external-evaluation result, we froze amended protocol v2 (`cv_external_validation_20260901_v2`) at 2026-09-01T17:33:29Z with canonical SHA-256 `ec619b392221037bf63628919b22bba0b45d9f1c0ccd720b9ec52531e0bf0ba6`. Amendment A001 resolved prespecified ancestry, genome-build and case-fraction metadata stops using authoritative source records; it did not alter the factorial design or selection rules, and the amendment ledger records that no external result had been viewed.
```

The protocol filename, both hashes and amendment ledger must be in the anonymous supplementary package.

### P1-3 — “External” can still be mistaken for sample-independent replication

**Locations:** lines 14, 31, 40, 316 and 318.

**Problem:** the limitations are commendably explicit, but “external evaluation” is still the dominant short label. The CAD files overlap through UK Biobank; the four GTEx tissues share a donor framework; per-tissue ancestry is unconfirmed; and the tissue labels appeared in exploratory code. “External” should not carry the headline meaning.

**Required terminology replacement:**

- Section heading at line 31: `\subsection{Prospectively frozen multi-tissue, multi-resource technical evaluation}`
- Abstract at line 14: replace “external technical evaluation” with “prospectively frozen multi-resource technical evaluation”.
- Keep the detailed externality paragraph at line 40 and the limitation at line 318.
- Reserve the phrase “external technical” for the formal protocol label, immediately qualified as resource/workflow-level and not participant-independent.

The current title is acceptable because it already says “multi-tissue technical evaluation” and does not claim independent validation.

### P1-4 — New-resource sample-size and case-fraction handling is under-specified

**Locations:** line 32, with policy context at lines 43--50.

**Problem:** the CAD file has fixed study totals, whereas HERMES2 supplies per-variant counts and the quoted 139,533/1,568,809 values are maximum aggregate counts. The text says per-variant counts are authoritative but does not state which values enter the gate and which single study-level case fraction enters a case-control model. Readers cannot reproduce denominator handling.

**Replace the HERMES2 sentence at line 32 with:**

```tex
The heart-failure resource is HERMES2 European overall heart-failure phenotype 1 on GRCh37/b37. Its maximum aggregate counts are 139,533 cases and 1,568,809 controls ($N_{\max}=1{,}708{,}342$). Deposited per-variant case and total sample sizes were retained and audited row by row; where a downstream case-control model required one case fraction, the pre-run metadata lock supplied the study-level value 0.0816774392949421 derived from the aggregate maxima rather than recomputing a varying row-level fraction.
```

Make the analogous CAD sentence state that 34,541/261,984 applies to the exact UK Biobank-only file. Do not call either denominator an independent cohort sample size.

### P1-5 — External-data preprocessing is not reproducible enough from the manuscript

**Locations:** lines 34--38.

**Missing items:**

1. The exact hg38-to-hg19 chain filename, retrieval source, checksum and lift-over software/version.
2. The eQTL Catalogue `an` handling: effective row-level $n=\mathrm{an}/2$ only after a positive-even audit.
3. Separate construction of eQTL MAF, GWAS MAF and LD-reference MAF; no single MAF should silently stand in for all three.
4. Duplicate-key detection and precedence before any join.
5. The exact ordered-overlap fraction denominator.
6. The LD extraction/correlation software and version, biallelic/variance filters, dosage orientation, and how the 2% missingness ceiling was evaluated.
7. The exact condition-number threshold ($10^{12}$ in protocol v2).
8. The role of complete-case-sample LD: sensitivity diagnostic versus model-dispatch stop, and the metric used to compare it with primary LD.
9. The terminal-code precedence when more than one defect occurs.

**Required action:** add a compact “External preprocessing and dispatch” paragraph covering these points and cite the complete protocol/manifest as an Online Resource. Values must be read from the locked configuration and environment record, not reconstructed from memory.

### P1-6 — The region registry contains overlapping outcome-anchored windows, but independence is not addressed

**Locations:** lines 36, 38 and 40; future Results line 94.

**Problem:** CAD- and HF-selected windows can overlap and are intentionally retained as distinct outcome anchors. Likewise, tissue attempts share GTEx donors and outcome summary data. Treating 192 attempts as independent observations would be invalid.

**Append to line 36:**

```tex
Windows were selected separately within outcomes; overlapping CAD- and HF-anchored windows were retained as distinct prespecified anchors and were not treated as independent regions.
```

In Results, report exact descriptive counts and distributions. Do not attach binomial confidence intervals, independence-based $P$ values, or tissue-ranking claims to the 192 correlated attempts.

### P1-7 — The old “external biological validation not performed” row needs a sharper boundary, not deletion

**Locations:** Table 2 line 173 and Fig. 5 caption line 180; related text at lines 29, 81 and 334.

**Problem:** “External biological validation — Not performed” remains scientifically true, but beside a new multi-resource technical evaluation it reads as a contradiction. The distinction should be explicit at the point of use.

**Replace Table 2 line 173 with:**

```tex
Biological or sample-independent cardiovascular validation & Not performed & No biological/replication endpoint & The multi-resource CAD/HF run evaluates workflow transportability only \\
```

**Replace the last clause of the Fig. 5 caption with:**

```tex
...; the prospectively frozen CAD/HF analysis is a workflow-level technical evaluation, and no external biological or sample-independent validation was performed.
```

After the run, add the technical evaluation as its own row in the evidence hierarchy at Table 4 (lines 338--363), with observed counts populated only from the locked summary and “replicated association, tissue specificity and sample independence” in the “Not supported” column.

### P1-8 — The abstract contains undefined model/status terms and meta-drafting language

**Location:** line 14.

**Problem:** the journal requires a 150--250-word self-contained abstract without undefined abbreviations. The present draft is approximately 188 words, but `READY`, SuSiE-RSS and `coloc.susie` are not defined, and “The final abstract will report...” is not manuscript prose.

**Required action after results lock:** keep the abstract at 150--250 words and either define or remove these names. Preferred concise wording: “Only inputs passing the technical gate were eligible for downstream fine-mapping and colocalization.” Keep 4--6 keywords; the current six satisfy the numerical requirement.

### P1-9 — Core test counts are already stale once the new external pipeline is included

**Locations:** line 84; Table 2 line 163; line 144; conclusion line 391; companion cover letter line 22.

**Problem:** “full 93-test suite” refers to the pre-existing core package. The external pipeline has a separate test suite and must not be folded into 93 without rerunning and recounting. Environment versions also need to cover the server-side pipeline and R packages, not only the local Python core.

**Replacement pattern:**

```tex
The core implementation suite contained 93 tests and passed 93/93 in the recorded environment. The external-evaluation pipeline contained [N_EXTERNAL_TESTS] tests and passed [N_EXTERNAL_PASSED]/[N_EXTERNAL_TESTS] in [ENVIRONMENT].
```

Update every repeated count from one generated verification summary. Report 500 input-condition replicates as “500 SuSiE-RSS fits and 500 coloc.abf evaluations,” not the ambiguous “500 consumer runs.”

### P1-10 — Resource/software citations and provenance are incomplete for the added analysis

**Locations:** lines 32--38, 84 and bibliography line 406; `references_external_validation.bib`.

**Problem:** the papers are cited, but the exact HERMES2 file/landing record, eQTL access/licence pages, chain file, 1000 Genomes VCF release, and versioned software are not all citable from the manuscript. The CVD Knowledge Portal paper is not a substitute for the exact HERMES2 data-file record. HERMES banner authorship/reuse instructions also require author confirmation before submission.

**Required action:**

- Cite the exact HERMES2 dataset landing page or stable repository record in addition to the Nature Genetics paper.
- Cite/accession the four QTD files and exact CAD file; give retrieval dates and SHA-256 hashes in Online Resource metadata.
- Add version/citation records for lift-over tooling, genotype/LD tooling, `susieR`, `coloc`, and the external Python/R runtime actually used.
- Verify and obey HERMES data-banner acknowledgement/authorship instructions; this cannot be inferred from the manuscript.
- In the compiled bibliography, confirm DOI entries render as full `https://doi.org/...` links as requested by the journal.

### P1-11 — Gate table is presented as comprehensive but omits new terminal states

**Locations:** Table 1 lines 101--129; future external Results line 94.

**Problem:** the new workflow can terminate before or at ordered overlap (at least `E_NO_ELIGIBLE_GENE` and `E_NO_ORDERED_OVERLAP`), but Table 1 omits these states. A reader cannot map all external counts to the displayed semantics.

**Required action:** either change the caption to “Selected gate decisions...” and direct readers to the full codebook, or add every code that can appear in the 192-attempt registry, including selection, lift-over/allele, zero-overlap, LD and model-dispatch codes. Keep selection-terminal codes distinct from gate-terminal codes.

### P1-12 — The evidence hierarchy and narrative will be incomplete after adding the new results

**Locations:** Table 4 lines 338--363; Discussion lines 316--336; Conclusions line 391.

**Problem:** Table 4 currently stops at the old generated/stored evidence. When real technical-evaluation results are inserted at line 94, the evidence table must add the same component and the Discussion must interpret both successes and failures without selectively emphasizing READY cases.

**Required action:** add one row for the frozen multi-resource technical evaluation and synchronize its observed counts with the result registry. In Discussion, explain the dominant terminal states, whether failures clustered by resource/tissue/region, and what each pattern says about file/workflow compatibility. Explicitly state that tissue differences may arise from sample size, donor correlation, ancestry mismatch or variant availability.

### P1-13 — The submission companion files are materially stale

**Locations:** `02_Title_Page.tex` lines 12, 22--24 and 53--60; `03_Cover_Letter.tex` lines 18--24; `_ESM_1_stage/README.md` line 3.

**Problem:** these files retain the old title, word counts, display counts, resource list and study summary. The title-page availability statement names only QTD000216 and GCST90132314, not the four QTD datasets, exact CAD file, HERMES2 file, amended protocol, or external outputs. The cover letter also claims only the pre-existing evidence.

**Required action:** update these files only after the manuscript/result registry is locked. Use generated word/display counts. Preserve author identifiers and declarations only on the separate title page. Keep the anonymous manuscript, figures and review-stage supplementary archive free of author names, affiliations, e-mail addresses and identifying PDF/ZIP metadata.

## Strongly recommended revisions

### P1-14 — Remove one premature table citation and reduce repeated boundary language

**Locations:** line 21; repeated caveats at lines 25, 29, 40, 81, 94, 173, 180, 216, 228, 253, 280, 310, 316, 318, 326--336, 363, 388 and 391.

**Problem:** line 21 cites the final scope table before Tables 1--4 are introduced, contrary to the journal request for consecutive table citation. The same “not biological/causal/clinical/superiority” boundary is then repeated in nearly every subsection, figure/table note and conclusion. The caution is appropriate but the repetition makes the paper longer and can obscure the actual contribution.

**Replace the last sentence of line 21 with:**

```tex
The relevant methods operate at adjacent workflow layers; the descriptive scope comparison is presented after the empirical evidence and is not a head-to-head performance comparison.
```

Retain the full boundary statement in four places: abstract, externality Methods paragraph, evidence-hierarchy table/note, and final Discussion/Conclusion. Shorten the remaining notes to the single limitation specific to that figure or table.

### P1-15 — Avoid an untested uniqueness claim in the Background

**Location:** line 21.

**Problem:** “The present contribution is the integration ... into one ... contract” can be read as claiming uniqueness, although no systematic review or head-to-head implementation benchmark was done.

**Replace with:**

```tex
Our implementation places provenance, deterministic row selection, exact LD identity and order, trait metadata, and caller-owned numerical policy within one machine-readable pre-inference contract.
```

### P1-16 — Mathematical notation should be typeset unambiguously

**Locations:** lines 53--55.

**Problem:** `p`, `lambda`, “max absolute lambda” and “min absolute lambda” are written as prose, which is less precise than the rest of the numerical methods.

**Suggested replacement core:**

```tex
For dimension $p$ and eigenvalues $\lambda_1,\ldots,\lambda_p$, the tolerance is $p\,\epsilon_{\mathrm{mach}}\max_i|\lambda_i|$. Numerical rank is $\sum_i I(|\lambda_i|>\mathrm{tol})$; for a numerically full-rank symmetric matrix, $\kappa_2=\max_i|\lambda_i|/\min_i|\lambda_i|$.
```

### P1-17 — Figure delivery should use the editable/vector assets already available

**Locations:** figure includes at lines 133, 139, 179, 221, 227 and 285; final new external-validation figure.

**Problem:** the TeX directory currently includes only PNG versions. The journal prefers vector files for line art and requires 1200 dpi for bitmap line art or 600 dpi for combination art, RGB for color, legible final-size lettering and accessibility beyond color alone.

**Required action:** include EPS/PDF vector versions for diagrams and plots where possible, or document that a raster is combination art at 600 dpi. Check embedded fonts, final 174-mm width, contrast and pattern/shape redundancy. Add descriptive captions/alt text for the new figure.

## Items that currently pass

- The current title is specific, cautious and cardiovascularly relevant; it does not claim biological discovery or participant-independent replication.
- The current abstract length is within the journal's 150--250-word range, and six keywords satisfy the 4--6 requirement. It still needs the P0/P1 content fixes above.
- The anonymous manuscript uses the Springer Nature LaTeX class, referee spacing, line numbers and numbered-reference style, and it compiled without unresolved citations/references or overfull boxes in this preflight.
- The separate title page architecture is appropriate for the journal's double-blind procedure.
- The manuscript repeatedly avoids association, causal, therapeutic and clinical claims, and it explicitly discloses donor/resource dependence and ancestry uncertainty.
- The AI-use statement at lines 86--87 follows the journal instruction to document substantive LLM use in Methods and preserves human accountability.
- Deterministic fixture counts are correctly described as exact implementation checks rather than sensitivity/specificity estimates; no confidence interval is implied for a non-random fixture set.

## Final release gate

Do not create the submission ZIP until all of the following are true:

1. Frozen v2 manifest and input hashes verified.
2. All 192 attempts accounted for exactly once, with zero-overlap and no-gene states retained.
3. Exact rerun comparison passed or every discrepancy explained.
4. Abstract, Results, Discussion, Conclusions, evidence hierarchy and external figure/table generated from the same locked summary.
5. Core and external test suites rerun with separate counts and complete Python/R/package versions.
6. HERMES acknowledgement/authorship/reuse instructions confirmed by the authors.
7. Clean final `main.tex` and PDF built from the revised source; no old `main.tex`/PDF or placeholder remains.
8. Title page, cover letter, Online Resources and manifests updated to the same title, counts and resource identifiers.
9. Anonymous PDF/source/figures/ESM checked for author names, e-mail addresses, filesystem paths and identifying metadata.
10. Final PDF visually inspected page by page and the source archive test-compiled in an isolated directory.

## Journal-format basis

The journal's official submission page (checked 2026-09-01) specifies double-blind review, a separate title page, editable source files, a 150--250-word abstract, 4--6 keywords, numbered citations, Springer Nature LaTeX for mathematically oriented manuscripts, consecutive table/figure citation, and explicit supplementary-material captions. Source: <https://link.springer.com/journal/12539/submission-guidelines>.
