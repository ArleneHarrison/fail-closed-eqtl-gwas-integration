args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: run_coloc_susie_rank_diagnostic.R <input_dir> <case_fraction> <out_dir>")
}

input_dir <- args[[1]]
case_fraction <- as.numeric(args[[2]])
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

EXPECTED_COLOC_VERSION <- "5.2.3"
EXPECTED_SUSIER_VERSION <- "0.14.2"
IMPLEMENTATION_SEED <- 20260901L
SUSIE_L <- 10L
SUSIE_MAX_ITER <- 1000L
SUSIE_REPEAT_UNTIL_CONVERGENCE <- FALSE
SUSIE_ESTIMATE_RESIDUAL_VARIANCE <- FALSE
SUSIE_TOL <- 0.001
SUSIE_COVERAGE <- 0.95
SUSIE_MIN_ABS_CORR <- 0.5
SUSIE_SCALED_PRIOR_VARIANCE <- 0.2
SUSIE_CHECK_PRIOR <- TRUE
SUSIE_Z_LD_WEIGHT <- 0
COLOC_P1 <- 1e-4
COLOC_P2 <- 1e-4
COLOC_P12 <- 5e-6
COLOC_OVERLAP_MIN <- 0.5
COLOC_TRIM_BY_POSTERIOR <- TRUE
COLOC_BACK_CALCULATE_LBF <- FALSE
set.seed(IMPLEMENTATION_SEED)

actual_coloc_version <- if (requireNamespace("coloc", quietly = TRUE)) {
  as.character(packageVersion("coloc"))
} else {
  "NOT_INSTALLED"
}
actual_susier_version <- if (requireNamespace("susieR", quietly = TRUE)) {
  as.character(packageVersion("susieR"))
} else {
  "NOT_INSTALLED"
}

write.csv(data.frame(
  expected_coloc_version = EXPECTED_COLOC_VERSION,
  actual_coloc_version = actual_coloc_version,
  expected_susieR_version = EXPECTED_SUSIER_VERSION,
  actual_susieR_version = actual_susier_version,
  susie_L = SUSIE_L,
  seed = IMPLEMENTATION_SEED,
  max_iter = SUSIE_MAX_ITER,
  repeat_until_convergence = SUSIE_REPEAT_UNTIL_CONVERGENCE,
  estimate_residual_variance = SUSIE_ESTIMATE_RESIDUAL_VARIANCE,
  tol = SUSIE_TOL,
  coverage = SUSIE_COVERAGE,
  min_abs_corr = SUSIE_MIN_ABS_CORR,
  scaled_prior_variance = SUSIE_SCALED_PRIOR_VARIANCE,
  check_prior = SUSIE_CHECK_PRIOR,
  z_ld_weight = SUSIE_Z_LD_WEIGHT,
  coloc_p1 = COLOC_P1,
  coloc_p2 = COLOC_P2,
  coloc_p12 = COLOC_P12,
  coloc_overlap_min = COLOC_OVERLAP_MIN,
  coloc_trim_by_posterior = COLOC_TRIM_BY_POSTERIOR,
  coloc_back_calculate_lbf = COLOC_BACK_CALCULATE_LBF,
  prior_variance_api_default = 50,
  prior_variance_activity = "inactive_when_N_is_supplied_in_susieR_0.14.2",
  back_calculate_lbf_activity = "formal_accepted_but_unreferenced_in_coloc_5.2.3_body",
  model = "SuSiE-RSS_then_coloc.susie",
  implementation_correction = "final_pre_primary_result_conformance_lock"
), file.path(out_dir, "model_implementation_policy.csv"), row.names = FALSE)

stable_error_code <- function(reason) {
  if (grepl("version-pinned coloc and susieR packages are required", reason, fixed = TRUE)) {
    return("E_R_PACKAGES_MISSING")
  }
  if (grepl("package version mismatch", reason, fixed = TRUE)) {
    return("E_R_PACKAGE_VERSION_MISMATCH")
  }
  if (grepl("did not converge", reason, fixed = TRUE)) return("E_SUSIE_NONCONVERGENCE")
  if (grepl("orders are not identical", reason, fixed = TRUE)) return("E_LD_VARIANT_ORDER_MISMATCH")
  if (grepl("not numerically PSD", reason, fixed = TRUE)) return("E_LD_NOT_PSD")
  if (grepl("LD matrix is asymmetric", reason, fixed = TRUE)) return("E_LD_ASYMMETRIC")
  if (grepl("LD matrix diagonal is not one", reason, fixed = TRUE)) return("E_LD_DIAGONAL_INVALID")
  if (grepl("sample metadata are invalid", reason, fixed = TRUE)) return("E_SAMPLE_METADATA_INVALID")
  "E_MODEL_RUNNER"
}

write_status <- function(status, code, reason = "") {
  write.csv(data.frame(
    status = status, stable_code = code, reason = reason,
    coloc_version = actual_coloc_version,
    susieR_version = actual_susier_version,
    susie_L = SUSIE_L,
    seed = IMPLEMENTATION_SEED,
    max_iter = SUSIE_MAX_ITER,
    repeat_until_convergence = SUSIE_REPEAT_UNTIL_CONVERGENCE,
    estimate_residual_variance = SUSIE_ESTIMATE_RESIDUAL_VARIANCE,
    tol = SUSIE_TOL,
    coverage = SUSIE_COVERAGE,
    min_abs_corr = SUSIE_MIN_ABS_CORR,
    scaled_prior_variance = SUSIE_SCALED_PRIOR_VARIANCE,
    check_prior = SUSIE_CHECK_PRIOR,
    z_ld_weight = SUSIE_Z_LD_WEIGHT,
    coloc_p1 = COLOC_P1,
    coloc_p2 = COLOC_P2,
    coloc_p12 = COLOC_P12,
    coloc_overlap_min = COLOC_OVERLAP_MIN,
    coloc_trim_by_posterior = COLOC_TRIM_BY_POSTERIOR,
    coloc_back_calculate_lbf = COLOC_BACK_CALCULATE_LBF,
    prior_variance_api_default = 50,
    prior_variance_activity = "inactive_when_N_is_supplied_in_susieR_0.14.2",
    back_calculate_lbf_activity = "formal_accepted_but_unreferenced_in_coloc_5.2.3_body"
  ),
            file.path(out_dir, "status.csv"), row.names = FALSE)
}

tryCatch({
  if (!requireNamespace("coloc", quietly = TRUE) || !requireNamespace("susieR", quietly = TRUE)) {
    stop("version-pinned coloc and susieR packages are required")
  }
  if (!identical(actual_coloc_version, EXPECTED_COLOC_VERSION) ||
      !identical(actual_susier_version, EXPECTED_SUSIER_VERSION)) {
    stop(sprintf(
      "package version mismatch: expected coloc %s and susieR %s; observed coloc %s and susieR %s",
      EXPECTED_COLOC_VERSION, EXPECTED_SUSIER_VERSION,
      actual_coloc_version, actual_susier_version
    ))
  }
  summary_data <- read.csv(file.path(input_dir, "summary.csv"), check.names = FALSE)
  ld <- as.matrix(read.csv(file.path(input_dir, "ld.csv"), row.names = 1, check.names = FALSE))
  variant_ids <- as.character(summary_data$target_variant)
  if (!identical(rownames(ld), variant_ids) || !identical(colnames(ld), variant_ids)) {
    stop("summary and LD variant orders are not identical")
  }
  if (any(!is.finite(ld))) stop("LD matrix contains non-finite values")
  symmetry_error <- max(abs(ld - t(ld)))
  diagonal_error <- max(abs(diag(ld) - 1))
  if (symmetry_error > 1e-8) stop("LD matrix is asymmetric")
  if (diagonal_error > 1e-8) stop("LD matrix diagonal is not one")

  eigenvalues <- eigen(ld, symmetric = TRUE, only.values = TRUE)$values
  largest_abs <- max(abs(eigenvalues))
  rank_tolerance <- nrow(ld) * .Machine$double.eps * largest_abs
  numerical_rank <- sum(abs(eigenvalues) > rank_tolerance)
  rank_deficient <- numerical_rank < nrow(ld)
  min_eigenvalue <- min(eigenvalues)
  if (min_eigenvalue < -rank_tolerance) {
    stop(sprintf("LD is not numerically PSD: minimum eigenvalue %.9g, tolerance %.9g",
                 min_eigenvalue, rank_tolerance))
  }
  condition_number <- if (rank_deficient) Inf else largest_abs / min(abs(eigenvalues))
  write.csv(data.frame(
    n_variants = nrow(ld), numerical_rank = numerical_rank,
    rank_tolerance = rank_tolerance, rank_deficient = rank_deficient,
    min_eigenvalue = min_eigenvalue, condition_number_2 = condition_number,
    symmetry_error = symmetry_error, diagonal_error = diagonal_error,
    rank_policy = "diagnostic", conditioning_applied = FALSE,
    jitter_applied = FALSE
  ), file.path(out_dir, "ld_model_policy.csv"), row.names = FALSE)

  eqtl_n <- median(as.numeric(summary_data$eqtl_n), na.rm = TRUE)
  if (!is.finite(eqtl_n) && "an" %in% names(summary_data)) {
    eqtl_n <- median(as.numeric(summary_data$an), na.rm = TRUE) / 2
  }
  gwas_n <- median(as.numeric(summary_data$gwas_n), na.rm = TRUE)
  if (!is.finite(gwas_n) && "N" %in% names(summary_data)) {
    gwas_n <- median(as.numeric(summary_data$N), na.rm = TRUE)
  }
  beta1 <- if ("eqtl_beta" %in% names(summary_data)) summary_data$eqtl_beta else summary_data$beta
  se1 <- if ("eqtl_se" %in% names(summary_data)) summary_data$eqtl_se else summary_data$se_eqtl
  beta2 <- if ("gwas_beta" %in% names(summary_data)) summary_data$gwas_beta else summary_data$gwas_beta_aligned
  se2 <- if ("gwas_se" %in% names(summary_data)) summary_data$gwas_se else summary_data$se_gwas
  position <- if ("position" %in% names(summary_data)) summary_data$position else summary_data$BP
  maf1 <- if ("eqtl_maf" %in% names(summary_data)) summary_data$eqtl_maf else summary_data$maf
  maf2 <- if ("gwas_maf" %in% names(summary_data)) summary_data$gwas_maf else summary_data$maf
  if (!all(is.finite(c(eqtl_n, gwas_n))) || !is.finite(case_fraction) ||
      case_fraction <= 0 || case_fraction >= 1) stop("sample metadata are invalid")

  dataset1 <- list(beta = as.numeric(beta1), varbeta = as.numeric(se1)^2,
                   snp = variant_ids, position = as.numeric(position), type = "quant",
                   N = eqtl_n, MAF = as.numeric(maf1), LD = ld)
  dataset2 <- list(beta = as.numeric(beta2), varbeta = as.numeric(se2)^2,
                   snp = variant_ids, position = as.numeric(position), type = "cc",
                   s = case_fraction, N = gwas_n, MAF = as.numeric(maf2), LD = ld)
  result <- coloc::coloc.susie(
    dataset1, dataset2,
    back_calculate_lbf = COLOC_BACK_CALCULATE_LBF,
    p1 = COLOC_P1,
    p2 = COLOC_P2,
    p12 = COLOC_P12,
    overlap.min = COLOC_OVERLAP_MIN,
    trim_by_posterior = COLOC_TRIM_BY_POSTERIOR,
    susie.args = list(
      L = SUSIE_L,
      max_iter = SUSIE_MAX_ITER,
      repeat_until_convergence = SUSIE_REPEAT_UNTIL_CONVERGENCE,
      estimate_residual_variance = SUSIE_ESTIMATE_RESIDUAL_VARIANCE,
      tol = SUSIE_TOL,
      coverage = SUSIE_COVERAGE,
      min_abs_corr = SUSIE_MIN_ABS_CORR,
      scaled_prior_variance = SUSIE_SCALED_PRIOR_VARIANCE,
      check_prior = SUSIE_CHECK_PRIOR,
      z_ld_weight = SUSIE_Z_LD_WEIGHT
    )
  )
  saveRDS(result, file.path(out_dir, "coloc_susie_result.rds"))
  write.csv(if (is.null(result$summary)) data.frame() else result$summary,
            file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(if (is.null(result$results)) data.frame() else result$results,
            file.path(out_dir, "pairwise_results.csv"), row.names = FALSE)
  if (is.null(result$summary) || NROW(result$summary) == 0) {
    write_status("completed_no_comparable_credible_sets", "OK_NO_COMPARABLE_CREDIBLE_SETS",
                 "no pair of SuSiE credible sets was available")
  } else {
    write_status("completed", "OK")
  }
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
}, error = function(error) {
  reason <- conditionMessage(error)
  write_status("failed", stable_error_code(reason), reason)
  quit(status = 1)
})
