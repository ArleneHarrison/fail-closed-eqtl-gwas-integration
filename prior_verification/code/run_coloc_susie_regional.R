args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: run_coloc_susie_regional.R <input_dir> <case_fraction> <out_dir>")
}

input_dir <- args[[1]]
case_fraction <- as.numeric(args[[2]])
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_status <- function(status, reason = "") {
  write.csv(data.frame(status = status, reason = reason),
            file.path(out_dir, "status.csv"), row.names = FALSE)
}

tryCatch({
  summary_data <- read.csv(file.path(input_dir, "summary.csv"), check.names = FALSE)
  ld <- as.matrix(read.csv(file.path(input_dir, "ld.csv"), row.names = 1, check.names = FALSE))
  variant_ids <- as.character(summary_data$target_variant)
  if (!identical(rownames(ld), variant_ids) || !identical(colnames(ld), variant_ids)) {
    stop("summary and LD variant orders are not identical")
  }
  if (any(!is.finite(ld)) || max(abs(ld - t(ld))) > 1e-8) {
    stop("LD matrix contains non-finite or asymmetric values")
  }
  eigenvalues <- eigen(ld, symmetric = TRUE, only.values = TRUE)$values
  eigen_min <- min(eigenvalues)
  condition_number <- abs(max(eigenvalues)) / max(abs(eigen_min), .Machine$double.xmin)
  if (eigen_min <= -1e-6 || !is.finite(condition_number) || condition_number > 1e12) {
    stop(sprintf("LD matrix fails the fixed no-repair preflight (minimum eigenvalue %.6g; condition number %.6g)", eigen_min, condition_number))
  }
  write.csv(data.frame(min_eigenvalue_before = eigen_min, condition_number_abs_eigen = condition_number,
                       conditioning_applied = FALSE),
            file.path(out_dir, "ld_numeric_preflight.csv"), row.names = FALSE)
  eqtl_n <- median(summary_data$an) / 2
  gwas_n <- median(summary_data$N)
  dataset1 <- list(
    beta = summary_data$beta,
    varbeta = summary_data$se_eqtl ^ 2,
    snp = variant_ids,
    position = summary_data$BP,
    type = "quant",
    N = eqtl_n,
    MAF = summary_data$maf,
    LD = ld
  )
  dataset2 <- list(
    beta = summary_data$gwas_beta_aligned,
    varbeta = summary_data$se_gwas ^ 2,
    snp = variant_ids,
    position = summary_data$BP,
    type = "cc",
    s = case_fraction,
    N = gwas_n,
    MAF = summary_data$maf,
    LD = ld
  )
  result <- coloc::coloc.susie(
    dataset1, dataset2,
    susie.args = list(estimate_residual_variance = FALSE, max_iter = 1000)
  )
  saveRDS(result, file.path(out_dir, "coloc_susie_result.rds"))
  if (is.null(result$summary)) {
    write.csv(data.frame(), file.path(out_dir, "summary.csv"), row.names = FALSE)
  } else {
    write.csv(result$summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
  }
  if (is.null(result$results)) {
    write.csv(data.frame(), file.path(out_dir, "pairwise_results.csv"), row.names = FALSE)
  } else {
    write.csv(result$results, file.path(out_dir, "pairwise_results.csv"), row.names = FALSE)
  }
  if (is.null(result$summary) || NROW(result$summary) == 0) {
    write_status("completed_no_comparable_credible_sets",
                 "no pair of SuSiE credible sets was available for colocalization")
  } else {
    write_status("completed")
  }
}, error = function(error) {
  write_status("failed", conditionMessage(error))
  quit(status = 1)
})
