#!/usr/bin/env Rscript

# Audit external MR comparators without silently discarding sparse candidates.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: run_external_comparator_audit.R <variants.csv> <audit.csv> <availability.csv>")
}
variants_path <- args[[1]]
audit_path <- args[[2]]
availability_path <- args[[3]]

required_packages <- c("TwoSampleMR", "MRPRESSO", "MendelianRandomization", "cause", "BMA")
availability <- data.frame(
  package = required_packages,
  available = vapply(required_packages, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
write.csv(availability, availability_path, row.names = FALSE)

if (!availability$available[availability$package == "TwoSampleMR"]) {
  stop("TwoSampleMR is required for the comparator audit")
}

variants <- read.csv(variants_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("candidate_id", "beta", "se", "gwas_beta", "gwas_se")
missing_columns <- setdiff(required_columns, names(variants))
if (length(missing_columns)) {
  stop(paste("missing columns:", paste(missing_columns, collapse = ", ")))
}

safe_error <- function(expr) {
  tryCatch(list(value = force(expr), error = ""), error = function(e) list(value = NULL, error = conditionMessage(e)))
}

value_or <- function(df, column, default) {
  if (column %in% names(df)) df[[column]] else rep(default, nrow(df))
}

variant_label <- function(df) {
  out <- if ("rsid" %in% names(df)) as.character(df$rsid) else rep(NA_character_, nrow(df))
  fallback <- if ("variant" %in% names(df)) as.character(df$variant) else paste0("variant_", seq_len(nrow(df)))
  out[is.na(out) | !nzchar(out)] <- fallback[is.na(out) | !nzchar(out)]
  out
}

method_rows <- list()
index <- 1L
for (candidate_id in unique(variants$candidate_id)) {
  d <- variants[variants$candidate_id == candidate_id, , drop = FALSE]
  valid <- is.finite(d$beta) & is.finite(d$se) & is.finite(d$gwas_beta) & is.finite(d$gwas_se) &
    d$beta != 0 & d$se > 0 & d$gwas_se > 0
  d <- d[valid, , drop = FALSE]
  n_iv <- nrow(d)
  add <- function(method, status, reason = "", estimate = NA_real_, standard_error = NA_real_, p_value = NA_real_) {
    method_rows[[index]] <<- data.frame(
      candidate_id = candidate_id,
      method = method,
      n_iv = n_iv,
      status = status,
      reason = reason,
      estimate = estimate,
      standard_error = standard_error,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
    index <<- index + 1L
  }
  if (n_iv < 3L) {
    for (method in c("IVW (TwoSampleMR)", "MR-Egger (TwoSampleMR)", "Weighted median (TwoSampleMR)", "MR-PRESSO")) {
      add(method, "ineligible", "fewer_than_three_valid_signals")
    }
    next
  }
  dat <- data.frame(
    SNP = variant_label(d),
    beta.exposure = d$beta,
    se.exposure = d$se,
    beta.outcome = d$gwas_beta,
    se.outcome = d$gwas_se,
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    effect_allele.outcome = "A",
    other_allele.outcome = "G",
    eaf.exposure = 0.20,
    eaf.outcome = 0.20,
    pval.exposure = pmax(value_or(d, "pvalue", 1e-8), .Machine$double.xmin),
    pval.outcome = pmax(value_or(d, "gwas_p", 0.5), .Machine$double.xmin),
    exposure = "eQTL",
    outcome = "cardiovascular_outcome",
    id.exposure = "exposure",
    id.outcome = "outcome",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  requested_methods <- c("mr_ivw", "mr_weighted_median")
  if (n_iv >= 4L) requested_methods <- c(requested_methods, "mr_egger_regression")
  mr_result <- safe_error(TwoSampleMR::mr(dat, method_list = requested_methods))
  method_map <- c(
    "Inverse variance weighted" = "IVW (TwoSampleMR)",
    "MR Egger" = "MR-Egger (TwoSampleMR)",
    "Weighted median" = "Weighted median (TwoSampleMR)"
  )
  if (nzchar(mr_result$error)) {
    for (method in unname(method_map)) add(method, "failed", mr_result$error)
  } else {
    result <- mr_result$value
    for (label in names(method_map)) {
      if (identical(label, "MR Egger") && n_iv < 4L) {
        add(method_map[[label]], "ineligible", "fewer_than_four_valid_signals")
        next
      }
      record <- result[result$method == label, , drop = FALSE]
      if (!nrow(record)) {
        add(method_map[[label]], "not_returned", "method_not_returned_by_TwoSampleMR")
      } else {
        add(method_map[[label]], "completed", "", record$b[[1]], record$se[[1]], record$pval[[1]])
      }
    }
  }
  if (n_iv < 4L) {
    add("MR-PRESSO", "ineligible", "fewer_than_four_valid_signals")
  } else if (!availability$available[availability$package == "MRPRESSO"]) {
    add("MR-PRESSO", "not_available", "package_not_installed")
  } else {
    presso <- safe_error(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = dat,
      NbDistribution = 1000,
      SignifThreshold = 0.05,
      seed = 20260815
    ))
    if (nzchar(presso$error)) {
      add("MR-PRESSO", "failed", presso$error)
    } else {
      add("MR-PRESSO", "completed", "global_and_outlier_results_written_to_console_log")
    }
  }
}
write.csv(do.call(rbind, method_rows), audit_path, row.names = FALSE)
