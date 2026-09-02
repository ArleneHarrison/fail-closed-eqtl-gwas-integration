args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: test_runner_minimal.R <runner>")
runner <- normalizePath(args[[1]], mustWork = TRUE)
if (!identical(as.character(packageVersion("coloc")), "5.2.3") ||
    !identical(as.character(packageVersion("susieR")), "0.14.2")) {
  stop("minimal runner test requires the locked coloc 5.2.3 and susieR 0.14.2 environment")
}

root <- tempfile("runner_minimal_")
input <- file.path(root, "input")
output <- file.path(root, "output")
dir.create(input, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

variants <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")
summary <- data.frame(
  target_variant = variants,
  position = c(100, 200, 300),
  eqtl_beta = c(0.20, -0.10, 0.05),
  eqtl_se = c(0.05, 0.05, 0.05),
  gwas_beta = c(0.10, -0.05, 0.02),
  gwas_se = c(0.04, 0.04, 0.04),
  eqtl_maf = c(0.20, 0.30, 0.40),
  gwas_maf = c(0.25, 0.35, 0.45),
  eqtl_n = rep(200, 3),
  gwas_n = rep(1000, 3)
)
write.csv(summary, file.path(input, "summary.csv"), row.names = FALSE)
ld <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), nrow = 3, byrow = TRUE)
rownames(ld) <- variants
colnames(ld) <- variants
write.csv(ld, file.path(input, "ld.csv"))

status <- system2("Rscript", c(runner, input, "0.2", output))
stopifnot(status == 0)
model_status <- read.csv(file.path(output, "status.csv"), check.names = FALSE)
policy <- read.csv(file.path(output, "model_implementation_policy.csv"), check.names = FALSE)
stopifnot(model_status$stable_code[[1]] %in% c("OK", "OK_NO_COMPARABLE_CREDIBLE_SETS"))
stopifnot(model_status$coloc_version[[1]] == "5.2.3")
stopifnot(model_status$susieR_version[[1]] == "0.14.2")
stopifnot(model_status$susie_L[[1]] == 10)
stopifnot(model_status$max_iter[[1]] == 1000)
stopifnot(!model_status$repeat_until_convergence[[1]])
stopifnot(!model_status$estimate_residual_variance[[1]])
stopifnot(model_status$tol[[1]] == 0.001)
stopifnot(model_status$coverage[[1]] == 0.95)
stopifnot(model_status$min_abs_corr[[1]] == 0.5)
stopifnot(model_status$scaled_prior_variance[[1]] == 0.2)
stopifnot(model_status$check_prior[[1]])
stopifnot(model_status$z_ld_weight[[1]] == 0)
stopifnot(model_status$coloc_p1[[1]] == 1e-4)
stopifnot(model_status$coloc_p2[[1]] == 1e-4)
stopifnot(model_status$coloc_p12[[1]] == 5e-6)
stopifnot(model_status$coloc_overlap_min[[1]] == 0.5)
stopifnot(model_status$coloc_trim_by_posterior[[1]])
stopifnot(!model_status$coloc_back_calculate_lbf[[1]])
stopifnot(model_status$prior_variance_api_default[[1]] == 50)
stopifnot(model_status$prior_variance_activity[[1]] == "inactive_when_N_is_supplied_in_susieR_0.14.2")
stopifnot(model_status$back_calculate_lbf_activity[[1]] == "formal_accepted_but_unreferenced_in_coloc_5.2.3_body")
stopifnot(policy$implementation_correction[[1]] == "final_pre_primary_result_conformance_lock")
cat("PASS minimal locked runner\n")
