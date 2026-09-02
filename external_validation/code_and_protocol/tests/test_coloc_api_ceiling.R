if (!identical(as.character(packageVersion("coloc")), "5.2.3") ||
    !identical(as.character(packageVersion("susieR")), "0.14.2")) {
  stop("API ceiling test requires coloc 5.2.3 and susieR 0.14.2")
}

runsusie <- getFromNamespace("runsusie", "coloc")
stopifnot("repeat_until_convergence" %in% names(formals(runsusie)))
stopifnot(identical(formals(runsusie)$repeat_until_convergence, TRUE))
stopifnot(identical(formals(susieR::susie_rss)$z_ld_weight, 0))
stopifnot(identical(formals(susieR::susie_rss)$prior_variance, 50))
stopifnot(identical(formals(susieR::susie_rss)$check_prior, TRUE))
susie_suff_stat <- getFromNamespace("susie_suff_stat", "susieR")
stopifnot(identical(formals(susie_suff_stat)$tol, 0.001))
stopifnot(identical(formals(susie_suff_stat)$coverage, 0.95))
stopifnot(identical(formals(susie_suff_stat)$min_abs_corr, 0.5))
stopifnot(identical(formals(susie_suff_stat)$scaled_prior_variance, 0.2))
stopifnot(identical(formals(coloc::coloc.susie)$back_calculate_lbf, FALSE))
coloc_bf_bf <- getFromNamespace("coloc.bf_bf", "coloc")
stopifnot(identical(formals(coloc_bf_bf)$overlap.min, 0.5))
stopifnot(identical(formals(coloc_bf_bf)$trim_by_posterior, TRUE))

# Source-level version-path audit. With N supplied, susie_rss consumes its
# prior_variance formal but calls susie_suff_stat with scaled_prior_variance;
# therefore only scaled_prior_variance=0.2 is represented as effective in the
# formal lock. coloc 5.2.3 accepts back_calculate_lbf but does not reference it.
susie_rss_body <- paste(deparse(body(susieR::susie_rss)), collapse = "\n")
coloc_susie_body <- paste(deparse(body(coloc::coloc.susie)), collapse = "\n")
stopifnot(grepl("scaled_prior_variance", susie_rss_body, fixed = TRUE))
stopifnot(grepl("prior_variance", susie_rss_body, fixed = TRUE))
stopifnot(!grepl("back_calculate_lbf", coloc_susie_body, fixed = TRUE))

# Dynamic routing proof for the effective N-supplied prior. The trace observes
# the value actually received by susie_suff_stat, not just a matching token in
# a wrapper body.
.captured_suff_stat_args <- NULL
trace(
  "susie_suff_stat", where = asNamespace("susieR"), print = FALSE,
  tracer = quote(.GlobalEnv$.captured_suff_stat_args <- list(
    scaled_prior_variance = scaled_prior_variance,
    estimate_residual_variance = estimate_residual_variance,
    tol = tol,
    coverage = coverage,
    min_abs_corr = min_abs_corr
  ))
)
tryCatch(
  invisible(susieR::susie_rss(
    z = c(2, 1, -1), R = diag(3), n = 1000, L = 1L, max_iter = 1L,
    scaled_prior_variance = 0.2, estimate_residual_variance = FALSE,
    tol = 0.001, coverage = 0.95, min_abs_corr = 0.5,
    check_prior = TRUE, z_ld_weight = 0
  )),
  finally = untrace("susie_suff_stat", where = asNamespace("susieR"))
)
stopifnot(identical(.captured_suff_stat_args, list(
  scaled_prior_variance = 0.2,
  estimate_residual_variance = FALSE,
  tol = 0.001,
  coverage = 0.95,
  min_abs_corr = 0.5
)))

set.seed(1)
p <- 80
ids <- paste0("v", seq_len(p))
ld <- toeplitz(0.95 ^ (0:(p - 1)))
dimnames(ld) <- list(ids, ids)
z <- rep(c(6, -5, 4, -3), length.out = p) + rnorm(p)
dataset <- list(
  beta = z * 0.1,
  varbeta = rep(0.01, p),
  snp = ids,
  position = seq_len(p),
  type = "quant",
  N = 1000,
  MAF = rep(0.2, p),
  LD = ld
)

message <- tryCatch({
  coloc::coloc.susie(
    dataset, dataset,
    back_calculate_lbf = FALSE,
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 5e-6,
    overlap.min = 0.5,
    trim_by_posterior = TRUE,
    susie.args = list(
      L = 10L,
      max_iter = 1L,
      repeat_until_convergence = FALSE,
      estimate_residual_variance = FALSE,
      tol = 0.001,
      coverage = 0.95,
      min_abs_corr = 0.5,
      scaled_prior_variance = 0.2,
      check_prior = TRUE,
      z_ld_weight = 0
    )
  )
  "NO_ERROR"
}, error = function(error) conditionMessage(error))

# End-to-end proof that coloc.susie forwards repeat_until_convergence through
# susie.args to the coloc::runsusie wrapper: it must stop at one iteration and
# must not enter the wrapper's default x100 retry branch.
stopifnot(grepl("did not converge in 1 iterations", message, fixed = TRUE))
cat("PASS coloc 5.2.3/susieR 0.14.2 explicit API route and iteration ceiling\n")
