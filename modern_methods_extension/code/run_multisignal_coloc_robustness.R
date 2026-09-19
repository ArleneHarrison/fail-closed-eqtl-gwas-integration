#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(susieR)
  library(coloc)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("usage: run_multisignal_coloc_robustness.R <ld-dir> <out-dir> [replicates]")
}
ld_dir <- args[[1]]
if (!dir.exists(ld_dir)) stop("LD directory does not exist: ", ld_dir)
out_dir <- args[[2]]
replicates <- if (length(args) == 3L) as.integer(args[[3]]) else 50L
if (!is.finite(replicates) || replicates < 1L) stop("replicates must be a positive integer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed_base <- 20260920L
sample_sizes <- c(2000L, 5000L, 20000L)
architectures <- c("one_shared", "two_shared", "shared_plus_specific", "distinct_only")
ld_policies <- c("EUR_matched", "EAS_substituted", "AFR_substituted")

read_ld <- function(population) {
  dat <- read.delim(file.path(ld_dir, paste0("common_ld_", population, ".tsv")),
                    check.names = FALSE)
  variants <- dat[[1]]
  matrix <- as.matrix(dat[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- variants
  colnames(matrix) <- variants
  list(variants = variants, matrix = matrix)
}

nearest_psd <- function(x, floor = 1e-8) {
  labels <- dimnames(x)
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  y <- eig$vectors %*% (vals * t(eig$vectors))
  d <- sqrt(diag(y))
  y <- y / tcrossprod(d)
  y <- (y + t(y)) / 2
  diag(y) <- 1
  dimnames(y) <- labels
  y
}

draw_correlated <- function(r) {
  eig <- eigen((r + t(r)) / 2, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  as.vector(eig$vectors %*% (sqrt(vals) * rnorm(length(vals))))
}

make_dataset <- function(z, r, n, variants) {
  list(
    beta = z / sqrt(n), varbeta = rep(1 / n, length(z)), z = z,
    type = "quant", N = n, MAF = rep(0.25, length(z)),
    snp = variants, LD = r, sdY = 1
  )
}

fit_susie_dataset <- function(d) {
  suppressMessages(runsusie(
    d,
    maxit = 1000,
    repeat_until_convergence = FALSE,
    L = 5,
    estimate_residual_variance = FALSE,
    tol = 1e-3,
    coverage = 0.95,
    min_abs_corr = 0
  ))
}

regional_pph4_abf <- function(d1, d2) {
  d1$LD <- NULL
  d2$LD <- NULL
  invisible(capture.output(
    answer <- suppressMessages(suppressWarnings(coloc.abf(
      d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
    )))
  ))
  as.numeric(answer$summary[["PP.H4.abf"]])
}

summarize_coloc_susie <- function(fit1, fit2) {
  answer <- suppressMessages(suppressWarnings(coloc.susie(
    fit1, fit2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5,
    back_calculate_lbf = FALSE
  )))
  if (is.null(answer$summary) || nrow(answer$summary) == 0L) {
    return(c(comparable_pairs = 0, max_pp_h4 = NA_real_))
  }
  c(
    comparable_pairs = nrow(answer$summary),
    max_pp_h4 = max(answer$summary$PP.H4.abf, na.rm = TRUE)
  )
}

any_cs_contains <- function(fit, causal_indices) {
  cs <- susie_get_cs(fit, coverage = 0.95, min_abs_corr = 0)
  if (is.null(cs$cs) || length(cs$cs) == 0L) return(FALSE)
  all(vapply(causal_indices, function(index) {
    any(vapply(cs$cs, function(set) index %in% set, logical(1)))
  }, logical(1)))
}

causal_configuration <- function(replicate, architecture, p) {
  c1 <- 1L + ((replicate * 37L - 1L) %% p)
  c2 <- 1L + ((c1 + floor(p / 3) - 1L) %% p)
  c3 <- 1L + ((c1 + floor(2 * p / 3) - 1L) %% p)
  if (architecture == "one_shared") {
    return(list(t1 = c1, t2 = c1, shared = c1))
  }
  if (architecture == "two_shared") {
    return(list(t1 = c(c1, c2), t2 = c(c1, c2), shared = c(c1, c2)))
  }
  if (architecture == "shared_plus_specific") {
    return(list(t1 = c(c1, c2), t2 = c(c1, c3), shared = c1))
  }
  list(t1 = c1, t2 = c2, shared = integer(0))
}

simulate_z <- function(r, n, causal, effects) {
  mean_z <- sqrt(n) * as.vector(r[, causal, drop = FALSE] %*% effects)
  mean_z + draw_correlated(r)
}

eur <- read_ld("EUR")
eas <- read_ld("EAS")
afr <- read_ld("AFR")
stopifnot(identical(eur$variants, eas$variants), identical(eur$variants, afr$variants))

r_generate <- nearest_psd(eur$matrix)
r_fit <- list(
  EUR_matched = r_generate,
  EAS_substituted = nearest_psd(eas$matrix),
  AFR_substituted = nearest_psd(afr$matrix)
)
p <- nrow(r_generate)

rows <- vector("list", replicates * length(sample_sizes) * length(architectures) * length(ld_policies))
row_i <- 1L
for (replicate in seq_len(replicates)) {
  for (n in sample_sizes) {
    for (architecture in architectures) {
      set.seed(seed_base + replicate * 1000L + match(n, sample_sizes) * 100L +
                 match(architecture, architectures))
      truth <- causal_configuration(replicate, architecture, p)
      effects1 <- rep(0.085 / sqrt(length(truth$t1)), length(truth$t1))
      effects2 <- rep(0.075 / sqrt(length(truth$t2)), length(truth$t2))
      z1 <- simulate_z(r_generate, n, truth$t1, effects1)
      z2 <- simulate_z(r_generate, n, truth$t2, effects2)

      for (ld_policy in ld_policies) {
        fit_r <- r_fit[[ld_policy]]
        d1 <- make_dataset(z1, fit_r, n, eur$variants)
        d2 <- make_dataset(z2, fit_r, n, eur$variants)
        started <- proc.time()[["elapsed"]]
        observed <- tryCatch({
          fit1 <- fit_susie_dataset(d1)
          fit2 <- fit_susie_dataset(d2)
          coloc_multi <- summarize_coloc_susie(fit1, fit2)
          list(
            completed = TRUE,
            fit1_converged = isTRUE(fit1$converged),
            fit2_converged = isTRUE(fit2$converged),
            comparable_pairs = as.integer(coloc_multi[["comparable_pairs"]]),
            max_pp_h4_susie = as.numeric(coloc_multi[["max_pp_h4"]]),
            pp_h4_abf = regional_pph4_abf(d1, d2),
            trait1_all_causal_in_cs = any_cs_contains(fit1, truth$t1),
            trait2_all_causal_in_cs = any_cs_contains(fit2, truth$t2),
            error = ""
          )
        }, error = function(e) {
          list(
            completed = FALSE, fit1_converged = NA, fit2_converged = NA,
            comparable_pairs = NA_integer_, max_pp_h4_susie = NA_real_,
            pp_h4_abf = NA_real_, trait1_all_causal_in_cs = NA,
            trait2_all_causal_in_cs = NA, error = conditionMessage(e)
          )
        })
        elapsed <- proc.time()[["elapsed"]] - started
        rows[[row_i]] <- data.frame(
          replicate = replicate,
          sample_size = n,
          architecture = architecture,
          ld_policy = ld_policy,
          n_trait1_causal = length(truth$t1),
          n_trait2_causal = length(truth$t2),
          n_shared_causal = length(truth$shared),
          completed = observed$completed,
          fit1_converged = observed$fit1_converged,
          fit2_converged = observed$fit2_converged,
          comparable_pairs = observed$comparable_pairs,
          max_pp_h4_susie = observed$max_pp_h4_susie,
          pp_h4_abf = observed$pp_h4_abf,
          trait1_all_causal_in_cs = observed$trait1_all_causal_in_cs,
          trait2_all_causal_in_cs = observed$trait2_all_causal_in_cs,
          runtime_seconds = elapsed,
          error = observed$error,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      }
    }
  }
}
results <- do.call(rbind, rows)
write.table(results, file.path(out_dir, "multisignal_replicates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

wilson_interval <- function(successes, total, z = qnorm(0.975)) {
  if (total == 0L) return(c(NA_real_, NA_real_))
  phat <- successes / total
  denom <- 1 + z^2 / total
  centre <- (phat + z^2 / (2 * total)) / denom
  half <- z * sqrt(phat * (1 - phat) / total + z^2 / (4 * total^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}

groups <- split(results, interaction(results$architecture, results$sample_size,
                                     results$ld_policy, drop = TRUE))
summary_rows <- lapply(groups, function(x) {
  complete <- x[x$completed, , drop = FALSE]
  h4_susie <- sum(complete$max_pp_h4_susie >= 0.8, na.rm = TRUE)
  h4_abf <- sum(complete$pp_h4_abf >= 0.8, na.rm = TRUE)
  denom_susie <- sum(is.finite(complete$max_pp_h4_susie))
  denom_abf <- sum(is.finite(complete$pp_h4_abf))
  ci_susie <- wilson_interval(h4_susie, denom_susie)
  ci_abf <- wilson_interval(h4_abf, denom_abf)
  data.frame(
    architecture = x$architecture[[1]],
    sample_size = x$sample_size[[1]],
    ld_policy = x$ld_policy[[1]],
    attempted = nrow(x),
    completed = sum(x$completed),
    both_susie_converged = sum(x$fit1_converged & x$fit2_converged, na.rm = TRUE),
    median_comparable_pairs = if (nrow(complete)) median(complete$comparable_pairs, na.rm = TRUE) else NA,
    median_max_pp_h4_susie = if (denom_susie) median(complete$max_pp_h4_susie, na.rm = TRUE) else NA,
    h4_ge_0_8_susie = h4_susie,
    h4_ge_0_8_susie_denominator = denom_susie,
    h4_ge_0_8_susie_rate = if (denom_susie) h4_susie / denom_susie else NA,
    h4_ge_0_8_susie_ci_low = ci_susie[[1]],
    h4_ge_0_8_susie_ci_high = ci_susie[[2]],
    median_pp_h4_abf = if (denom_abf) median(complete$pp_h4_abf, na.rm = TRUE) else NA,
    h4_ge_0_8_abf = h4_abf,
    h4_ge_0_8_abf_denominator = denom_abf,
    h4_ge_0_8_abf_rate = if (denom_abf) h4_abf / denom_abf else NA,
    h4_ge_0_8_abf_ci_low = ci_abf[[1]],
    h4_ge_0_8_abf_ci_high = ci_abf[[2]],
    trait1_all_causal_cs_rate = mean(complete$trait1_all_causal_in_cs, na.rm = TRUE),
    trait2_all_causal_cs_rate = mean(complete$trait2_all_causal_in_cs, na.rm = TRUE),
    median_runtime_seconds = median(x$runtime_seconds),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(match(summary$architecture, architectures), summary$sample_size,
                         match(summary$ld_policy, ld_policies)), ]
write.table(summary, file.path(out_dir, "multisignal_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_data <- summary
plot_data$architecture <- factor(
  plot_data$architecture, levels = architectures,
  labels = c("One shared", "Two shared", "Shared + trait-specific", "Distinct only")
)
plot_data$ld_policy <- factor(
  plot_data$ld_policy, levels = ld_policies,
  labels = c("EUR matched", "EAS substituted", "AFR substituted")
)
plot_data$sample_size <- factor(plot_data$sample_size, levels = sample_sizes,
                                labels = paste0("N = ", format(sample_sizes, big.mark = ",")))

panel_a <- ggplot(plot_data, aes(ld_policy, median_max_pp_h4_susie, colour = sample_size,
                                 group = sample_size)) +
  geom_point(size = 1.8) + geom_line(linewidth = 0.55) +
  facet_wrap(~ architecture, ncol = 2) +
  scale_colour_manual(values = c("#4C78A8", "#F58518", "#54A24B")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Median maximum PP.H4 (coloc.susie)", colour = NULL, tag = "a")

panel_b <- ggplot(plot_data, aes(ld_policy, h4_ge_0_8_susie_rate, fill = sample_size)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_errorbar(aes(ymin = h4_ge_0_8_susie_ci_low, ymax = h4_ge_0_8_susie_ci_high,
                    group = sample_size),
                position = position_dodge(width = 0.75), width = 0.18, linewidth = 0.35) +
  facet_wrap(~ architecture, ncol = 2) +
  scale_fill_manual(values = c("#9ECAE1", "#6BAED6", "#2171B5")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion with PP.H4 >= 0.8", fill = NULL, tag = "b")

theme_journal <- theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.tag = element_text(face = "bold"), legend.position = "top")
panel_a <- panel_a + theme_journal
panel_b <- panel_b + theme_journal

if (requireNamespace("patchwork", quietly = TRUE)) {
  figure <- panel_a / panel_b
  ggsave(file.path(out_dir, "figure_multisignal_robustness.png"), figure,
         width = 7.2, height = 8.8, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "figure_multisignal_robustness.svg"), figure,
         width = 7.2, height = 8.8, units = "in", bg = "white")
} else {
  ggsave(file.path(out_dir, "figure_multisignal_robustness.png"), panel_a,
         width = 7.2, height = 4.3, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "figure_multisignal_robustness.svg"), panel_a,
         width = 7.2, height = 4.3, units = "in", bg = "white")
}

manifest <- list(
  schema_version = "1.0",
  seed_base = seed_base,
  replicates = replicates,
  variants = p,
  sample_sizes = sample_sizes,
  architectures = architectures,
  ld_policies = ld_policies,
  data_generating_ld = "1000 Genomes EUR common-variant LD; nearest-PSD projection for simulation only",
  consumers = c(
    paste0("susieR ", as.character(packageVersion("susieR"))),
    paste0("coloc ", as.character(packageVersion("coloc")))
  ),
  model_settings = "SuSiE-RSS L=5, max_iter=1000, tol=1e-3; coloc priors p1=1e-4, p2=1e-4, p12=1e-5",
  interpretation_boundary = paste(
    "Synthetic robustness benchmark only. It compares consumer behavior across declared architectures,",
    "sample sizes and LD substitutions; it does not estimate empirical cardiovascular association,",
    "biological colocalisation, clinical validity or superiority over another software package."
  )
)
dput(manifest, file = file.path(out_dir, "multisignal_manifest.R"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))

cat(sprintf("completed %d/%d condition rows\n", sum(results$completed), nrow(results)))
if (any(!results$completed)) {
  cat("non-completed rows by error:\n")
  print(sort(table(results$error[!results$completed]), decreasing = TRUE))
}
