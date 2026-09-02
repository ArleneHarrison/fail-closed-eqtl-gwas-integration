#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(susieR)
  library(coloc)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_downstream_input_stress_test.R <ld-dir> <out-dir>")
ld_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_ld <- function(population) {
  dat <- read.delim(file.path(ld_dir, paste0("common_ld_", population, ".tsv")), check.names = FALSE)
  variants <- dat[[1]]
  matrix <- as.matrix(dat[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- variants
  colnames(matrix) <- variants
  list(variants = variants, matrix = matrix)
}

nearest_psd <- function(x, floor = 1e-8) {
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  y <- eig$vectors %*% (vals * t(eig$vectors))
  d <- sqrt(diag(y))
  y <- y / tcrossprod(d)
  y <- (y + t(y)) / 2
  diag(y) <- 1
  y
}

draw_correlated <- function(r) {
  eig <- eigen((r + t(r)) / 2, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  as.vector(eig$vectors %*% (sqrt(vals) * rnorm(length(vals))))
}

fit_susie <- function(z, r, n) {
  susie_rss(z = z, R = r, n = n, L = 3, estimate_residual_variance = FALSE,
            max_iter = 200, tol = 1e-3, verbose = FALSE)
}

cs_contains <- function(fit, causal) {
  cs <- susie_get_cs(fit, coverage = 0.95, min_abs_corr = 0)
  if (length(cs$cs) == 0L) return(FALSE)
  any(vapply(cs$cs, function(x) causal %in% x, logical(1)))
}

regional_pph4 <- function(z1, z2, n) {
  p <- length(z1)
  snp <- paste0("v", seq_len(p))
  make_dataset <- function(z) list(
    beta = z / sqrt(n), varbeta = rep(1 / n, p), type = "quant",
    N = n, MAF = rep(0.25, p), snp = snp, sdY = 1
  )
  result <- tryCatch(
    {
      invisible(capture.output(
        answer <- coloc.abf(make_dataset(z1), make_dataset(z2),
                            p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
      ))
      answer
    },
    error = function(e) NULL
  )
  if (is.null(result) || is.null(result$summary)) return(NA_real_)
  as.numeric(result$summary[["PP.H4.abf"]])
}

eur <- read_ld("EUR")
afr <- read_ld("AFR")
eas <- read_ld("EAS")
stopifnot(identical(eur$variants, afr$variants), identical(eur$variants, eas$variants))

# Matrices remain labelled by their source population.  Numerical projection is
# used only to make covariance simulation and consumer fitting stable; the raw
# matrices and their diagnostics remain archived by the Python benchmark.
r_eur <- nearest_psd(eur$matrix)
r_afr <- nearest_psd(afr$matrix)
r_eas <- nearest_psd(eas$matrix)
p <- nrow(r_eur)
n <- 5000
replicates <- 100L

conditions <- c("correct", "summary_order_permuted", "ten_percent_sign_flipped",
                "EAS_LD_substituted", "AFR_LD_substituted")
run_replicate <- function(replicate) {
  set.seed(20260831 + replicate)
  rows <- vector("list", length(conditions))
  row_i <- 1L
  causal <- 1L + ((replicate * 37L - 1L) %% p)
  effect1 <- 0.12
  effect2 <- 0.10
  z1 <- sqrt(n) * effect1 * r_eur[, causal] + draw_correlated(r_eur)
  z2 <- sqrt(n) * effect2 * r_eur[, causal] + draw_correlated(r_eur)
  fit2_correct <- fit_susie(z2, r_eur, n)
  correct_pip <- fit2_correct$pip
  correct_lead <- which.max(correct_pip)

  permutation <- ((seq_len(p) * 53L - 1L) %% p) + 1L
  flip_index <- seq(1L + (replicate %% 10L), p, by = 10L)
  z2_flipped <- z2
  z2_flipped[flip_index] <- -z2_flipped[flip_index]
  condition_inputs <- list(
    correct = list(z = z2, r = r_eur, detectable = FALSE),
    summary_order_permuted = list(z = z2[permutation], r = r_eur, detectable = TRUE),
    ten_percent_sign_flipped = list(z = z2_flipped, r = r_eur, detectable = FALSE),
    # Output drift under ancestry substitution is measured here, but protection
    # is conditional on truthful caller-supplied provenance plus an enforced
    # ancestry-policy match; the matrix contents do not authenticate ancestry.
    EAS_LD_substituted = list(z = z2, r = r_eas, detectable = NA),
    AFR_LD_substituted = list(z = z2, r = r_afr, detectable = NA)
  )

  for (condition in conditions) {
    item <- condition_inputs[[condition]]
    fit2 <- tryCatch(fit_susie(item$z, item$r, n), error = function(e) e)
    if (inherits(fit2, "error")) {
      rows[[row_i]] <- data.frame(
        replicate = replicate, condition = condition, causal_index = causal,
        consumer_completed = FALSE, contract_detectable = item$detectable,
        causal_pip = NA_real_, causal_in_95_cs = NA, lead_changed = NA,
        pip_l1_drift = NA_real_, max_pp_h4 = NA_real_, stringsAsFactors = FALSE
      )
    } else {
      rows[[row_i]] <- data.frame(
        replicate = replicate, condition = condition, causal_index = causal,
        consumer_completed = TRUE, contract_detectable = item$detectable,
        causal_pip = fit2$pip[causal], causal_in_95_cs = cs_contains(fit2, causal),
        lead_changed = which.max(fit2$pip) != correct_lead,
        pip_l1_drift = sum(abs(fit2$pip - correct_pip)),
        max_pp_h4 = regional_pph4(z1, item$z, n), stringsAsFactors = FALSE
      )
    }
    row_i <- row_i + 1L
  }
  do.call(rbind, rows)
}

worker_count <- min(10L, parallel::detectCores(logical = FALSE))
results <- do.call(rbind, parallel::mclapply(seq_len(replicates), run_replicate,
                                             mc.cores = worker_count,
                                             mc.preschedule = TRUE))
stopifnot(nrow(results) == replicates * length(conditions))
stopifnot(all(results$consumer_completed))
stopifnot(all(is.finite(results$causal_pip)))
stopifnot(all(is.finite(results$pip_l1_drift)))
stopifnot(all(is.finite(results$max_pp_h4)))
write.table(results, file.path(out_dir, "downstream_stress_test_replicates.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

summary <- do.call(rbind, lapply(split(results, results$condition), function(x) {
  data.frame(
    condition = x$condition[[1]], replicates = nrow(x), completed = sum(x$consumer_completed),
    contract_detectable = x$contract_detectable[[1]],
    median_causal_pip = median(x$causal_pip, na.rm = TRUE),
    causal_in_95_cs_rate = mean(x$causal_in_95_cs, na.rm = TRUE),
    lead_changed_rate_vs_correct = mean(x$lead_changed, na.rm = TRUE),
    median_pip_l1_drift_vs_correct = median(x$pip_l1_drift, na.rm = TRUE),
    median_max_pp_h4 = median(x$max_pp_h4, na.rm = TRUE),
    nonmissing_coloc_replicates = sum(is.finite(x$max_pp_h4)), stringsAsFactors = FALSE
  )
}))
summary <- summary[match(conditions, summary$condition), ]
write.table(summary, file.path(out_dir, "downstream_stress_test_summary.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

labels <- c(correct = "Correct contract", summary_order_permuted = "Summary order\npermuted",
            ten_percent_sign_flipped = "10% sign flips\n(content corruption)",
            EAS_LD_substituted = "EAS LD\nsubstituted", AFR_LD_substituted = "AFR LD\nsubstituted")
plot_data <- results
plot_data$condition <- factor(plot_data$condition, levels = conditions, labels = labels)

panel_a <- ggplot(plot_data, aes(condition, pip_l1_drift)) +
  geom_boxplot(width = 0.65, outlier.size = 0.7, fill = "#BFD7EA", colour = "#17324D") +
  labs(x = NULL, y = "PIP L1 drift versus correct input", tag = "a")
panel_b <- ggplot(plot_data, aes(condition, causal_pip)) +
  geom_boxplot(width = 0.65, outlier.size = 0.7, fill = "#D9EAD3", colour = "#254117") +
  labs(x = NULL, y = "PIP assigned to simulated causal variant", tag = "b")
panel_c <- ggplot(plot_data, aes(condition, max_pp_h4)) +
  geom_boxplot(width = 0.65, outlier.size = 0.7, fill = "#FCE5CD", colour = "#6B3A1E") +
  labs(x = NULL, y = "Regional PP.H4 (coloc.abf)", tag = "c")
theme_journal <- theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.tag = element_text(face = "bold"))
panel_a <- panel_a + theme_journal
panel_b <- panel_b + theme_journal
panel_c <- panel_c + theme_journal

if (requireNamespace("patchwork", quietly = TRUE)) {
  figure <- panel_a / panel_b / panel_c
  ggsave(file.path(out_dir, "figure_6_downstream_stress_test.png"), figure,
         width = 7.1, height = 9.0, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "figure_6_downstream_stress_test.svg"), figure,
         width = 7.1, height = 9.0, units = "in", bg = "white")
} else {
  ggsave(file.path(out_dir, "figure_6_downstream_stress_test.png"), panel_a,
         width = 7.1, height = 3.2, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "figure_6_downstream_stress_test.svg"), panel_a,
         width = 7.1, height = 3.2, units = "in", bg = "white")
}

manifest <- list(
  seed_rule = "20260831 + replicate index", replicates_per_condition = replicates, variants = p,
  parallel_workers = worker_count,
  nominal_sample_size = n, simulated_truth = "one shared causal variant; EUR LD data-generating matrix",
  consumers = c(paste0("susieR ", as.character(packageVersion("susieR"))),
                paste0("coloc.abf ", as.character(packageVersion("coloc")))),
  coloc_priors = "p1=1e-4, p2=1e-4, p12=1e-5; quantitative traits; sdY=1; synthetic MAF=0.25",
  important_boundary = paste(
    "The implemented exact-order check can stop the summary-order mismatch.",
    "Ancestry/LD substitution is preventable only with truthful declared provenance and an enforced caller policy; matrix contents do not authenticate ancestry.",
    "An internally corrupted effect sign with otherwise self-consistent metadata is intentionally treated",
    "as an undetectable content-level failure boundary, not as a prevented error."
  )
)
dput(manifest, file = file.path(out_dir, "downstream_stress_test_manifest.R"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
cat(sprintf("completed %d consumer fits across %d conditions\n", nrow(results), length(conditions)))
