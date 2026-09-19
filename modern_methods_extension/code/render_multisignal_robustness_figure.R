#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: render_multisignal_robustness_figure.R <summary.tsv> <out-dir>")
summary_file <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

architectures <- c("one_shared", "two_shared", "shared_plus_specific", "distinct_only")
ld_policies <- c("EUR_matched", "EAS_substituted", "AFR_substituted")
sample_sizes <- c(2000L, 5000L, 20000L)
plot_data <- read.delim(summary_file, check.names = FALSE)
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

if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
figure <- panel_a / panel_b
ggsave(file.path(out_dir, "figure_multisignal_robustness.png"), figure,
       width = 7.2, height = 8.8, units = "in", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "figure_multisignal_robustness.svg"), figure,
       width = 7.2, height = 8.8, units = "in", bg = "white")
