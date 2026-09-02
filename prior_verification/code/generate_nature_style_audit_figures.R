#!/usr/bin/env Rscript

# Generate journal-style scientific figures directly from locked audit outputs.
# The figures describe workflow eligibility only. They do not display association,
# posterior, biological, clinical, causal, or comparative-performance results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(jsonlite)
  library(digest)
  library(grid)
})

options(stringsAsFactors = FALSE)

ROOT <- gsub("\\\\", "/", Sys.getenv("HCSMR_ROOT", unset = getwd()))
ROOT <- sub("/$", "", ROOT)
OUT <- file.path(ROOT, "paper", "final_submission_package_20260826", "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

INK <- "#1A1A1A"
BLUE <- "#0072B2"
ORANGE <- "#D55E00"
TEAL <- "#009E73"
SKY <- "#56B4E9"
YELLOW <- "#E69F00"
GREY <- "#9E9E9E"
LIGHT <- "#E6E6E6"
PALE <- "#F5F5F5"

theme_journal <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(size = rel(0.88), colour = "#4D4D4D", hjust = 0),
      axis.title = element_text(face = "plain", colour = INK),
      axis.text = element_text(colour = INK),
      axis.line = element_line(linewidth = 0.35, colour = INK),
      axis.ticks = element_line(linewidth = 0.35, colour = INK),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(4.5, "mm"),
      legend.position = "top",
      legend.justification = "left",
      panel.grid.major.y = element_line(linewidth = 0.25, colour = LIGHT),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

theme_void_journal <- function(base_size = 8.5) {
  theme_void(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(size = rel(0.88), colour = "#4D4D4D", hjust = 0),
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 6)
    )
}

sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

assert_contains <- function(path, strings) {
  text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  missing <- strings[!vapply(strings, function(x) grepl(x, text, fixed = TRUE), logical(1))]
  if (length(missing)) stop("Locked evidence mismatch in ", path, ": ", paste(missing, collapse = "; "))
  invisible(text)
}

save_pair <- function(plot, number, filename, width, height, sources, caption, title) {
  png_path <- file.path(OUT, paste0(filename, ".png"))
  svg_path <- file.path(OUT, paste0(filename, ".svg"))
  if (file.exists(png_path)) unlink(png_path)
  grDevices::png(png_path, width = width, height = height, units = "in", res = 600, bg = "white", type = "cairo")
  print(plot)
  grDevices::dev.off()
  ggsave(svg_path, plot, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
  script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  script_path <- if (grepl("^[A-Za-z]:", script_arg)) script_arg else file.path(ROOT, script_arg)
  rel <- function(path) sub(paste0("^", gsub("([.\\+*?^$(){}|\\[\\]\\\\])", "\\\\\\1", ROOT), "/?"), "", gsub("\\\\", "/", path))
  payload <- list(
    schema_version = "2.0",
    figure_number = number,
    title = title,
    caption = caption,
    figure_contract = "workflow-only; no association, posterior, biological, clinical, causal, or comparative-performance claim",
    visual_style = "R/ggplot2 journal figure: white background, axis-led panels, thin black rules, Arial-compatible typography, colour-blind-safe accents, editable SVG and 600-dpi PNG",
    software = list(
      R = R.version.string,
      ggplot2 = as.character(packageVersion("ggplot2")),
      patchwork = as.character(packageVersion("patchwork")),
      svglite = as.character(packageVersion("svglite"))
    ),
    sources = lapply(sources, function(path) list(
      path = rel(path),
      sha256 = sha256(path),
      byte_count = unname(file.info(path)$size)
    )),
    script = rel(script_path),
    script_sha256 = sha256(script_path),
    outputs = list(
      png = list(path = rel(png_path), sha256 = sha256(png_path), byte_count = unname(file.info(png_path)$size)),
      svg = list(path = rel(svg_path), sha256 = sha256(svg_path), byte_count = unname(file.info(svg_path)$size))
    )
  )
  manifest_path <- file.path(OUT, paste0(if (number == "1B") "figure_1b" else paste0("figure_", tolower(number)), "_manifest.json"))
  write_json(payload, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(payload)
}

`%||%` <- function(x, y) if (length(x) == 0 || is.null(x) || is.na(x)) y else x

# Locked inputs.
gate_spec <- file.path(ROOT, "paper", "new_manuscript_preflight_gate_spec_2026-08-26.md")
publication_tsv <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260826", "publication_table_data.tsv")
fixed_report <- file.path(ROOT, "paper", "reanalysis_outputs", "local_verification_20260823", "ifitm2_cad_minimal_regional_run_20260824.md")
coordinate_report <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260826", "coordinate_control", "benchmark_run_report_20260826.md")
mutation_tsv <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "preflight_mutation_benchmark", "mutation_summary.tsv")
boundary_tsv <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "preflight_mutation_benchmark", "tolerance_boundary.tsv")
scaling_tsv <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "preflight_mutation_benchmark", "scaling_summary.tsv")
transfer_tsv <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "multiregion_transfer", "executed", "benchmark_summary.tsv")
transfer_json <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "multiregion_transfer", "executed", "benchmark_summary.json")
transfer_verification <- file.path(ROOT, "paper", "reanalysis_outputs", "new_manuscript_20260828", "multiregion_transfer", "executed", "input_verification.json")

assert_contains(fixed_report, c("295,093 indexed source rows", "5,182 IFITM2 rows", "5,047 mapped and 135 unmapped", "4,657 eQTL--CAD rows", "2,863 aligned and 1,794 flipped", "4,365 unique harmonized variants", "condition number `7.63e14`"))
assert_contains(coordinate_report, c("251 retained selected-identifier rows", "133 aligned, 110 explicitly reversed", "243 rows entered the summary audit", "225 were present in the locked LD archive", "condition number 810317569661433.5"))

# Figure 1: minimal workflow architecture.
stages <- data.frame(
  x = 1:7,
  label = c("Source\nidentity", "Coordinate\nmapping", "Allele\naudit", "VCF\ncoverage", "Trait\nmetadata", "LD identity\n& order", "Numerical\ngate")
)
p1 <- ggplot(stages, aes(x, y = 1)) +
  geom_segment(data = data.frame(x = 1:6, xend = 2:7), aes(x = x + 0.34, xend = xend - 0.34, y = 1, yend = 1), linewidth = 0.45, colour = INK, arrow = arrow(length = unit(1.6, "mm"), type = "closed"), inherit.aes = FALSE) +
  geom_tile(width = 0.66, height = 0.42, fill = "white", colour = INK, linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.6, lineheight = 0.95, colour = INK) +
  annotate("segment", x = 7.34, xend = 7.78, y = 1, yend = 1.23, linewidth = 0.45, colour = INK, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("segment", x = 7.34, xend = 7.78, y = 1, yend = 0.77, linewidth = 0.45, colour = INK, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("rect", xmin = 7.78, xmax = 8.75, ymin = 1.08, ymax = 1.38, fill = alpha(TEAL, 0.10), colour = TEAL, linewidth = 0.55) +
  annotate("text", x = 8.265, y = 1.23, label = "READY\nmodel eligible", size = 2.3, fontface = "bold", colour = TEAL, lineheight = 0.95) +
  annotate("rect", xmin = 7.78, xmax = 8.75, ymin = 0.62, ymax = 0.92, fill = alpha(ORANGE, 0.10), colour = ORANGE, linewidth = 0.55) +
  annotate("text", x = 8.265, y = 0.77, label = "INELIGIBLE\nstop + audit", size = 2.3, fontface = "bold", colour = ORANGE, lineheight = 0.95) +
  annotate("text", x = 4.5, y = 0.46, label = "Raw and transformed fields, row-level status, source hashes, and terminal code are retained at every stage", size = 2.45, colour = "#4D4D4D") +
  coord_cartesian(xlim = c(0.55, 8.82), ylim = c(0.35, 1.5), clip = "off") +
  theme_void_journal(8.5)
save_pair(p1, "1", "figure_1", 7.2, 2.35, c(gate_spec), "Figure 1. Fail-closed workflow from immutable source record to technical gate. Workflow architecture only; a READY status is not an association or model result, and an INELIGIBLE status preserves a stop record.", "Fail-closed workflow architecture")

# Figure 1B: executable outcomes represented as scientific panels rather than cards.
pub <- read.delim(publication_tsv, check.names = FALSE)
pub$case_label <- factor(pub$case_id, levels = rev(pub$case_id), labels = rev(c("Synthetic READY", "Synthetic numerical stop", "Fixed processing case", "Coordinate control")))
pub$status <- factor(pub$status, levels = c("READY", "INELIGIBLE"))

p1b_a <- ggplot(pub, aes(x = 1, y = case_label, colour = status)) +
  geom_point(size = 3.0, shape = 16) +
  geom_text(aes(x = 1.08, label = status), hjust = 0, size = 2.45, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c(READY = TEAL, INELIGIBLE = ORANGE), guide = "none") +
  coord_cartesian(xlim = c(0.95, 1.48)) +
  labs(x = NULL, y = NULL, title = "Decision") +
  theme_journal() + theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())

counts <- pub %>% select(case_label, summary_rows, ld_variant_ids) %>% pivot_longer(-case_label, names_to = "measure", values_to = "count")
p1b_b <- ggplot(counts, aes(x = count, y = case_label, shape = measure)) +
  geom_point(size = 2.2, colour = BLUE) +
  scale_x_log10(labels = label_number(big.mark = ","), breaks = c(1, 10, 100, 1000)) +
  scale_shape_manual(values = c(summary_rows = 16, ld_variant_ids = 1), labels = c(summary_rows = "Summary rows", ld_variant_ids = "LD variants")) +
  labs(x = "Count (log scale)", y = NULL, shape = NULL, title = "Input size") +
  theme_journal() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

cond <- pub %>% mutate(condition_numeric = suppressWarnings(as.numeric(gsub(" .*", "", condition_number))), log_condition = log10(condition_numeric))
p1b_c <- ggplot(cond, aes(x = log_condition, y = case_label, colour = status)) +
  geom_vline(xintercept = 12, linetype = 2, linewidth = 0.45, colour = INK) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(READY = TEAL, INELIGIBLE = ORANGE), guide = "none") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 15), labels = c("1", expression(10^4), expression(10^8), expression(10^12), expression(10^15))) +
  labs(x = "Condition number", y = NULL, title = "LD preflight") +
  theme_journal() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

p1b <- (p1b_a | p1b_b | p1b_c) + plot_layout(widths = c(1.1, 1.25, 1.25), guides = "collect") + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial"), legend.position = "bottom"))
save_pair(p1b, "1B", "figure_1b_executable_gate_outcomes", 7.2, 3.15, c(publication_tsv), "Figure 1B. Executable gate outcomes rendered solely from the locked publication table. READY and INELIGIBLE are technical workflow states; neither is an association, posterior, biological result, calibration estimate, or performance comparison.", "Executable gate outcomes")

# Figure 2: fixed-case accounting.
flow <- data.frame(
  step = factor(c("Regional source", "Selected rows", "Mapped", "Traceable allele", "Unique LD"), levels = rev(c("Regional source", "Selected rows", "Mapped", "Traceable allele", "Unique LD"))),
  count = c(295093, 5182, 5047, 4657, 4365)
)
p2_a <- ggplot(flow, aes(x = count, y = step)) +
  geom_segment(aes(x = 1, xend = count, yend = step), linewidth = 0.45, colour = GREY) +
  geom_point(size = 2.6, colour = BLUE) +
  geom_text(aes(label = comma(count)), hjust = 1.18, nudge_x = 0, size = 2.5, colour = INK) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5), labels = c("1k", "10k", "100k"), expand = expansion(mult = c(0.02, 0.30))) +
  labs(x = "Rows or variants (log scale)", y = NULL, title = "Processing cascade") +
  theme_journal() + theme(panel.grid.major.y = element_blank())

recon <- data.frame(
  category = factor(c("Aligned", "Explicitly flipped", "No exact allele match", "GWAS coordinate absent", "Unmapped liftover", "VCF absent", "VCF allele mismatch"), levels = rev(c("Aligned", "Explicitly flipped", "No exact allele match", "GWAS coordinate absent", "Unmapped liftover", "VCF absent", "VCF allele mismatch"))),
  count = c(2863, 1794, 1, 96, 135, 171, 122),
  class = c("Traceable", "Traceable", rep("Retained failure", 5))
)
p2_b <- ggplot(recon, aes(x = count, y = category, fill = class)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = comma(count)), hjust = -0.15, size = 2.45) +
  scale_fill_manual(values = c(Traceable = BLUE, `Retained failure` = ORANGE), guide = guide_legend(title = NULL)) +
  scale_x_continuous(limits = c(0, 3500), expand = expansion(mult = c(0, 0)), labels = label_number(big.mark = ",")) +
  labs(x = "Selected eQTL rows", y = NULL, title = "Reconciliation of 5,182 rows") +
  theme_journal() + theme(panel.grid.major.y = element_blank())

vcf <- data.frame(class = factor(c("Eligible", "Retained ineligible"), levels = c("Retained ineligible", "Eligible")), count = c(42192, 222))
p2_c <- ggplot(vcf, aes(x = "VCF", y = count, fill = class)) +
  geom_col(width = 0.48) +
  annotate("text", x = 1, y = 21500, label = "42,192\neligible", size = 2.45, colour = "white", fontface = "bold", lineheight = 0.95) +
  annotate("text", x = 1, y = 43700, label = "222 retained\nineligible", size = 2.15, colour = ORANGE, fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = c(Eligible = TEAL, `Retained ineligible` = ORANGE), guide = "none") +
  scale_y_continuous(limits = c(0, 45500), labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Regional VCF records", title = "VCF audit") +
  theme_journal() + theme(panel.grid.major.x = element_blank())

p2 <- (p2_a | p2_b | p2_c) + plot_layout(widths = c(1.05, 1.35, 0.72), guides = "collect") + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial"), legend.position = "bottom"))
save_pair(p2, "2", "figure_2", 7.2, 3.75, c(fixed_report), "Figure 2. Fixed-case full-row accounting from the locked audit record. All failed and ambiguous rows remain explicit. Counts terminate at the technical preflight; no model, association, posterior, biological, or treatment result is shown.", "Fixed-case full-row audit")

# Figure 3: normalized numerical stop diagnostic.
diag_fixed <- data.frame(
  metric = factor(c("Symmetry error", "Diagonal error", "Condition number"), levels = rev(c("Symmetry error", "Diagonal error", "Condition number"))),
  observed = c(2.22e-16, 2.22e-16, 7.63e14),
  limit = c(1e-8, 1e-8, 1e12),
  decision = c("within limit", "within limit", "stop")
) %>% mutate(ratio = observed / limit)
p3 <- ggplot(diag_fixed, aes(x = ratio, y = metric, colour = decision)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5, colour = INK) +
  geom_segment(aes(x = pmin(ratio, 1), xend = pmax(ratio, 1), yend = metric), linewidth = 0.55) +
  geom_point(size = 3.0) +
  scale_colour_manual(values = c(`within limit` = BLUE, stop = ORANGE), guide = "none") +
  scale_x_log10(breaks = c(1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2), labels = parse(text = c("10^-8", "10^-6", "10^-4", "10^-2", "1", "10^2")), expand = expansion(mult = c(0.06, 0.18))) +
  labs(x = "Observed value / prospectively supplied limit", y = NULL) +
  annotate("text", x = 1, y = 3.48, label = "Limit", size = 2.5, fontface = "bold", vjust = 0, colour = INK) +
  annotate("text", x = 1e-4, y = 3.35, label = "within supplied limit", size = 2.35, colour = BLUE) +
  annotate("text", x = 40, y = 3.35, label = "stop", size = 2.35, colour = ORANGE) +
  annotate("label", x = 250, y = 1.15, label = "E_LD_ILL_CONDITIONED\ncondition number = 7.63 x 10^14\nno repair; no model call", hjust = 1, vjust = 0, size = 2.4, colour = ORANGE, fill = "white") +
  coord_cartesian(ylim = c(0.55, 3.5), clip = "off") +
  theme_journal() + theme(panel.grid.major.y = element_blank(), legend.position = "none")
save_pair(p3, "3", "figure_3", 7.2, 3.25, c(fixed_report, gate_spec), "Figure 3. Fixed-case numerical preflight diagnostic normalized to the supplied limits. The condition-number check crossed the prospective stop boundary and produced an INELIGIBLE record. No matrix repair or model was run.", "Fixed-case numerical stop diagnostic")

# Figure 4: independent coordinate-control benchmark.
coord_counts <- data.frame(category = factor(c("Selected rows", "Aligned", "Flipped", "Retained failure", "Unique summary", "LD archive"), levels = rev(c("Selected rows", "Aligned", "Flipped", "Retained failure", "Unique summary", "LD archive"))), count = c(251, 133, 110, 8, 234, 225), class = c("input", "traceable", "traceable", "failure", "input", "input"))
p4_a <- ggplot(coord_counts, aes(x = count, y = category, colour = class)) +
  geom_segment(aes(x = 0, xend = count, yend = category), linewidth = 0.45, colour = GREY) +
  geom_point(size = 2.5) +
  geom_text(aes(label = count), hjust = -0.18, size = 2.45, colour = INK) +
  scale_colour_manual(values = c(input = BLUE, traceable = TEAL, failure = ORANGE), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Rows or variants", y = NULL, title = "Coordinate-control accounting") +
  theme_journal() + theme(panel.grid.major.y = element_blank())

coord_codes <- data.frame(code = factor(c("Duplicate", "MAF invalid", "Case fraction", "LD set mismatch", "LD ill-conditioned"), levels = rev(c("Duplicate", "MAF invalid", "Case fraction", "LD set mismatch", "LD ill-conditioned"))), affected = c(9, 9, 243, 9, 1))
p4_b <- ggplot(coord_codes, aes(x = affected, y = code)) +
  geom_col(width = 0.58, fill = ORANGE) +
  geom_text(aes(label = affected), hjust = -0.15, size = 2.45) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Affected rows or checks", y = NULL, title = "Terminal codes") +
  theme_journal() + theme(panel.grid.major.y = element_blank())

coord_diag <- data.frame(metric = factor(c("Symmetry", "Diagonal", "Condition"), levels = rev(c("Symmetry", "Diagonal", "Condition"))), ratio = c(2.220446049250313e-16 / 1e-8, 2.220446049250313e-16 / 1e-8, 810317569661433.5 / 1e12), decision = c("within limit", "within limit", "stop"))
p4_c <- ggplot(coord_diag, aes(x = ratio, y = metric, colour = decision)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, colour = INK) +
  geom_point(size = 2.6) +
  scale_x_log10(breaks = c(1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2), labels = parse(text = c("10^-8", "10^-6", "10^-4", "10^-2", "1", "10^2"))) +
  scale_colour_manual(values = c(`within limit` = BLUE, stop = ORANGE), guide = "none") +
  labs(x = "Observed / limit", y = NULL, title = "Numerical ratios") +
  theme_journal() + theme(panel.grid.major.y = element_blank())

p4 <- (p4_a | p4_b | p4_c) + plot_layout(widths = c(1.05, 1.05, 0.90)) + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial")))
save_pair(p4, "4", "figure_4", 7.2, 3.65, c(coordinate_report, gate_spec), "Figure 4. Independent coordinate-control workflow benchmark. The panels show full-row accounting, retained technical codes, and normalized preflight diagnostics. This coordinate-only benchmark has no association, biological, causal, cell-type, or clinical interpretation.", "Independent coordinate-control benchmark")

# Figure 5: systematic contract validation, boundaries and scaling.
mut <- read.delim(mutation_tsv, check.names = FALSE)
mut$label <- gsub("_", " ", mut$case)
mut$label <- factor(mut$label, levels = rev(mut$label))
p5_a <- ggplot(mut, aes(x = detection_rate, y = label)) +
  geom_errorbarh(aes(xmin = exact_95_ci_low, xmax = exact_95_ci_high), height = 0, linewidth = 0.45, colour = BLUE) +
  geom_point(size = 1.7, colour = BLUE) +
  scale_x_continuous(limits = c(0.95, 1.005), breaks = c(0.96, 0.98, 1.00), labels = percent_format(accuracy = 1)) +
  labs(x = "Expected-code detection (95% exact CI)", y = NULL, title = "Single-fault mutation checks") +
  theme_journal(7.2) + theme(panel.grid.major.y = element_blank())

bound <- read.delim(boundary_tsv, check.names = FALSE)
bound$short <- factor(c("Asymmetry 0.5x", "Diagonal 0.5x", "Asymmetry 2x", "Diagonal 2x", "Condition 0.40x", "Condition 4.00x"), levels = rev(c("Asymmetry 0.5x", "Diagonal 0.5x", "Asymmetry 2x", "Diagonal 2x", "Condition 0.40x", "Condition 4.00x")))
bound$ratio <- c(bound$symmetry_max_error[1] / 1e-8, bound$diagonal_max_error[2] / 1e-8, bound$symmetry_max_error[3] / 1e-8, bound$diagonal_max_error[4] / 1e-8, bound$condition_number[5] / 1e12, bound$condition_number[6] / 1e12)
p5_b <- ggplot(bound, aes(x = ratio, y = short, colour = observed_status)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, colour = INK) +
  geom_point(size = 2.4) +
  scale_x_log10(breaks = c(0.4, 0.5, 1, 2, 4), labels = label_number(accuracy = 0.1)) +
  scale_colour_manual(values = c(READY = TEAL, INELIGIBLE = ORANGE), name = NULL) +
  labs(x = "Observed / boundary", y = NULL, title = "Prospective boundary behavior") +
  theme_journal(7.5) + theme(panel.grid.major.y = element_blank())

scale_df <- read.delim(scaling_tsv, check.names = FALSE)
p5_c <- ggplot(scale_df, aes(x = n_variants, y = median_runtime_ms)) +
  geom_ribbon(aes(ymin = min_runtime_ms, ymax = max_runtime_ms), fill = alpha(SKY, 0.22), colour = NA) +
  geom_line(linewidth = 0.55, colour = BLUE) +
  geom_point(size = 2.0, colour = BLUE) +
  scale_x_log10(breaks = scale_df$n_variants, labels = label_number(big.mark = ",")) +
  scale_y_log10(labels = label_number()) +
  labs(x = "Variants", y = "Gate runtime (ms; log scale)", title = "Descriptive local scaling") +
  theme_journal(7.5)

p5 <- (p5_a | (p5_b / p5_c)) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial")))
save_pair(p5, "5", "figure_5", 7.2, 6.25, c(mutation_tsv, boundary_tsv, scaling_tsv), "Figure 5. Executable contract validation and descriptive scaling. Synthetic checks verify the implemented status-code contract and prospective boundaries; timings describe one recorded workstation. No biological model or method comparison is represented.", "Executable contract validation")

# Figure 6: predeclared transfer benchmark.
tr <- read.delim(transfer_tsv, check.names = FALSE)
tr$region_id <- factor(tr$region_id, levels = tr$region_id)
tr$start_mb <- seq(5, 115, by = 10)
tr$class <- ifelse(tr$source_status == "E_NO_ELIGIBLE_GENE", "No eligible identifier", "Reached LD preflight")
p6_a <- ggplot(tr, aes(x = start_mb, y = 1, colour = class)) +
  geom_line(colour = INK, linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_text(aes(label = region_id), nudge_y = -0.12, size = 2.25, colour = INK) +
  scale_colour_manual(values = c(`Reached LD preflight` = BLUE, `No eligible identifier` = YELLOW), name = NULL) +
  scale_x_continuous(breaks = seq(5, 115, by = 10), labels = paste0(seq(5, 115, by = 10), " Mb")) +
  coord_cartesian(ylim = c(0.78, 1.08)) +
  labs(x = "Predeclared chr11 window start", y = NULL, title = "Deterministic region transfer") +
  theme_journal(7.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

codes <- c("E_NO_ELIGIBLE_GENE", "E_VARIANT_DUPLICATE", "E_MAF_INVALID", "E_CASE_FRACTION_HETEROGENEOUS", "E_LD_VARIANT_SET_MISMATCH", "E_LD_ILL_CONDITIONED")
heat <- expand.grid(region_id = levels(tr$region_id), code = codes) %>%
  left_join(tr %>% select(region_id, status_codes), by = "region_id") %>%
  mutate(present = mapply(function(code, status) grepl(code, status, fixed = TRUE), code, status_codes), code = factor(code, levels = rev(codes)), region_id = factor(region_id, levels = levels(tr$region_id)))
p6_b <- ggplot(heat, aes(x = region_id, y = code, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_manual(values = c(`FALSE` = PALE, `TRUE` = ORANGE), guide = "none") +
  labs(x = NULL, y = NULL, title = "Retained terminal-code matrix") +
  theme_journal(7.0) + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(size = 6.8))

cond_tr <- tr %>% filter(is.finite(ld_condition_number))
p6_c <- ggplot(cond_tr, aes(x = region_id, y = ld_condition_number)) +
  geom_hline(yintercept = 1e12, linetype = 2, linewidth = 0.45, colour = INK) +
  geom_point(size = 2.3, colour = BLUE) +
  scale_y_log10(breaks = 10^(12:15), labels = parse(text = c("10^12", "10^13", "10^14", "10^15"))) +
  labs(x = "Region", y = "Condition number", title = "LD numerical stop") +
  annotate("text", x = 1, y = 1e12, label = "supplied limit", hjust = 0, vjust = -0.5, size = 2.25) +
  theme_journal(7.2)

p6 <- (p6_a / p6_b / p6_c) + plot_layout(heights = c(0.72, 1.15, 1.0)) + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10, family = "Arial")))
save_pair(p6, "6", "figure_6", 7.2, 7.2, c(transfer_tsv, transfer_json, transfer_verification), "Figure 6. Predeclared 12-region workflow transfer benchmark. Every declared window retains a terminal status. Status-code frequencies and condition numbers are restricted to this locked benchmark; no model was run and no association, biological, causal, clinical, or method-superiority inference is supported.", "Predeclared real-data transfer benchmark")

cat("Generated journal-style figures 1, 1B, 2, 3, 4, 5, and 6 in", OUT, "\n")
