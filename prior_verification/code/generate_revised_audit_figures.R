suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(digest)
  library(ggrepel)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript generate_revised_audit_figures.R <summary.tsv> <out_dir>")
source_tsv <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(source_tsv, check.names = FALSE, stringsAsFactors = FALSE)
d$case_label <- ifelse(d$case_id == "coordinate_case", "C2", d$case_id)

ink <- "#1F2933"; blue <- "#2E74B5"; teal <- "#168C84"; orange <- "#D97706"; grey <- "#718096"; pale <- "#EEF3F7"
theme_journal <- theme_classic(base_size = 10, base_family = "sans") +
  theme(axis.title = element_text(colour = ink), axis.text = element_text(colour = ink),
        plot.title = element_text(face = "bold", colour = ink, size = 12),
        plot.subtitle = element_text(colour = grey, size = 9), legend.position = "bottom",
        plot.margin = margin(8, 12, 8, 8))

save_figure <- function(plot, number, stem, width, height, caption, contract) {
  png <- file.path(out_dir, paste0(stem, ".png"))
  svg <- file.path(out_dir, paste0(stem, ".svg"))
  ggsave(png, plot, width = width, height = height, dpi = 600, bg = "white", device = grDevices::png, type = "cairo-png")
  ggsave(svg, plot, width = width, height = height, device = svglite, bg = "white")
  manifest <- list(
    schema_version = "2.0", figure_number = number, caption = caption,
    figure_contract = contract,
    sources = list(list(path = source_tsv, sha256 = digest(file = source_tsv, algo = "sha256"))),
    outputs = list(
      png = list(path = png, sha256 = digest(file = png, algo = "sha256"), dpi = 600),
      svg = list(path = svg, sha256 = digest(file = svg, algo = "sha256"))
    ), R = R.version.string, ggplot2 = as.character(packageVersion("ggplot2"))
  )
  write_json(manifest, file.path(out_dir, paste0(stem, "_manifest.json")), auto_unbox = TRUE, pretty = TRUE)
}

flow <- data.frame(x = 1:5, y = 1, label = c("Immutable\ninputs", "Full row\naudit", "Deterministic\nanalysis view", "Declared\ngate policy", "Separate\nmodel plan"))
p1 <- ggplot(flow, aes(x, y)) +
  geom_segment(aes(x = 1.3, xend = 4.7, yend = 1), colour = grey, linewidth = 0.7, arrow = arrow(length = unit(0.12, "inches"))) +
  geom_label(aes(label = label), fill = pale, colour = ink, linewidth = 0.3, size = 3.1, lineheight = 0.95) +
  annotate("text", x = 4, y = 0.72, label = "READY = declared-policy conformance only", colour = teal, size = 3) +
  annotate("text", x = 5, y = 1.27, label = "Model eligibility is not inferred", colour = orange, size = 3) +
  coord_cartesian(xlim = c(0.6, 5.4), ylim = c(0.55, 1.45), clip = "off") +
  labs(title = "Audit-first interface", subtitle = "Integrity checks and model-owned numerical policies are explicit") + theme_void(base_size = 10)
save_figure(p1, "1", "figure_1_revised_contract", 7.2, 2.5,
  "Figure 1. Corrected audit-first interface. READY denotes conformance to the named gate policy only; the gate does not infer model eligibility or a biological result.",
  "workflow-only; no association, posterior, biological, clinical, or performance claim")

rows_long <- rbind(
  data.frame(case_id = d$case_label, category = "Retained analysis rows", value = d$analysis_rows),
  data.frame(case_id = d$case_label, category = "Exact duplicates excluded", value = d$exact_duplicates_excluded),
  data.frame(case_id = d$case_label, category = "Not in locked LD excluded", value = d$rows_not_in_locked_ld_excluded)
)
rows_long <- subset(rows_long, case_id != "R06")
p2 <- ggplot(rows_long, aes(x = case_id, y = value, fill = category)) +
  geom_col(width = 0.72) + scale_fill_manual(values = c(teal, blue, orange)) +
  labs(title = "Audited construction of the analysis view", subtitle = "All exclusions remain in row-level audit files", x = NULL, y = "Rows", fill = NULL) +
  theme_journal + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p2, "2", "figure_2_revised_analysis_view", 7.2, 4.2,
  "Figure 2. Deterministic construction of the unique analysis view. Exact duplicate and LD-absent rows are excluded with explicit row-level decisions; no discordant duplicate was found in these stored artifacts.",
  "workflow-only; no association, posterior, biological, clinical, or performance claim")

rank_d <- subset(d, !is.na(numerical_rank) & ld_variants > 0)
p3 <- ggplot(rank_d, aes(x = ld_variants, y = numerical_rank, label = case_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = grey) +
  geom_hline(yintercept = 502, linetype = 3, colour = orange) +
  geom_point(size = 2.8, colour = blue) + geom_text_repel(size = 2.7, colour = ink, min.segment.length = 0, seed = 1, box.padding = 0.25) +
  coord_cartesian(ylim = c(0, max(rank_d$ld_variants) * 1.03)) +
  labs(title = "Stored LD matrices are rank deficient", subtitle = "The 11 p>503 windows have a structural ceiling of 502; C2 deficiency is empirical", x = "LD variants", y = "Numerical rank") + theme_journal
save_figure(p3, "3", "figure_3_revised_rank", 7.2, 4.5,
  "Figure 3. Numerical rank versus LD dimension in the stored additional-coordinate and transfer artifacts. All matrices are rank deficient. The 11 transfer matrices exceed the 503-person reference dimension and have a centered rank ceiling of 502; C2 contains 225 variants, so its empirical deficiency requires a separate construction-level explanation.",
  "workflow-only numerical diagnostic; no association, posterior, biological, clinical, or performance claim")

status_d <- d
status_d$status_label <- ifelse(status_d$gate_status == "READY", "Policy conformance\nREADY/OK", "No eligible gene\nINELIGIBLE")
p4 <- ggplot(status_d, aes(x = case_label, y = 1, fill = status_label)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(rank_deficient == "True", paste0(numerical_rank, "\n/ ", ld_variants), "NO\nGENE")), size = 2.4, colour = "white", lineheight = 0.9) +
  scale_fill_manual(values = c("No eligible gene\nINELIGIBLE" = orange, "Policy conformance\nREADY/OK" = teal)) +
  labs(title = "Post-hoc policy re-audit", subtitle = "Study-level case fraction; rank recorded diagnostically; model eligibility remains unset", x = NULL, y = NULL, fill = NULL) +
  theme_journal + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p4, "4", "figure_4_revised_policy_reaudit", 7.2, 3.6,
  "Figure 4. Post-hoc re-audit under the named study-level case-fraction and diagnostic-rank policies. READY is policy conformance only; the figure is not an eligibility, association, or validation result.",
  "workflow-only; no association, posterior, biological, clinical, or performance claim")

evidence <- data.frame(
  tier = factor(
    c("Unit/contract tests", "Synthetic mutations", "Cross-decomposition check",
      "Stored real artifacts", "Synthetic downstream stress test",
      "External biological validation"),
    levels = rev(c("Unit/contract tests", "Synthetic mutations", "Cross-decomposition check",
                   "Stored real artifacts", "Synthetic downstream stress test",
                   "External biological validation"))
  ),
  completed = c(1, 1, 1, 1, 1, 0),
  label = c("93/93 tests", "2,500/2,500 named codes", "80/80 rank comparisons",
            "Post-hoc traceability only", "100 simulations x 5 conditions", "Not performed")
)
p5 <- ggplot(evidence, aes(x = completed, y = tier, fill = factor(completed))) +
  geom_col(width = 0.62) + geom_text(aes(label = label), hjust = ifelse(evidence$completed == 1, 1.03, -0.03), colour = ifelse(evidence$completed == 1, "white", ink), size = 3.0) +
  scale_fill_manual(values = c("0" = "#D7DEE5", "1" = blue), guide = "none") + coord_cartesian(xlim = c(0, 1.12), clip = "off") +
  labs(title = "Evidence hierarchy and remaining boundary", subtitle = "Synthetic consumer stress tests do not establish biological validity", x = NULL, y = NULL) +
  theme_journal + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(p5, "5", "figure_5_revised_evidence_boundary", 7.2, 5.0,
  "Figure 5. Evidence hierarchy. Contract tests, fixed mutations, cross-decomposition checks, stored-artifact traceability, and a synthetic downstream-consumer stress test were completed; external biological validation was not performed.",
  "software and synthetic downstream-consumer evidence only; no association, biological, clinical, or comparative-performance claim")
