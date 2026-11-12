# ============================================================
# 20_Figure_workflow_overview.R
#
# Purpose: draw the overview figure of the analysis workflow
#          (10 plain boxes, one style, no software versions).
# Input:   none (the text of the boxes is written below; edit it here).
# Output:  Merged_analysis/thesis_figures_final/
#            Figure_workflow_overview_<timestamp>.png
# Needs:   ggplot2 only.
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

ma      <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
fig_dir <- file.path(ma, "thesis_figures_final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "serif"      # Times New Roman on Windows
add_title   <- TRUE         # set to FALSE to leave the title to the thesis caption

# ---- 1. Boxes (edit the text here) ----------------------------
# left / right edges: wide box 0.5-9.5; half boxes 0.5-4.8 and 5.2-9.5
box <- function(xmin, xmax, ymin, title, detail, main = FALSE)
  data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymin + 1,
             title = title, detail = detail, main = main)

boxes <- rbind(
  box(0.5, 4.8, 11.0, "Paryani et al. 2024 (GSE242195)",
      "76 post-mortem samples, HD and control\nstriatum and cingulate cortex"),
  box(5.2, 9.5, 11.0, "Lin et al. 2016 (GSE79666)",
      "14 post-mortem samples, 7 HD and 7 control\nmotor cortex"),
  box(0.5, 9.5,  9.5, "Read processing, identical for both datasets",
      "trimming (fastp), transcript quantification (Salmon), gene-level counts (tximport)"),
  box(0.5, 9.5,  8.0, "Verification of sample identity",
      "4 sex and 7 tissue labels corrected, 1 sample excluded, 3 flagged,\n2 donor assignments corrected; 75 + 14 samples analysed"),
  box(0.5, 9.5,  6.5, "Differential expression, one DESeq2 model per study",
      "HD versus control in the striatum, cingulate cortex and motor cortex"),
  box(0.5, 4.8,  5.0, "Calibration for repeated donors",
      "random-group test; donor-aware limma-voom\ngives the high-confidence genes"),
  box(5.2, 9.5,  5.0, "Robustness",
      "cell-type marker scores;\n16 sensitivity analyses"),
  box(0.5, 9.5,  3.5, "Cortex versus cortex",
      "cingulate and motor cortex agree and are combined into a pooled cortex"),
  box(0.5, 9.5,  2.0, "Striatum versus pooled cortex (main result)",
      "genes whose HD effect differs between the regions", main = TRUE),
  box(0.5, 9.5,  0.5, "Biological interpretation",
      "Gene Ontology analysis; genes changed in all three regions; comparison with the literature")
)
boxes$xmid <- (boxes$xmin + boxes$xmax) / 2

# ---- 2. Arrows (from the bottom of one row to the top of the next) ----
arrow_at <- function(x, y_from, y_to) data.frame(x = x, y = y_from, yend = y_to)
arrows <- rbind(
  arrow_at(2.65, 11.0, 10.5), arrow_at(7.35, 11.0, 10.5),   # datasets -> read processing
  arrow_at(5.00,  9.5,  9.0),                               # -> verification
  arrow_at(5.00,  8.0,  7.5),                               # -> differential expression
  arrow_at(2.65,  6.5,  6.0), arrow_at(7.35,  6.5,  6.0),   # -> calibration, robustness
  arrow_at(2.65,  5.0,  4.5), arrow_at(7.35,  5.0,  4.5),   # -> cortex versus cortex
  arrow_at(5.00,  3.5,  3.0),                               # -> striatum versus pooled cortex
  arrow_at(5.00,  2.0,  1.5)                                # -> interpretation
)

# ---- 3. Plot ---------------------------------------------------
p <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = main),
            colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "grey88"), guide = "none") +
  geom_text(data = boxes, aes(x = xmid, y = ymax - 0.27, label = title),
            family = font_family, fontface = "bold", size = 3.7) +
  geom_text(data = boxes, aes(x = xmid, y = ymin + 0.36, label = detail),
            family = font_family, size = 3.2, lineheight = 0.95) +
  geom_segment(data = arrows, aes(x = x, xend = x, y = y, yend = yend),
               arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
               linewidth = 0.4) +
  coord_cartesian(xlim = c(0.3, 9.7), ylim = c(0.3, 12.2), expand = FALSE) +
  theme_void(base_family = font_family) +
  theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5,
                                     margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 8)),
        plot.margin   = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", colour = NA))

if (add_title) {
  p <- p + labs(title    = "Overview of the analysis",
                subtitle = "Each box is one step; arrows show the order. The grey box is the main result.")
}

# ---- 4. Save (time-stamped, never overwrites) ------------------
png_file <- file.path(fig_dir, format(Sys.time(), "Figure_workflow_overview_%Y%m%d_%H%M%S.png"))
ggsave(png_file, p, width = 7, height = 9.5, dpi = 300)
cat("Saved:", png_file, "\n")
