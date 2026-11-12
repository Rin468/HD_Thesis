# ============================================================
# 18b_Figure10_GO_robust_region_dependent_genes.R
#
# Purpose: draw Figure 10 (bar chart) from the GO results of script 18.
#   Gene Ontology biological-process terms over-represented among the
#   704 robust region-dependent genes whose HD effect is more negative in
#   the striatum than in the pooled cortex (top 15 non-redundant terms).
#
# Input:  final_runs/GO_regional_20261101_163253/GO_BP_rob_neg_nonredundant.csv (script 18)
# Output: Merged_analysis/thesis_figures_final/
#           Figure10_GO_robust_region_dependent_genes_<timestamp>.png
#           Figure10_GO_robust_region_dependent_genes_<timestamp>_plotted_terms.csv
#
# Rules (no hand-picking):
#   terms         = the 15 non-redundant terms with the smallest adjusted p-value
#   example genes = the first five genes of each term as listed by clusterProfiler
#   colour groups = assigned below by term name (the script stops if a term is unassigned)
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

ma      <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
gdir    <- file.path(ma, "final_runs", "GO_regional_20261101_163253")
fig_dir <- file.path(ma, "thesis_figures_final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

n_terms    <- 15
n_examples <- 5
add_title  <- TRUE   # set to FALSE to leave the title to the thesis caption

# ---- 1. Top terms --------------------------------------------
go <- read.csv(file.path(gdir, "GO_BP_rob_neg_nonredundant.csv"))
go <- go[order(go$p.adjust), ][seq_len(n_terms), ]

# ---- 2. Process groups (colours) -----------------------------
grp_levels <- c("Protein synthesis (ribosomal proteins)",
                "Dopamine/cAMP signalling of striatal neurons and behaviour",
                "Synaptic transmission, ion and calcium transport",
                "Energy metabolism")
groups <- c(
  "cytoplasmic translation"                                                 = grp_levels[1],
  "learning or memory"                                                      = grp_levels[2],
  "adenylate cyclase-modulating G protein-coupled receptor signaling pathway" = grp_levels[2],
  "regulation of G protein-coupled receptor signaling pathway"              = grp_levels[2],
  "locomotory behavior"                                                     = grp_levels[2],
  "response to amphetamine"                                                 = grp_levels[2],
  "modulation of chemical synaptic transmission"                            = grp_levels[3],
  "regulation of monoatomic ion transport"                                  = grp_levels[3],
  "regulation of monoatomic cation transmembrane transport"                 = grp_levels[3],
  "regulation of membrane potential"                                        = grp_levels[3],
  "calcium ion transport"                                                   = grp_levels[3],
  "calcium ion transmembrane import into cytosol"                           = grp_levels[3],
  "signal release"                                                          = grp_levels[3],
  "calcium ion transmembrane transport"                                     = grp_levels[3],
  "ATP biosynthetic process"                                                = grp_levels[4])
missing <- setdiff(go$Description, names(groups))
if (length(missing)) stop("No colour group for: ", paste(missing, collapse = "; "))
go$group <- factor(groups[go$Description], levels = grp_levels)
grp_cols <- setNames(c("#2F6DB0", "#C8622A", "#3E9B74", "#8C8C8C"), grp_levels)

# ---- 3. Labels -------------------------------------------------
lab <- paste0(toupper(substring(go$Description, 1, 1)), substring(go$Description, 2))
lab <- gsub("signaling pathway", "signalling", lab)
lab <- gsub("signaling", "signalling", lab)
lab <- gsub("behavior", "behaviour", lab)
go$label    <- vapply(lab, function(x) paste(strwrap(x, 34), collapse = "\n"), "")
go$neglog   <- -log10(go$p.adjust)
go$examples <- vapply(strsplit(go$geneID, "/"),
                      function(g) paste(head(g, n_examples), collapse = ", "), "")
go$y        <- rev(seq_len(nrow(go)))                       # most significant at the top

inside       <- go$neglog > 0.6 * max(go$neglog)              # long bars: text inside the bar
go$text_x    <- ifelse(inside, go$neglog - 0.1, go$neglog + 0.1)
go$hjust     <- ifelse(inside, 1, 0)
go$col_count <- ifelse(inside, "white", "black")
go$col_genes <- ifelse(inside, "white", "grey35")

# ---- 4. Plot ---------------------------------------------------
p <- ggplot(go) +
  geom_rect(aes(xmin = 0, xmax = neglog, ymin = y - 0.31, ymax = y + 0.31, fill = group)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  annotate("text", x = -log10(0.05) + 0.1, y = max(go$y) + 0.65,
           label = "adjusted p = 0.05", hjust = 0, size = 3, colour = "grey30") +
  geom_text(aes(x = text_x, y = y + 0.17, label = paste(Count, "genes"),
                hjust = hjust, colour = col_count), size = 3.4) +
  geom_text(aes(x = text_x, y = y - 0.18, label = examples,
                hjust = hjust, colour = col_genes), size = 2.7, fontface = "italic") +
  scale_colour_identity() +
  scale_fill_manual(values = grp_cols, name = NULL, drop = FALSE) +
  scale_y_continuous(breaks = go$y, labels = go$label, expand = expansion(add = c(0.6, 0.9))) +
  scale_x_continuous(limits = c(0, max(go$neglog) * 1.07), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = expression("Significance of enrichment, " * -log[10] * "(adjusted p-value)"), y = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(colour = "black"),
        legend.position    = "bottom",
        plot.title         = element_text(face = "bold", size = 12),
        plot.subtitle      = element_text(size = 10, colour = "grey25"))

if (add_title) {
  p <- p + labs(
    title    = "Biological processes of genes with a more negative HD effect in the striatum than in the cortex",
    subtitle = paste0("Top ", n_terms, " non-redundant Gene Ontology biological-process terms among the ",
                      "704 robust region-dependent genes\nBars: significance; labels: number of genes ",
                      "and the first five genes of each term"))
}

# ---- 5. Save (time-stamped, never overwrites) ------------------
stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
png_file <- file.path(fig_dir, paste0("Figure10_GO_robust_region_dependent_genes_", stamp, ".png"))
ggsave(png_file, p, width = 10, height = if (add_title) 8.2 else 7.5, dpi = 300)
write.csv(go[, c("ID", "Description", "group", "Count", "p.adjust", "neglog", "examples")],
          sub("\\.png$", "_plotted_terms.csv", png_file), row.names = FALSE)

print(go[, c("Description", "group", "Count", "p.adjust")], row.names = FALSE)
cat("\nSaved:", png_file, "\n")

