# 11c_plot_striatum_vs_cortex_summary.R - summary figure of region-dependent HD effects
# Numbers are the confirmed values of Table 8:
#   589 (447 / 142): main run bystudy_20261015_151051, interaction contrast
#   1,867 (1,316 / 551), 1,477 (1,058 / 419) and the eight groups: pooled_cortex_20261101_123913 (scripts 11, 11b)
library(ggplot2); library(patchwork)

neg <- "#2F6DAE"; pos <- "#C8622A"
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

left <- data.frame(
  comparison = rep(c("Striatum vs cingulate cortex\n(within Paryani et al., one model)",
                     "Striatum vs pooled cortex\n(main comparison)",
                     "Striatum vs pooled cortex\n(equal weights, sensitivity)"), each = 2),
  direction  = rep(c("HD effect more negative in the striatum",
                     "HD effect more positive in the striatum"), 3),
  genes      = c(447, 142, 1316, 551, 1058, 419))
left$comparison <- factor(left$comparison, levels = rev(unique(left$comparison)))
left$direction  <- factor(left$direction, levels = rev(unique(left$direction)))   # negative drawn first
tot <- aggregate(genes ~ comparison, left, sum)

p1 <- ggplot(left, aes(genes, comparison, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = fmt(genes)), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.6) +
  geom_text(data = tot, aes(genes, comparison, label = fmt(genes)), inherit.aes = FALSE,
            hjust = -0.15, size = 4) +
  scale_fill_manual(values = c("HD effect more negative in the striatum" = neg,
                               "HD effect more positive in the striatum" = pos),
                    breaks = c("HD effect more negative in the striatum",
                               "HD effect more positive in the striatum"), name = NULL) +
  scale_x_continuous(limits = c(0, 2300), expand = c(0, 0)) +
  labs(title = "How many genes differ between the regions", y = NULL,
       x = "Genes whose HD effect differs significantly\nbetween striatum and cortex (adjusted p < 0.05)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")

right <- data.frame(
  group = c("Striatum only, lower in HD", "Striatum only, higher in HD",
            "Both regions, lower in HD", "Both regions, higher in HD",
            "Both regions, opposite directions",
            "Pooled cortex only, lower in HD", "Pooled cortex only, higher in HD",
            "Neither region on its own"),
  genes = c(1169, 438, 71, 18, 34, 41, 26, 70),
  col   = c(neg, pos, neg, pos, "#444444", neg, pos, "grey70"))
stopifnot(sum(right$genes) == 1867, sum(left$genes[3:4]) == 1867)
right$group <- factor(right$group, levels = rev(right$group))

p2 <- ggplot(right, aes(genes, group, fill = col)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = fmt(genes)), hjust = -0.2, size = 4) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1350), breaks = seq(0, 1200, 200), expand = c(0, 0)) +
  labs(title = "Where these 1,867 genes change in HD", y = NULL,
       x = "Region-dependent genes of the main comparison\n(striatum vs pooled cortex, 1,867 genes)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_line(colour = "grey90"))

fig <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
out <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis/thesis_figures"
ggsave(file.path(out, "Figure_striatum_vs_cortex_summary.png"), fig, width = 13, height = 5, dpi = 300)
ggsave(file.path(out, "Figure_striatum_vs_cortex_summary.pdf"), fig, width = 13, height = 5)
write.csv(left,  file.path(out, "Figure_striatum_vs_cortex_summary_left_panel.csv"),  row.names = FALSE)
write.csv(right[, c("group", "genes")], file.path(out, "Figure_striatum_vs_cortex_summary_right_panel.csv"), row.names = FALSE)
cat("Saved to", out, "\n")
