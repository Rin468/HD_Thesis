# =============================================================================
#  08_regional_effect_crossstudy.R   (written 31 Oct 2026)
#  Question: are the genes whose HD effect differs between striatum and cingulate
#  cortex in the Paryani data (interaction contrast, 589 genes) also different
#  when the cortex side is replaced by the independent motor-cortex data of
#  Lin et al. (GSE79666)?
#
#  Uses the saved main-run models only (no refit unless the object is unfitted).
#  Same DEG rule as the main run: padj < 0.05, excluding genes with zero counts in
#  a group of the contrast or |log2FC| > 10. Unshrunken (MLE) log2FC, as in the
#  between-dataset z-test of Methods 3.6.
#
#  Self-checks against the main run: striatum 6,267, cingulate 1,000,
#  interaction 589, motor cortex 714, striatum-vs-motor-cortex differences 1,635.
#  If any check does not match, the script still runs but prints a WARNING:
#  send the log to Claude before using the numbers.
#
#  NOTE: the striatal estimate is the same in both comparisons, so the two are
#  not independent; the check shows whether the region difference survives when
#  the cingulate data are swapped for another cortex dataset.
#
#  Output: NEW folder Merged_analysis/final_runs/regional_crossstudy_<date>_<time>/
#  Run:    Session > Restart R, open this file, Source.
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
out_dir <- file.path(base, "final_runs", format(Sys.time(), "regional_crossstudy_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
check <- function(label, got, expected) {
  say(sprintf("%-55s %6d  (main run: %d)%s", label, got, expected,
              if (got != expected) "   <-- WARNING: does not match" else ""))
}

# ---- 1. Load the two main-run models -------------------------------------------
obj_dir <- file.path(run, "objects")
fP <- list.files(obj_dir, pattern = "^dds_ModelP.*\\.rds$", full.names = TRUE)
if (length(fP) != 1) { print(list.files(obj_dir)); stop("Paryani model file not found uniquely - send Claude the list above.") }
ddsP <- readRDS(fP)
ddsL <- readRDS(file.path(obj_dir, "dds_ModelL_GSE79666.rds"))
if (is.null(mcols(ddsP)$dispersion)) { say("Paryani model not fitted - fitting now (several minutes)"); ddsP <- DESeq(ddsP) }
if (is.null(mcols(ddsL)$dispersion)) { say("Lin model not fitted - fitting now"); ddsL <- DESeq(ddsL) }
say("Paryani model: ", basename(fP), " | design: ", deparse(design(ddsP)), " | samples: ", ncol(ddsP))
say("Lin model: design ", deparse(design(ddsL)), " | samples: ", ncol(ddsL))
say("Paryani group levels: ", paste(levels(ddsP$group), collapse = ", "))

# ---- 2. Paryani contrasts as numeric vectors on the group coefficients ----------
rn  <- resultsNames(ddsP)
ref <- levels(ddsP$group)[1]
g <- function(lv) {
  v <- setNames(numeric(length(rn)), rn)
  if (lv == ref) return(v)
  nm <- paste0("group_", lv, "_vs_", ref)
  if (!nm %in% rn) stop("Coefficient ", nm, " not found. resultsNames: ", paste(rn, collapse = ", "))
  v[nm] <- 1; v
}
cnt <- counts(ddsP)
nonzero <- function(lv) rowSums(cnt[, ddsP$group == lv, drop = FALSE]) > 0
tab_P <- function(vec, groups) {
  r <- results(ddsP, contrast = vec, alpha = 0.05)
  rel <- Reduce(`&`, lapply(groups, nonzero)) & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10)
  data.frame(gene_id = rownames(r), b = r$log2FoldChange, se = r$lfcSE, padj = r$padj,
             tested = rel & !is.na(r$padj), sig = rel & !is.na(r$padj) & r$padj < 0.05)
}
A  <- tab_P(g("Striatum_HD") - g("Striatum_CTRL"), c("Striatum_HD", "Striatum_CTRL"))
BP <- tab_P(g("Cortex_HD")   - g("Cortex_CTRL"),   c("Cortex_HD", "Cortex_CTRL"))
EP <- tab_P(g("Striatum_HD") - g("Striatum_CTRL") - g("Cortex_HD") + g("Cortex_CTRL"),
            c("Striatum_HD", "Striatum_CTRL", "Cortex_HD", "Cortex_CTRL"))

# ---- 3. Lin motor cortex -------------------------------------------------------
rL <- results(ddsL, name = "condition_HD_vs_CTRL", alpha = 0.05)
cL <- counts(ddsL)
relL <- rowSums(cL[, ddsL$condition == "HD"]) > 0 & rowSums(cL[, ddsL$condition == "CTRL"]) > 0 &
  !(!is.na(rL$log2FoldChange) & abs(rL$log2FoldChange) > 10)
BT <- data.frame(gene_id = rownames(rL), b = rL$log2FoldChange, se = rL$lfcSE, padj = rL$padj,
                 tested = relL & !is.na(rL$padj), sig = relL & !is.na(rL$padj) & rL$padj < 0.05)

say("\n-- Self-checks against the main run --")
check("Striatum HD vs control DEGs", sum(A$sig), 6267)
check("Cingulate HD vs control DEGs", sum(BP$sig), 1000)
check("Interaction (striatal minus cingulate HD effect) DEGs", sum(EP$sig), 589)
check("Motor cortex HD vs control DEGs", sum(BT$sig), 714)

# ---- 4. Striatum vs motor cortex: difference of independent estimates ------------
ids <- intersect(A$gene_id[A$tested], BT$gene_id[BT$tested])
a <- A[match(ids, A$gene_id), ]; t <- BT[match(ids, BT$gene_id), ]
e <- EP[match(ids, EP$gene_id), ]; c <- BP[match(ids, BP$gene_id), ]
d_ST <- a$b - t$b
z    <- d_ST / sqrt(a$se^2 + t$se^2)
p_z  <- 2 * pnorm(-abs(z)); padj_z <- p.adjust(p_z, "BH")
sig_ST <- padj_z < 0.05
say("\nGenes tested in striatum and motor cortex: ", length(ids), " (main run text: 23,730 tested in both studies)")
check("Striatum-vs-motor-cortex HD-effect differences", sum(sig_ST), 1635)

# ---- 5. Do the region-dependent Paryani genes also differ against motor cortex? ---
inE    <- e$sig %in% TRUE
same_d <- sign(e$b) == sign(d_ST)
cl_cing <- abs(t$b - c$b) < abs(t$b - a$b)          # motor-cortex effect closer to cingulate than to striatum
ft <- fisher.test(table(interaction_gene = inE, differs_vs_motor_cortex = sig_ST))
summ <- data.frame(
  measure = c("Interaction genes also tested in motor cortex",
              "  of these: striatal HD effect also differs from motor cortex (padj < 0.05)",
              "  of these: same direction of the difference",
              "  of all interaction genes: same direction of the difference (any significance)",
              "  of all interaction genes: motor-cortex effect closer to cingulate than to striatal effect",
              "All other tested genes: striatal effect differs from motor cortex (padj < 0.05)",
              "Odds ratio (interaction genes vs other genes, differing vs motor cortex)",
              "Fisher test p-value",
              "Spearman rho, interaction estimate vs striatum-minus-motor-cortex, all genes",
              "Spearman rho, same, interaction genes only"),
  value = c(sum(inE),
            sum(inE & sig_ST), sum(inE & sig_ST & same_d, na.rm = TRUE),
            sum(inE & same_d, na.rm = TRUE), sum(inE & cl_cing, na.rm = TRUE),
            sum(!inE & sig_ST),
            unname(round(ft$estimate, 2)), signif(ft$p.value, 3),
            round(cor(e$b, d_ST, method = "spearman", use = "complete.obs"), 3),
            round(cor(e$b[inE], d_ST[inE], method = "spearman", use = "complete.obs"), 3)),
  percent = c(NA,
              round(100 * sum(inE & sig_ST) / sum(inE), 1),
              round(100 * sum(inE & sig_ST & same_d, na.rm = TRUE) / sum(inE), 1),
              round(100 * sum(inE & same_d, na.rm = TRUE) / sum(inE), 1),
              round(100 * sum(inE & cl_cing, na.rm = TRUE) / sum(inE), 1),
              round(100 * sum(!inE & sig_ST) / sum(!inE), 1),
              NA, NA, NA, NA))
say("\n-- Region-dependent genes against the independent motor-cortex data --")
cat(capture.output(print(summ, right = FALSE)), sep = "\n")
cat(capture.output(print(summ, right = FALSE)), sep = "\n", file = log_file, append = TRUE)
write.csv(summ, file.path(out_dir, "Regional_effect_crossstudy_summary.csv"), row.names = FALSE)

genes <- data.frame(gene_id = ids,
                    striatum_HD_effect = round(a$b, 3), cingulate_HD_effect = round(c$b, 3),
                    motor_cortex_HD_effect = round(t$b, 3),
                    interaction_striatum_minus_cingulate = round(e$b, 3), interaction_padj = signif(e$padj, 3),
                    striatum_minus_motor_cortex = round(d_ST, 3), striatum_vs_motor_cortex_padj = signif(padj_z, 3),
                    interaction_gene = inE, differs_vs_motor_cortex = sig_ST)
write.csv(genes[genes$interaction_gene, ], file.path(out_dir, "Regional_effect_crossstudy_interaction_genes.csv"), row.names = FALSE)

# ---- 6. Figure (plain labels) ---------------------------------------------------
genes$class <- factor(ifelse(genes$interaction_gene, "Region-dependent HD effect in Paryani data", "Other genes"),
                      levels = c("Other genes", "Region-dependent HD effect in Paryani data"))
genes <- genes[order(genes$class), ]
p <- ggplot(genes, aes(interaction_striatum_minus_cingulate, striatum_minus_motor_cortex, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c("Other genes" = "grey80", "Region-dependent HD effect in Paryani data" = "#D55E00")) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(x = "Striatal minus cingulate HD effect (log2 fold change, Paryani et al.)",
       y = "Striatal HD effect minus motor-cortex HD effect (log2 fold change, Lin et al.)",
       colour = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "Regional_effect_crossstudy_scatter.png"), p, width = 6.5, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "Regional_effect_crossstudy_scatter.pdf"), p, width = 6.5, height = 6.5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude: run_log.txt (or the console output) and Regional_effect_crossstudy_summary.csv")
say("Outputs: ", out_dir)
