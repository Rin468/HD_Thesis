# =============================================================================
#  11_pooled_cortex.R   (written 1 Nov 2026)
#  Question: which genes respond to HD differently in the striatum than in the
#  cortex, when the cortex is the POOLED cortex of both studies
#  (cingulate cortex, Paryani et al. + motor cortex, Lin et al.)?
#
#  Method (meta-analysis of the per-study estimates, no new model fit):
#   1. Pooled cortical HD effect = inverse-variance weighted mean of the cingulate
#      and motor-cortex HD effects (fixed-effect meta-analysis). Each study keeps
#      its own dispersion, so the problem of one joint DESeq2 model is avoided.
#   2. Region-dependent effect = striatal HD effect minus pooled cortical HD effect.
#      Its variance accounts for the covariance between the striatal and cingulate
#      estimates (same model, overlapping donors), derived exactly from the SEs of
#      the striatum, cingulate and interaction contrasts:
#        Cov(S, C) = (SE_S^2 + SE_C^2 - SE_interaction^2) / 2
#      Var(S - pooled) = SE_S^2 + SE_pooled^2 - 2 * a * Cov(S, C),  a = cingulate weight
#   3. Heterogeneity: do the two cortices disagree? (z-test of cingulate vs motor cortex)
#   4. Sensitivity: equal weights for the two cortices instead of inverse-variance.
#  Unshrunken log2 fold changes and Wald SEs, Benjamini-Hochberg over tested genes.
#
#  Output: NEW folder Merged_analysis/final_runs/pooled_cortex_<date>_<time>/
#  Run:    Session > Restart R, open this file, Source.  (about 1-2 minutes)
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
sens    <- file.path(base, "final_runs/sensitivity_20261017_163558")
out_dir <- file.path(base, "final_runs", format(Sys.time(), "pooled_cortex_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }

# ---- 1. Per-study estimates (same rule as the main run) ---------------------------
# Low-memory version: one model in memory at a time; assays that results() does not use are dropped.
slim <- function(d) { for (a in c("mu", "H", "replaceCounts", "replaceCooks")) assays(d)[[a]] <- NULL; d }
ddsP <- slim(readRDS(file.path(run, "objects", "dds_ModelP_Paryani.rds"))); invisible(gc())
rn <- resultsNames(ddsP); ref <- levels(ddsP$group)[1]
g <- function(lv) { v <- setNames(numeric(length(rn)), rn); if (lv == ref) return(v)
v[paste0("group_", lv, "_vs_", ref)] <- 1; v }
cnt <- counts(ddsP); grp <- ddsP$group
tab_P <- function(vec, groups) {
  r <- results(ddsP, contrast = vec, alpha = 0.05)
  nz <- function(lv) rowSums(cnt[, grp == lv, drop = FALSE]) > 0
  rel <- Reduce(`&`, lapply(groups, nz)) & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10)
  out <- data.frame(gene_id = rownames(r), b = r$log2FoldChange, se = r$lfcSE, padj = r$padj,
                    tested = rel & !is.na(r$padj), sig = rel & !is.na(r$padj) & r$padj < 0.05)
  rm(r); invisible(gc()); out
}
SH <- "Striatum_HD"; SC <- "Striatum_CTRL"; CH <- "Cortex_HD"; CC <- "Cortex_CTRL"
S <- tab_P(g(SH) - g(SC), c(SH, SC))
C <- tab_P(g(CH) - g(CC), c(CH, CC))
E <- tab_P(g(SH) - g(SC) - g(CH) + g(CC), c(SH, SC, CH, CC))
rm(ddsP, cnt); invisible(gc())
ddsL <- slim(readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds")))
rL <- results(ddsL, name = "condition_HD_vs_CTRL", alpha = 0.05); cL <- counts(ddsL)
relL <- rowSums(cL[, ddsL$condition == "HD"]) > 0 & rowSums(cL[, ddsL$condition == "CTRL"]) > 0 &
  !(!is.na(rL$log2FoldChange) & abs(rL$log2FoldChange) > 10)
Tm <- data.frame(gene_id = rownames(rL), b = rL$log2FoldChange, se = rL$lfcSE, padj = rL$padj,
                 tested = relL & !is.na(rL$padj), sig = relL & !is.na(rL$padj) & rL$padj < 0.05)
say("Self-check (main run: striatum 6267, cingulate 1000, interaction 589, motor cortex 714): ",
    sum(S$sig), ", ", sum(C$sig), ", ", sum(E$sig), ", ", sum(Tm$sig))

# genes estimable in all four contrasts
ids <- Reduce(intersect, list(S$gene_id[S$tested], C$gene_id[C$tested], E$gene_id[E$tested], Tm$gene_id[Tm$tested]))
s <- S[match(ids, S$gene_id), ]; c <- C[match(ids, C$gene_id), ]
e <- E[match(ids, E$gene_id), ]; t <- Tm[match(ids, Tm$gene_id), ]
say("Genes tested in striatum, cingulate, interaction and motor cortex: ", length(ids))
say("Check interaction = striatum - cingulate, max |difference|: ", signif(max(abs(e$b - (s$b - c$b))), 3))

bh <- function(z) p.adjust(2 * pnorm(-abs(z)), "BH")

# ---- 2. Pooled cortical HD effect ------------------------------------------------------
pool <- function(a) {               # a = weight of the cingulate estimate (vector or scalar)
  b_cx  <- a * c$b + (1 - a) * t$b
  se_cx <- sqrt(a^2 * c$se^2 + (1 - a)^2 * t$se^2)
  cov_SC <- (s$se^2 + c$se^2 - e$se^2) / 2
  d      <- s$b - b_cx
  var_d  <- s$se^2 + se_cx^2 - 2 * a * cov_SC
  list(b_cx = b_cx, se_cx = se_cx, padj_cx = bh(b_cx / se_cx),
       d = d, se_d = sqrt(pmax(var_d, 1e-12)), padj_d = bh(d / sqrt(pmax(var_d, 1e-12))), neg_var = sum(var_d <= 0))
}
w_c <- 1 / c$se^2; w_t <- 1 / t$se^2
a_iv <- w_c / (w_c + w_t)
P  <- pool(a_iv)      # main: inverse-variance weights
P5 <- pool(0.5)       # sensitivity: equal weights
say("Median weight of the cingulate estimate in the pooled cortex: ", round(median(a_iv), 3),
    " (motor cortex: ", round(1 - median(a_iv), 3), ")")
say("Genes with non-positive variance (should be 0): ", P$neg_var, " / ", P5$neg_var)

# heterogeneity between the two cortices
padj_het <- bh((c$b - t$b) / sqrt(c$se^2 + t$se^2))

summ <- function(lab, x) {
  sigC <- x$padj_cx < 0.05; sigD <- x$padj_d < 0.05; inE <- e$sig
  data.frame(weights = lab,
             pooled_cortex_DEGs = sum(sigC), pooled_cortex_higher = sum(sigC & x$b_cx > 0), pooled_cortex_lower = sum(sigC & x$b_cx < 0),
             striatum_vs_pooled_cortex_genes = sum(sigD), more_negative_in_striatum = sum(sigD & x$d < 0), more_positive_in_striatum = sum(sigD & x$d > 0),
             of_these_also_interaction_gene = sum(sigD & inE),
             interaction_genes_recovered_pct = round(100 * sum(sigD & inE & sign(x$d) == sign(e$b)) / sum(inE), 1),
             spearman_rho_with_interaction = round(cor(x$d, e$b, method = "spearman"), 3))
}
res <- rbind(summ("inverse-variance (main)", P), summ("equal weights (sensitivity)", P5))
show("Pooled cortex (cingulate + motor cortex) and striatum vs pooled cortex", res)
say("Genes whose HD effect differs between the two cortices (heterogeneity, padj < 0.05): ", sum(padj_het < 0.05))
write.csv(res, file.path(out_dir, "PooledCortex_summary.csv"), row.names = FALSE)

# categories as in the thesis (striatum vs pooled cortex)
sigS <- s$sig; sigX <- P$padj_cx < 0.05; sigD <- P$padj_d < 0.05
cat_tab <- data.frame(category = c("Significant in both, higher in HD", "Significant in both, lower in HD",
                                   "Significant in both, opposite directions", "Striatum only, higher in HD",
                                   "Striatum only, lower in HD", "Pooled cortex only, higher in HD",
                                   "Pooled cortex only, lower in HD", "Not significant in either"),
                      genes = c(sum(sigS & sigX & s$b > 0 & P$b_cx > 0), sum(sigS & sigX & s$b < 0 & P$b_cx < 0),
                                sum(sigS & sigX & sign(s$b) != sign(P$b_cx)), sum(sigS & !sigX & s$b > 0), sum(sigS & !sigX & s$b < 0),
                                sum(!sigS & sigX & P$b_cx > 0), sum(!sigS & sigX & P$b_cx < 0), sum(!sigS & !sigX)))
cat_tab$with_significant_region_difference <- c(
  sum(sigS & sigX & s$b > 0 & P$b_cx > 0 & sigD), sum(sigS & sigX & s$b < 0 & P$b_cx < 0 & sigD),
  sum(sigS & sigX & sign(s$b) != sign(P$b_cx) & sigD), sum(sigS & !sigX & s$b > 0 & sigD), sum(sigS & !sigX & s$b < 0 & sigD),
  sum(!sigS & sigX & P$b_cx > 0 & sigD), sum(!sigS & sigX & P$b_cx < 0 & sigD), sum(!sigS & !sigX & sigD))
show("HD effect categories, striatum vs pooled cortex (main weights)", cat_tab)
write.csv(cat_tab, file.path(out_dir, "PooledCortex_categories.csv"), row.names = FALSE)

# high-confidence flag for the striatum (donor-aware limma-voom), if available
vd <- tryCatch(read.csv(file.path(sens, "tables", "SENS_voom_donor_results.csv")), error = function(e) NULL)
sym <- if (!is.null(vd)) setNames(vd$gene_symbol, vd$gene_id)[ids] else rep(NA, length(ids))
genes <- data.frame(gene_id = ids, gene_symbol = unname(sym),
                    striatum_HD_effect = round(s$b, 3), cingulate_HD_effect = round(c$b, 3), motor_cortex_HD_effect = round(t$b, 3),
                    cingulate_weight = round(a_iv, 3), pooled_cortex_HD_effect = round(P$b_cx, 3), pooled_cortex_padj = signif(P$padj_cx, 3),
                    striatum_minus_pooled_cortex = round(P$d, 3), region_difference_padj = signif(P$padj_d, 3),
                    region_difference_padj_equal_weights = signif(P5$padj_d, 3), cortex_heterogeneity_padj = signif(padj_het, 3),
                    interaction_gene_Paryani = e$sig)
write.csv(genes[genes$region_difference_padj < 0.05, ], file.path(out_dir, "PooledCortex_region_dependent_genes.csv"), row.names = FALSE)
write.csv(genes[genes$pooled_cortex_padj < 0.05, ], file.path(out_dir, "PooledCortex_HD_vs_control_DEGs.csv"), row.names = FALSE)
saveRDS(genes, file.path(out_dir, "PooledCortex_all_genes.rds"))

# ---- 3. Figure (plain labels) -------------------------------------------------------
genes$class <- factor(ifelse(genes$region_difference_padj < 0.05,
                             ifelse(genes$striatum_minus_pooled_cortex < 0, "HD effect more negative in striatum", "HD effect more positive in striatum"),
                             "No significant regional difference"),
                      levels = c("No significant regional difference", "HD effect more negative in striatum", "HD effect more positive in striatum"))
genes <- genes[order(genes$class), ]
p <- ggplot(genes, aes(striatum_HD_effect, pooled_cortex_HD_effect, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_colour_manual(values = c("No significant regional difference" = "grey80",
                                 "HD effect more negative in striatum" = "#0072B2",
                                 "HD effect more positive in striatum" = "#D55E00")) +
  coord_cartesian(xlim = c(-5, 6), ylim = c(-3, 6)) +
  labs(x = "Striatal HD effect (log2 fold change, Paryani et al.)",
       y = "Pooled cortical HD effect (log2 fold change;\ncingulate cortex and motor cortex)", colour = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "PooledCortex_striatum_vs_pooled_cortex.png"), p, width = 6.5, height = 6.5, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude the console output (or run_log.txt). Outputs: ", out_dir)