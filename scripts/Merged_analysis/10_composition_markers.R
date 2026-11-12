# =============================================================================
#  10_composition_markers.R   (written 1 Nov 2026)
#  Question: how much of the HD signal reflects altered cell composition
#  (neuron loss, astrocytosis) rather than changed expression within cells?
#  What:
#    1. Marker-gene scores per sample for neurons, medium spiny neurons (MSNs),
#       astrocytes, microglia and oligodendrocytes (mean z-score of VST
#       expression, z-scored within tissue).
#    2. HD vs control test of each score per tissue. Paryani samples are first
#       averaged per donor within tissue (caudate + accumbens of one donor count
#       once), so repeated donors are not treated as independent.
#    3. Sensitivity model: the main DESeq2 design plus the neuronal score as a
#       covariate. Genes that stay significant change beyond what neuron loss
#       explains. CAUTION: neuron loss is part of HD itself, so this model also
#       removes some genuine disease signal (over-adjustment); report it as a
#       sensitivity analysis, not as the corrected result.
#  Run time: about 10-30 min (two DESeq2 refits).
#  Output: NEW folder Merged_analysis/final_runs/composition_<date>_<time>/
#  Run:    Session > Restart R, open this file, Source.
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(AnnotationDbi); library(org.Hs.eg.db) })

base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
out_dir <- file.path(base, "final_runs", format(Sys.time(), "composition_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }

ddsP <- readRDS(file.path(run, "objects", "dds_ModelP_Paryani.rds"))
ddsL <- readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds"))
if (!"donor" %in% names(colData(ddsP))) stop("No 'donor' column in the Paryani model - send Claude names(colData(ddsP)).")

# ---- 1. Marker sets (canonical brain cell-type markers) ---------------------------
markers <- list(
  "Neurons"                = c("RBFOX3", "SNAP25", "SYT1", "STMN2"),
  "Medium spiny neurons"   = c("PPP1R1B", "PDE10A", "DRD1", "DRD2", "PENK", "ADORA2A", "GPR88"),
  "Astrocytes"             = c("GFAP", "AQP4", "ALDH1L1", "SLC1A2"),
  "Microglia"              = c("AIF1", "CX3CR1", "P2RY12", "CSF1R"),
  "Oligodendrocytes"       = c("MBP", "PLP1", "MOG", "MOBP"))

score_tab <- function(dds, tissue) {
  vs  <- assay(vst(dds, blind = TRUE))
  ens <- sub("\\..*$", "", rownames(vs))
  sy  <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  out <- data.frame(sample = colnames(dds), tissue = tissue, condition = ifelse(grepl("HD", tissue_cond(dds)), "HD", "Control"))
  for (s in names(markers)) {
    rows <- which(sy %in% markers[[s]])
    found <- unique(sy[rows])
    say(sprintf("%-22s markers found: %s", s, paste(found, collapse = ", ")))
    z <- vs[rows, , drop = FALSE]
    # z-score each gene within tissue, then average genes
    for (tt in unique(tissue)) { i <- tissue == tt; z[, i] <- t(scale(t(z[, i, drop = FALSE]))) }
    out[[s]] <- colMeans(z, na.rm = TRUE)
  }
  out
}
tissue_cond <- function(dds) if ("group" %in% names(colData(dds))) as.character(dds$group) else as.character(dds$condition)

say("== Paryani ==")
tP <- ifelse(grepl("^Striatum", ddsP$group), "Striatum", "Cingulate cortex")
scP <- score_tab(ddsP, tP); scP$donor <- as.character(ddsP$donor); scP$study <- "Paryani et al."
say("== Lin et al. (motor cortex) ==")
scL <- score_tab(ddsL, rep("Motor cortex", ncol(ddsL))); scL$donor <- scL$sample; scL$study <- "Lin et al."
sc <- rbind(scP, scL)
write.csv(sc, file.path(out_dir, "Composition_scores_per_sample.csv"), row.names = FALSE)

# ---- 2. HD vs control per tissue, one value per donor and tissue -------------------
tests <- list()
for (tt in unique(sc$tissue)) for (s in names(markers)) {
  d <- sc[sc$tissue == tt, ]
  dd <- aggregate(d[[s]], by = list(donor = d$donor, condition = d$condition), FUN = mean)
  w <- suppressWarnings(wilcox.test(x ~ condition, data = dd, exact = FALSE))
  tests[[length(tests) + 1]] <- data.frame(tissue = tt, cell_type = s,
                                           donors_HD = sum(dd$condition == "HD"), donors_control = sum(dd$condition == "Control"),
                                           median_HD = round(median(dd$x[dd$condition == "HD"]), 2), median_control = round(median(dd$x[dd$condition == "Control"]), 2),
                                           wilcoxon_p = signif(w$p.value, 3))
}
tests <- do.call(rbind, tests); tests$padj_BH <- signif(p.adjust(tests$wilcoxon_p, "BH"), 3)
show("Cell-type marker scores, HD vs control (one value per donor and tissue)", tests)
write.csv(tests, file.path(out_dir, "Composition_scores_HD_vs_control.csv"), row.names = FALSE)

sc$condition <- factor(sc$condition, levels = c("Control", "HD"))
long <- do.call(rbind, lapply(names(markers), function(s) data.frame(tissue = sc$tissue, condition = sc$condition, cell_type = s, score = sc[[s]])))
long$cell_type <- factor(long$cell_type, levels = names(markers))
long$tissue <- factor(long$tissue, levels = c("Striatum", "Cingulate cortex", "Motor cortex"))
p <- ggplot(long, aes(condition, score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  facet_grid(tissue ~ cell_type) +
  scale_fill_manual(values = c(Control = "#0072B2", HD = "#D55E00"), guide = "none") +
  labs(x = NULL, y = "Marker-gene score (mean z-score within tissue)") + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "Composition_scores_by_tissue.png"), p, width = 10, height = 6.5, dpi = 300)

# ---- 3. Main design + neuronal score as covariate -----------------------------------
deg <- function(r, cnt, grp_ok) {
  rel <- grp_ok & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10)
  data.frame(b = r$log2FoldChange, padj = r$padj, sig = rel & !is.na(r$padj) & r$padj < 0.05, rel = rel)
}
cmp <- function(main, alt, label) {
  same <- sign(main$b) == sign(alt$b); ok <- main$rel & alt$rel & is.finite(main$b) & is.finite(alt$b)
  data.frame(comparison = label, DEGs_main = sum(main$sig), DEGs_with_neuronal_score = sum(alt$sig),
             main_DEGs_retained_pct = round(100 * sum(main$sig & alt$sig & same, na.rm = TRUE) / sum(main$sig), 1),
             spearman_rho_log2FC = round(cor(main$b[ok], alt$b[ok], method = "spearman"), 3))
}
# Paryani
ddsN <- ddsP; ddsN$neuronal_score <- scP[["Neurons"]]
design(ddsN) <- ~ batch_study + sex + neuronal_score + group
say("\nFitting Paryani model with neuronal score (several minutes)...")
ddsN <- DESeq(ddsN, quiet = TRUE)
rn <- resultsNames(ddsP); ref <- levels(ddsP$group)[1]
gv <- function(rnx, lv) { v <- setNames(numeric(length(rnx)), rnx); if (lv == ref) return(v); v[paste0("group_", lv, "_vs_", ref)] <- 1; v }
cnt <- counts(ddsP); nz <- function(lv) rowSums(cnt[, ddsP$group == lv, drop = FALSE]) > 0
cons <- list(A = c("Striatum_HD", "Striatum_CTRL"), B_P = c("Cortex_HD", "Cortex_CTRL"),
             E_P = c("Striatum_HD", "Striatum_CTRL", "Cortex_HD", "Cortex_CTRL"))
labs <- c(A = "Striatum, HD vs control", B_P = "Cingulate cortex, HD vs control", E_P = "Region-dependent HD effect (interaction)")
vec <- function(rnx, k) if (k == "E_P") gv(rnx, "Striatum_HD") - gv(rnx, "Striatum_CTRL") - gv(rnx, "Cortex_HD") + gv(rnx, "Cortex_CTRL") else
  gv(rnx, cons[[k]][1]) - gv(rnx, cons[[k]][2])
res <- list()
for (k in names(cons)) {
  ok <- Reduce(`&`, lapply(cons[[k]], nz))
  m  <- deg(results(ddsP, contrast = vec(rn, k), alpha = 0.05), cnt, ok)
  a  <- deg(results(ddsN, contrast = vec(resultsNames(ddsN), k), alpha = 0.05), cnt, ok)
  res[[k]] <- cmp(m, a, labs[[k]])
}
# Lin et al.
ddsLN <- ddsL; ddsLN$neuronal_score <- scL[["Neurons"]]
design(ddsLN) <- ~ sex + neuronal_score + condition
ddsLN <- DESeq(ddsLN, quiet = TRUE)
cL <- counts(ddsL); okL <- rowSums(cL[, ddsL$condition == "HD"]) > 0 & rowSums(cL[, ddsL$condition == "CTRL"]) > 0
res$BT <- cmp(deg(results(ddsL,  name = "condition_HD_vs_CTRL", alpha = 0.05), cL, okL),
              deg(results(ddsLN, name = "condition_HD_vs_CTRL", alpha = 0.05), cL, okL), "Motor cortex, HD vs control")
res <- do.call(rbind, res)
show("Main design vs main design + neuronal marker score (self-check: main DEGs 6267, 1000, 589, 714)", res)
write.csv(res, file.path(out_dir, "Composition_covariate_vs_main.csv"), row.names = FALSE)
saveRDS(ddsN, file.path(out_dir, "dds_ModelP_Paryani_neuronal_score.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude the console output (or run_log.txt) and Composition_scores_by_tissue.png. Outputs: ", out_dir)