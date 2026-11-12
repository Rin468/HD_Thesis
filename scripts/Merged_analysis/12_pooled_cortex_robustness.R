# =============================================================================
#  12_pooled_cortex_robustness.R   (written 1 Nov 2026)
#  Question: which of the 1,867 genes whose HD effect differs between the striatum
#  and the pooled cortex (script 11) are robust to
#   (a) the repeated-donor problem  -> repeat the pooled analysis with the
#       donor-aware limma-voom estimates for the Paryani contrasts
#       (calibrated in the null calibration); "high-confidence" = significant in
#       both versions, same direction (same idea as script 09);
#   (b) the HD-only sequencing batch 4 -> refit the Paryani DESeq2 model without
#       batch 4 and repeat the pooled analysis.
#  The motor-cortex estimate (Lin et al.) is the same DESeq2 estimate in all versions
#  (one sample per donor, no batch problem).
#
#  Input : latest final_runs/pooled_cortex_*/PooledCortex_all_genes.rds (script 11),
#          SENS_voom_donor_results.csv (sensitivity run), the two main DESeq2 models.
#  Output: NEW folder Merged_analysis/final_runs/pooled_robustness_<date>_<time>/
#  Run:    close other programs, Session > Restart R, open this file, Source.
#          About 5-10 minutes (the model without batch 4 is refitted).
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2) })

base <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run  <- file.path(base, "final_runs/bystudy_20261015_151051")
sens <- file.path(base, "final_runs/sensitivity_20261017_163558")
pc_dirs <- sort(list.dirs(file.path(base, "final_runs"), recursive = FALSE, full.names = TRUE))
pc_dirs <- pc_dirs[grepl("^pooled_cortex_[0-9]", basename(pc_dirs)) &
                     file.exists(file.path(pc_dirs, "PooledCortex_all_genes.rds"))]
stopifnot(length(pc_dirs) > 0)
pc_dir  <- tail(pc_dirs, 1)
out_dir <- file.path(base, "final_runs", format(Sys.time(), "pooled_robustness_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }
slim <- function(d) { for (a in c("mu", "H", "replaceCounts", "replaceCooks")) assays(d)[[a]] <- NULL; d }

# ---- 1. Main pooled result (script 11) ---------------------------------------------
G <- readRDS(file.path(pc_dir, "PooledCortex_all_genes.rds"))
ids <- G$gene_id
main_sig <- G$region_difference_padj < 0.05
d_main   <- G$striatum_minus_pooled_cortex
say("Pooled run used: ", basename(pc_dir))
say("Self-check (script 11: 22371 genes, 1867 region-dependent): ", nrow(G), ", ", sum(main_sig))

# ---- 2. Motor cortex (Lin et al.), same estimate in every version ------------------------
ddsL <- slim(readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds")))
rL <- results(ddsL, name = "condition_HD_vs_CTRL", alpha = 0.05)
t_b <- rL$log2FoldChange[match(ids, rownames(rL))]; t_se <- rL$lfcSE[match(ids, rownames(rL))]
rm(ddsL, rL); invisible(gc())

# striatum minus pooled cortex, same formulas as script 11
pool_diff <- function(sb, sse, cb, cse, ese) {
  a    <- (1 / cse^2) / (1 / cse^2 + 1 / t_se^2)
  b_cx <- a * cb + (1 - a) * t_b
  se_cx <- sqrt(a^2 * cse^2 + (1 - a)^2 * t_se^2)
  cov_SC <- (sse^2 + cse^2 - ese^2) / 2
  d <- sb - b_cx; v <- sse^2 + se_cx^2 - 2 * a * cov_SC
  ok <- is.finite(d) & is.finite(v) & v > 0
  p <- rep(NA_real_, length(d)); p[ok] <- 2 * pnorm(-abs(d[ok] / sqrt(v[ok])))
  list(d = d, padj = p.adjust(p, "BH"), n_tested = sum(ok), n_nonpos_var = sum(is.finite(v) & v <= 0))
}
robust <- function(x) main_sig & (x$padj < 0.05) %in% TRUE & (sign(x$d) == sign(d_main)) %in% TRUE

# ---- 3. (a) Donor-aware limma-voom version -------------------------------------------
vd <- read.csv(file.path(sens, "tables", "SENS_voom_donor_results.csv"))
v  <- vd[match(ids, vd$gene_id), ]
say("\nGenes found in the donor-aware limma-voom results: ", sum(!is.na(v$gene_id)), " of ", length(ids))
say("Check voom interaction = striatum - cingulate, max |difference|: ",
    signif(max(abs(v$log2FC_E_P - (v$log2FC_A - v$log2FC_B_P)), na.rm = TRUE), 3))
z_int <- 2 * pnorm(-abs(v$log2FC_E_P / v$lfcSE_E_P))
say("Check SE columns: Spearman rho of z-based and reported voom interaction p-values: ",
    round(cor(p.adjust(z_int, "BH"), v$padj_E_P, method = "spearman", use = "complete.obs"), 3),
    " (should be close to 1)")
PV <- pool_diff(v$log2FC_A, v$lfcSE_A, v$log2FC_B_P, v$lfcSE_B_P, v$lfcSE_E_P)
hc <- robust(PV)
say("voom version: genes tested ", PV$n_tested, ", non-positive variances ", PV$n_nonpos_var,
    ", region-dependent genes ", sum(PV$padj < 0.05, na.rm = TRUE))

# ---- 4. (b) Paryani model refitted without the HD-only batch ---------------------------
ddsP <- readRDS(file.path(run, "objects", "dds_ModelP_Paryani.rds"))
hd_only <- names(which(tapply(grepl("HD", ddsP$group), ddsP$batch_study, all)))
say("\nHD-only batch removed: ", paste(hd_only, collapse = ", "), " (", sum(ddsP$batch_study %in% hd_only), " samples)")
d1 <- slim(ddsP[, !(ddsP$batch_study %in% hd_only)]); rm(ddsP); invisible(gc())
colData(d1) <- droplevels(colData(d1))
SH <- "Striatum_HD"; SC <- "Striatum_CTRL"; CH <- "Cortex_HD"; CC <- "Cortex_CTRL"
stopifnot(all(c(SH, SC, CH, CC) %in% levels(d1$group)))
say("Samples per group without batch 4: ", paste(names(table(d1$group)), table(d1$group), collapse = ", "))
say("Refitting DESeq2 (a few minutes)...")
d1 <- DESeq(d1, quiet = TRUE); invisible(gc())
rn1 <- resultsNames(d1); ref1 <- levels(d1$group)[1]
g1 <- function(lv) { x <- setNames(numeric(length(rn1)), rn1); if (lv == ref1) return(x)
x[paste0("group_", lv, "_vs_", ref1)] <- 1; x }
cnt1 <- counts(d1); grp1 <- d1$group
tab1 <- function(vec, groups) {   # same estimability rule as the main run
  r <- results(d1, contrast = vec, alpha = 0.05)
  nz  <- function(lv) rowSums(cnt1[, grp1 == lv, drop = FALSE]) > 0
  rel <- Reduce(`&`, lapply(groups, nz)) & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10) & !is.na(r$padj)
  i <- match(ids, rownames(r))
  data.frame(b = ifelse(rel[i] %in% TRUE, r$log2FoldChange[i], NA), se = ifelse(rel[i] %in% TRUE, r$lfcSE[i], NA))
}
S1 <- tab1(g1(SH) - g1(SC), c(SH, SC)); C1 <- tab1(g1(CH) - g1(CC), c(CH, CC))
E1 <- tab1(g1(SH) - g1(SC) - g1(CH) + g1(CC), c(SH, SC, CH, CC))
rm(d1, cnt1); invisible(gc())
P1 <- pool_diff(S1$b, S1$se, C1$b, C1$se, E1$se)
wo4 <- robust(P1)
say("Without batch 4: genes tested ", P1$n_tested, ", non-positive variances ", P1$n_nonpos_var,
    ", region-dependent genes ", sum(P1$padj < 0.05, na.rm = TRUE))

# ---- 5. Summary ---------------------------------------------------------------------
inE <- G$interaction_gene_Paryani %in% TRUE
neg <- d_main < 0
row <- function(lab, x) data.frame(gene_set = lab, genes = sum(x),
                                   more_negative_in_striatum = sum(x & neg), more_positive_in_striatum = sum(x & !neg),
                                   pct_of_1867 = round(100 * sum(x) / sum(main_sig), 1), also_interaction_gene = sum(x & inE))
summ <- rbind(
  row("All region-dependent genes (striatum vs pooled cortex)", main_sig),
  row("High-confidence (also in donor-aware limma-voom version)", hc),
  row("Significant without batch 4", wo4),
  row("High-confidence and significant without batch 4", hc & wo4))
show("Robustness of the striatum vs pooled cortex genes", summ)
write.csv(summ, file.path(out_dir, "PooledRobustness_summary.csv"), row.names = FALSE)

rho <- data.frame(comparison = c("Main vs donor-aware limma-voom version", "Main vs without batch 4"),
                  spearman_rho_all_genes = c(round(cor(d_main, PV$d, method = "spearman", use = "complete.obs"), 3),
                                             round(cor(d_main, P1$d, method = "spearman", use = "complete.obs"), 3)),
                  spearman_rho_1867_genes = c(round(cor(d_main[main_sig], PV$d[main_sig], method = "spearman", use = "complete.obs"), 3),
                                              round(cor(d_main[main_sig], P1$d[main_sig], method = "spearman", use = "complete.obs"), 3)),
                  same_direction_pct_1867 = c(round(100 * mean(sign(PV$d[main_sig]) == sign(d_main[main_sig]), na.rm = TRUE), 1),
                                              round(100 * mean(sign(P1$d[main_sig]) == sign(d_main[main_sig]), na.rm = TRUE), 1)))
show("Agreement of the region-dependent effect (striatum minus pooled cortex)", rho)
write.csv(rho, file.path(out_dir, "PooledRobustness_agreement.csv"), row.names = FALSE)

flags <- data.frame(gene_id = ids, gene_symbol = G$gene_symbol,
                    striatum_HD_effect = G$striatum_HD_effect, pooled_cortex_HD_effect = G$pooled_cortex_HD_effect,
                    striatum_minus_pooled_cortex = d_main, region_difference_padj = G$region_difference_padj,
                    voom_striatum_minus_pooled_cortex = round(PV$d, 3), voom_padj = signif(PV$padj, 3),
                    without_batch4_striatum_minus_pooled_cortex = round(P1$d, 3), without_batch4_padj = signif(P1$padj, 3),
                    high_confidence = hc, significant_without_batch4 = wo4, interaction_gene_Paryani = inE)[main_sig, ]
flags <- flags[order(!(flags$high_confidence & flags$significant_without_batch4), flags$region_difference_padj), ]
write.csv(flags, file.path(out_dir, "PooledRobustness_region_dependent_genes_flags.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude the console output (or run_log.txt). Outputs: ", out_dir)