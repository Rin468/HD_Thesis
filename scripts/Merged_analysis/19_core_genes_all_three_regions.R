# ============================================================
# 19_core_genes_all_three_regions.R
#
# Purpose: confirm the genes significant in all three regions
# (striatum and cingulate cortex, Paryani et al.; motor cortex, Lin et al.),
# their direction of change, the opposite-direction exception (AHNAK2),
# and the overlap with the high-confidence lists; save one core-gene table
# (source for Section 4.9 and Appendix Table 14).
#
# Inputs:
#   main run              final_runs/bystudy_20261015_151051
#   high-confidence lists final_runs/voomexact_20261101_153230 (corrected donors)
# Outputs (new folder final_runs/coregenes_<timestamp>):
#   Genes_significant_all_three_regions.csv  (all genes significant in all three regions)
#   Core_genes_same_direction.csv            (the core genes, i.e. same direction; Table 14)
#   Core_genes_summary.txt                   (all counts and checks)
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ma   <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
runs <- file.path(ma, "final_runs")
main <- file.path(runs, "bystudy_20261015_151051")
hc   <- file.path(runs, "voomexact_20261101_153230")
out  <- file.path(runs, format(Sys.time(), "coregenes_%Y%m%d_%H%M%S"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

strip <- function(x) sub("\\..*$", "", x)   # remove the Ensembl version suffix

# ---- 1. Per-region results of the main run ------------------
rbc <- readRDS(file.path(main, "objects", "results_by_contrast.rds"))
regions <- c(A = "striatum", B_P = "cingulate", B_T = "motor_cortex")

per_region <- lapply(names(regions), function(k) {
  d <- as.data.frame(rbc[[k]])
  d <- d[!grepl("PAR_Y", rownames(d)), ]                 # drop Y-PAR copies (same ID as the X copy)
  d$gene <- strip(rownames(d))
  stopifnot(!anyDuplicated(d$gene))
  d <- d[d$reliable %in% TRUE & !is.na(d$padj),
         c("gene", "log2FC", "log2FC_shr", "padj", "sig")]
  names(d)[-1] <- paste0(names(d)[-1], "_", regions[[k]])
  d
})
m <- Reduce(function(a, b) merge(a, b, by = "gene"), per_region)   # genes tested in both studies

# ---- 2. Genes significant in all three regions --------------
s3 <- m[m$sig_striatum & m$sig_cingulate & m$sig_motor_cortex, ]

lfc_cols <- c("log2FC_striatum", "log2FC_cingulate", "log2FC_motor_cortex")
shr_cols <- c("log2FC_shr_striatum", "log2FC_shr_cingulate", "log2FC_shr_motor_cortex")
pad_cols <- c("padj_striatum", "padj_cingulate", "padj_motor_cortex")

same_mle <- apply(sign(s3[, lfc_cols]), 1, function(x) length(unique(x)) == 1)
same_shr <- apply(sign(s3[, shr_cols]), 1, function(x) length(unique(x)) == 1)
s3$direction <- ifelse(!same_mle, "opposite between regions",
                       ifelse(s3$log2FC_striatum > 0, "higher in HD", "lower in HD"))
s3$weakest_padj <- apply(s3[, pad_cols], 1, max)

# gene symbols (org.Hs.eg.db) and Paryani DESeq2 base mean
s3$symbol <- unname(mapIds(org.Hs.eg.db, s3$gene, "SYMBOL", "ENSEMBL", multiVals = "first"))
s3$symbol[is.na(s3$symbol)] <- s3$gene[is.na(s3$symbol)]
dds <- readRDS(file.path(main, "objects", "dds_ModelP_Paryani.rds"))
bm  <- mcols(dds)$baseMean
names(bm) <- strip(rownames(dds))
bm  <- bm[!duplicated(names(bm))]
s3$baseMean_Paryani <- unname(bm[s3$gene])

# high-confidence membership (corrected-donor run)
hcA  <- strip(read.csv(file.path(hc, "HighConfidence_genes_exact_true_donors_A.csv"))$gene_id)
hcBP <- strip(read.csv(file.path(hc, "HighConfidence_genes_exact_true_donors_B_P.csv"))$gene_id)
s3$highconf_striatum  <- s3$gene %in% hcA
s3$highconf_cingulate <- s3$gene %in% hcBP

s3 <- s3[order(s3$direction, -abs(rowMeans(s3[, shr_cols]))), ]
keep <- c("gene", "symbol", "direction", shr_cols, lfc_cols, pad_cols, "weakest_padj",
          "baseMean_Paryani", "highconf_striatum", "highconf_cingulate")
s3 <- s3[, keep]
core <- s3[s3$direction != "opposite between regions", ]

write.csv(s3,   file.path(out, "Genes_significant_all_three_regions.csv"), row.names = FALSE)
write.csv(core, file.path(out, "Core_genes_same_direction.csv"),          row.names = FALSE)

# ---- 3. Checks against the saved gene lists ------------------
up <- readLines(file.path(main, "gene_lists", "all_three_regions_up.txt"))
dn <- readLines(file.path(main, "gene_lists", "all_three_regions_down.txt"))

check_genes <- c("PITX1", "HOXA10", "SHOX2")
summ <- c(
  paste("Genes tested in both studies:", nrow(m)),
  paste("Significant in all three regions:", nrow(s3)),
  paste("Same direction (core genes):", nrow(core),
        "| higher:", sum(core$direction == "higher in HD"),
        "| lower:",  sum(core$direction == "lower in HD")),
  paste("Direction agrees between MLE and shrunken LFC for all genes:", all(same_mle == same_shr)),
  paste("Core 'higher' identical to all_three_regions_up.txt:",  setequal(core$gene[core$direction == "higher in HD"], up)),
  paste("Core 'lower' identical to all_three_regions_down.txt:", setequal(core$gene[core$direction == "lower in HD"], dn)),
  "Opposite-direction gene(s) (shrunken LFC striatum / cingulate / motor cortex):",
  apply(s3[s3$direction == "opposite between regions", c("symbol", shr_cols)], 1,
        function(x) paste(" ", x[1], paste(round(as.numeric(x[-1]), 2), collapse = " / "))),
  paste("Core genes in striatal high-confidence list:",  sum(core$highconf_striatum)),
  paste("Core genes in cingulate high-confidence list:", sum(core$highconf_cingulate)),
  paste("Core genes in both high-confidence lists:",
        sum(core$highconf_striatum & core$highconf_cingulate)),
  "Check genes (shrunken LFC striatum / cingulate / motor cortex; base mean):",
  apply(core[core$symbol %in% check_genes, c("symbol", shr_cols, "baseMean_Paryani")], 1,
        function(x) paste(" ", x[1], paste(round(as.numeric(x[2:4]), 2), collapse = " / "),
                          ";", round(as.numeric(x[5]), 1)))
)
writeLines(summ, file.path(out, "Core_genes_summary.txt"))
cat(summ, sep = "\n")
cat("\nSaved to:", out, "\n")

# ---- 4. Direction of the core genes in the sensitivity analyses ----
# Final sensitivity run sensitivity_20261017_163558 (per-gene results saved for these analyses only).
# S6 uses the old donor assignment; the corrected juvenile-onset run saved summary numbers only.
sens <- file.path(runs, "sensitivity_20261017_163558", "objects")
sfiles <- c("Without the HD-only batch (S1)"                      = "results_S1_without_batch4.rds",
            "Without relabelled or marker-flagged samples (S2)"   = "results_S2_without_corrected_flagged.rds",
            "Without the possible outlier control (S3)"           = "results_S3_without_C5596_Accum.rds",
            "Without juvenile-onset samples, old donors (S6)"     = "results_S6_without_juvenile_onset.rds",
            "Motor cortex without the only female control (S7)"  = "results_S7_GSE79666_without_CTRL3.rds",
            "Without the four sex-corrected samples (S9)"         = "results_S9_without_sex_corrected.rds")
main_col  <- c(A = "log2FC_striatum", B_P = "log2FC_cingulate", B_T = "log2FC_motor_cortex")
reg_label <- c(A = "striatum", B_P = "cingulate cortex", B_T = "motor cortex")

rows <- list()
for (s in names(sfiles)) {
  r <- readRDS(file.path(sens, sfiles[[s]]))
  for (k in intersect(names(main_col), names(r))) {
    d <- as.data.frame(r[[k]])
    d <- d[!grepl("PAR_Y", rownames(d)), ]
    rownames(d) <- strip(rownames(d))
    x    <- d[core$gene, ]                                   # NA rows for genes not in this analysis
    ok   <- !is.na(x$log2FC) & x$reliable %in% TRUE
    same <- sign(x$log2FC) == sign(core[[main_col[[k]]]])
    rows[[paste(s, k)]] <- data.frame(
      analysis           = s,
      region             = reg_label[[k]],
      core_genes_tested  = sum(ok),
      same_direction     = sum(ok & same),
      changed_direction  = paste(core$symbol[ok & !same], collapse = "; "),
      still_significant  = sum(ok & x$sig %in% TRUE))
  }
}
sens_tab <- do.call(rbind, rows); rownames(sens_tab) <- NULL
write.csv(sens_tab, file.path(out, "Core_genes_direction_in_sensitivity_analyses.csv"), row.names = FALSE)
print(sens_tab[, c("analysis", "region", "core_genes_tested", "same_direction", "still_significant")])
cat("Genes that changed direction:", if (all(sens_tab$changed_direction == "")) "none" else
  paste(unique(sens_tab$changed_direction[sens_tab$changed_direction != ""]), collapse = " | "), "\n")
