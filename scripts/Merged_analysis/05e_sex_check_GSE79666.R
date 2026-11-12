# =============================================================================
#  05e_sex_check_GSE79666.R
#  Sensitivity check: GSE79666 (Lin et al. 2016, motor cortex BA4) HD vs CTRL
#  Main design ~ sex + condition  vs  design without sex ~ condition.
#  Same input object, same genes and same DEG rule as the main run and S10 (05d).
#  Output folder: results/sensitivity_sex_GSE79666
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
out_dir <- file.path(base, "results/sensitivity_sex_GSE79666")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load the main-run GSE79666 object ---------------------------------------
dds_L <- readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds"))
print(table(sex = dds_L$sex, condition = dds_L$condition))

# ---- Same DEG rule as the main run and S10 (copied from 05d) ------------------
# padj < 0.05, excluding genes with zero counts in a group or |log2FC| > 10
deg_tab <- function(dds) {
  res <- results(dds, name = "condition_HD_vs_CTRL", alpha = 0.05)
  cnt <- counts(dds)
  reliable <- rowSums(cnt[, dds$condition == "HD"]) > 0 &
              rowSums(cnt[, dds$condition == "CTRL"]) > 0 &
              !(!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 10)
  data.frame(gene_id = rownames(res), log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
             pvalue = res$pvalue, padj = res$padj, reliable = reliable,
             sig = reliable & !is.na(res$padj) & res$padj < 0.05)
}

# ---- Refit both designs on the same genes ------------------------------------
dds_ref   <- dds_L; design(dds_ref)   <- ~ sex + condition; dds_ref   <- DESeq(dds_ref,   quiet = TRUE)
dds_nosex <- dds_L; design(dds_nosex) <- ~ condition;       dds_nosex <- DESeq(dds_nosex, quiet = TRUE)

ref <- deg_tab(dds_ref)
nos <- deg_tab(dds_nosex)
stopifnot(identical(ref$gene_id, nos$gene_id))
cat("Main-model DEGs (should be 714, as in the main run):", sum(ref$sig), "\n")

# ---- Gene symbol and chromosome (X / Y flag) ---------------------------------
ens <- sub("\\..*$", "", ref$gene_id)                       # strip Ensembl version
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
cyt <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "MAP",    multiVals = "first")
chr <- sub("[pq].*$", "", cyt)                              # "Xq13.2" -> "X"

# ---- Comparison metrics (same as S10) ----------------------------------------
same <- sign(ref$log2FC) == sign(nos$log2FC)
ok   <- ref$reliable & nos$reliable & is.finite(ref$log2FC) & is.finite(nos$log2FC)
only_nosex <- nos$sig & !ref$sig          # DEGs that appear only without sex
lost       <- ref$sig & !nos$sig          # main DEGs lost without sex

comp <- data.frame(
  analysis               = c("Main model: sex + condition", "Without sex: condition only"),
  DEGs                   = c(sum(ref$sig), sum(nos$sig)),
  up_in_HD               = c(sum(ref$sig & ref$log2FC > 0), sum(nos$sig & nos$log2FC > 0)),
  down_in_HD             = c(sum(ref$sig & ref$log2FC < 0), sum(nos$sig & nos$log2FC < 0)),
  main_DEGs_retained_pct = c(NA, round(100 * sum(ref$sig & nos$sig & same, na.rm = TRUE) / sum(ref$sig), 1)),
  new_DEGs_only_without_sex = c(NA, sum(only_nosex)),
  of_these_on_X_or_Y     = c(NA, sum(only_nosex & chr %in% c("X", "Y"), na.rm = TRUE)),
  main_DEGs_lost         = c(NA, sum(lost)),
  direction_flips        = c(NA, sum(ref$sig & nos$sig & !same, na.rm = TRUE)),
  spearman_rho_log2FC    = c(NA, round(cor(ref$log2FC[ok], nos$log2FC[ok], method = "spearman"), 3)))
print(comp)
write.csv(comp, file.path(out_dir, "SexCheck_GSE79666_without_sex_vs_main.csv"), row.names = FALSE)

# ---- Gene lists: new without sex, and lost without sex -----------------------
gene_list <- function(idx, label) {
  if (!any(idx)) return(NULL)
  d <- data.frame(change = label, gene_id = ref$gene_id[idx], gene_symbol = sym[idx],
                  chromosome = chr[idx],
                  log2FC_main = round(ref$log2FC[idx], 3), padj_main = signif(ref$padj[idx], 3),
                  log2FC_without_sex = round(nos$log2FC[idx], 3), padj_without_sex = signif(nos$padj[idx], 3),
                  row.names = NULL)
  d[order(d$padj_without_sex), ]
}
genes <- rbind(gene_list(only_nosex, "Significant only without sex"),
               gene_list(lost,       "Significant only with sex"))
print(head(genes, 30))
write.csv(genes, file.path(out_dir, "SexCheck_GSE79666_changed_genes.csv"), row.names = FALSE)

# ---- Same comparison table as S10, for the robustness table ------------------
cat("\nChromosome of the genes that appear only without sex:\n")
print(table(chromosome = chr[only_nosex], useNA = "ifany"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
