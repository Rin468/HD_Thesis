# =============================================================================
#  05d_S10_age_GSE79666.R
#  RECOVERED on 2026-10-31 from RStudio history_database (18 Oct 2026 console session).
#  S10 sensitivity: GSE79666 (Lin et al. 2016, motor cortex BA4) HD vs CTRL,
#  main design ~ sex + condition vs age-adjusted ~ sex + age_c + condition.
#  Outputs: results/tables/GSE79666_GEO_characteristics.csv
#           results/sensitivity_age_GSE79666/  (age table, Wilcoxon test, S10 comparison,
#           all-gene table, diagnostics, dispersion ratios, dds_ModelL_GSE79666_age_adjusted.rds)
#  Code is verbatim; one-off inspection lines were left out. Indentation lost in history.
# =============================================================================

# ---- Part 1: GEO sample characteristics of GSE79666 ---------------------------
if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
library(GEOquery)
ma  <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
gse <- getGEO("GSE79666", GSEMatrix = TRUE)
ph  <- pData(gse[[1]])
print(names(ph)[grepl(":ch1$", names(ph))])     # every characteristic GEO provides
get <- function(col) if (col %in% names(ph)) ph[[col]] else NA
lin_meta <- data.frame(
  sample_id    = ph$geo_accession,
  sample_label = ph$title,
  tissue       = get("tissue:ch1"),
  age          = suppressWarnings(as.numeric(get("age:ch1"))),
  diagnosis    = get("diagnosis:ch1"),
  stringsAsFactors = FALSE)
lin_meta$condition <- ifelse(grepl("^HD", lin_meta$sample_label), "HD", "CTRL")
print(lin_meta)
out <- file.path(ma, "results", "tables", "GSE79666_GEO_characteristics.csv")
write.csv(lin_meta, out, row.names = FALSE)
# Is age balanced between the groups?
print(tapply(lin_meta$age, lin_meta$condition, summary))
print(wilcox.test(age ~ condition, data = lin_meta, exact = FALSE))

# ---- Part 2: age-adjusted refit (S10) -------------------------------------------
suppressPackageStartupMessages(library(DESeq2))
base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
out_dir <- file.path(base, "results/sensitivity_age_GSE79666")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dds_L <- readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds"))
# Age and grade from GEO/BioSample (check grades against your GEOquery table)
ages  <- c("HD-57"=65,"HD-58"=44,"HD-61"=69,"HD-68"=56,"HD-283"=51,"HD-288"=57,"HD-290"=52,
           "CTRL-3"=46,"CTRL-4"=42,"CTRL-8"=50,"CTRL-9"=48,"CTRL-11"=57,"CTRL-12"=64,"CTRL-13"=59)
grade <- c("HD-57"=2,"HD-58"=3,"HD-61"=3,"HD-68"=2,"HD-283"=4,"HD-288"=3,"HD-290"=3)
stopifnot(setequal(names(ages), dds_L$sample_label))
dds_L$age   <- unname(ages[dds_L$sample_label])
dds_L$age_c <- dds_L$age - mean(dds_L$age)        # centred; does not change the HD vs CTRL estimate
# 1) Save the per-sample table and the age-balance test
age_tab <- data.frame(sample = dds_L$sample_label, condition = dds_L$condition,
                      sex_inferred = dds_L$sex, age_years = dds_L$age,
                      Vonsattel_grade = unname(grade[dds_L$sample_label]))
age_tab <- age_tab[order(age_tab$condition, age_tab$sample), ]
write.csv(age_tab, file.path(out_dir, "GSE79666_age_grade_by_sample.csv"), row.names = FALSE)
wt <- wilcox.test(age_years ~ condition, data = age_tab, exact = FALSE)
write.csv(data.frame(test = "Wilcoxon rank-sum, age by condition", W = unname(wt$statistic),
                     p_value = wt$p.value), file.path(out_dir, "GSE79666_age_balance_test.csv"), row.names = FALSE)
# 2) Refit the main design (reference) and the age-adjusted design on the same genes
dds_ref <- dds_L; design(dds_ref) <- ~ sex + condition;         dds_ref <- DESeq(dds_ref, quiet = TRUE)
dds_age <- dds_L; design(dds_age) <- ~ sex + age_c + condition; dds_age <- DESeq(dds_age, quiet = TRUE)
cat("Non-converged genes, age model:", sum(!(mcols(dds_age)$betaConv %in% TRUE)), "\n")
# Same DEG rule as the main run: padj < 0.05, excluding genes with zero counts in a group or |log2FC| > 10
deg_tab <- function(dds) {
  res <- results(dds, name = "condition_HD_vs_CTRL", alpha = 0.05)
  cnt <- counts(dds)
  reliable <- rowSums(cnt[, dds$condition == "HD"]) > 0 & rowSums(cnt[, dds$condition == "CTRL"]) > 0 &
    !(!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 10)
  data.frame(gene_id = rownames(res), log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
             pvalue = res$pvalue, padj = res$padj, reliable = reliable,
             sig = reliable & !is.na(res$padj) & res$padj < 0.05)
}
ref <- deg_tab(dds_ref); agd <- deg_tab(dds_age)
stopifnot(identical(ref$gene_id, agd$gene_id))
same <- sign(ref$log2FC) == sign(agd$log2FC)
ok   <- ref$reliable & agd$reliable & is.finite(ref$log2FC) & is.finite(agd$log2FC)
age_res <- results(dds_age, name = "age_c", alpha = 0.05)
comp <- data.frame(
  analysis          = c("Main model: sex + condition", "Age-adjusted: sex + age + condition"),
  DEGs              = c(sum(ref$sig), sum(agd$sig)),
  up_in_HD          = c(sum(ref$sig & ref$log2FC > 0), sum(agd$sig & agd$log2FC > 0)),
  down_in_HD        = c(sum(ref$sig & ref$log2FC < 0), sum(agd$sig & agd$log2FC < 0)),
  main_DEGs_retained_pct = c(NA, round(100 * sum(ref$sig & agd$sig & same, na.rm = TRUE) / sum(ref$sig), 1)),
  new_DEGs_only_with_age = c(NA, sum(agd$sig & !ref$sig)),
  spearman_rho_log2FC    = c(NA, round(cor(ref$log2FC[ok], agd$log2FC[ok], method = "spearman"), 3)),
  genes_associated_with_age_padj05 = c(NA, sum(!is.na(age_res$padj) & age_res$padj < 0.05)))
print(comp)
write.csv(comp, file.path(out_dir, "S10_GSE79666_age_adjusted_vs_main.csv"), row.names = FALSE)
write.csv(agd,  file.path(out_dir, "S10_GSE79666_age_adjusted_HD_vs_CTRL_all_genes.csv"), row.names = FALSE)
saveRDS(dds_age, file.path(out_dir, "dds_ModelL_GSE79666_age_adjusted.rds"))
lost <- ref$sig & !agd$sig
diag_tab <- data.frame(
  check = c("Main DEGs lost with age adjustment",
            "  of these: padj with age 0.05-0.10",
            "  of these: padj with age > 0.20",
            "Direction flips among main DEGs",
            "Median SE ratio (age / main), main DEGs",
            "Median |log2FC| ratio (age / main), main DEGs",
            "Correlation of age with HD status",
            "Genes tested, main model", "Genes tested, age model"),
  value = c(sum(lost),
            sum(lost & agd$padj >= 0.05 & agd$padj < 0.10, na.rm = TRUE),
            sum(lost & agd$padj > 0.20, na.rm = TRUE),
            sum(ref$sig & !same, na.rm = TRUE),
            round(median(agd$lfcSE[ref$sig] / ref$lfcSE[ref$sig]), 3),
            round(median(abs(agd$log2FC[ref$sig]) / abs(ref$log2FC[ref$sig])), 3),
            round(cor(dds_L$age, as.numeric(dds_L$condition == "HD")), 3),
            sum(!is.na(ref$padj)), sum(!is.na(agd$padj))))
print(diag_tab)
print(table(main_padj = cut(ref$padj[ref$sig], c(0, 0.001, 0.01, 0.05)),
            still_DEG_with_age = agd$sig[ref$sig]))
write.csv(diag_tab, file.path(out_dir, "S10_GSE79666_age_diagnostics.csv"), row.names = FALSE)
cat("Median dispersion ratio (age / main), main DEGs:",
    round(median(dispersions(dds_age)[ref$sig] / dispersions(dds_ref)[ref$sig]), 3), "\n")

# ---- Part 3: dispersion ratio, age model vs main-run model ------------------------
library(DESeq2)
base    <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run     <- file.path(base, "final_runs/bystudy_20261015_151051")
out_dir <- file.path(base, "results/sensitivity_age_GSE79666")
dds_L   <- readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds"))
dds_age <- readRDS(file.path(out_dir, "dds_ModelL_GSE79666_age_adjusted.rds"))
stopifnot(identical(rownames(dds_L), rownames(dds_age)))
main_sig <- !is.na(results(dds_L, name = "condition_HD_vs_CTRL", alpha = 0.05)$padj) &
  results(dds_L, name = "condition_HD_vs_CTRL", alpha = 0.05)$padj < 0.05
ratio <- dispersions(dds_age) / dispersions(dds_L)
disp_tab <- data.frame(genes = c("all tested genes", "main DEGs"),
                       n = c(sum(is.finite(ratio)), sum(main_sig)),
                       median_ratio = round(c(median(ratio, na.rm = TRUE), median(ratio[main_sig])), 3),
                       mean_ratio   = round(c(mean(ratio, na.rm = TRUE),   mean(ratio[main_sig])), 3))
print(disp_tab)
write.csv(disp_tab, file.path(out_dir, "S10_dispersion_ratio_all_vs_DEGs.csv"), row.names = FALSE)
