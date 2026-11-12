# =============================================================================
# 11b_pooled_cortex_heterogeneity.R
#
# Heterogeneity between the two cortical HD effects and a heterogeneity-adjusted
# sensitivity analysis of the pooled cortex.
#
#   1. lambda_het: inflation factor of the per-gene heterogeneity statistics
#      z_het = (b_cing - b_motor) / sqrt(SE_cing^2 + SE_motor^2)
#      over all genes that passed the count filter in both models (B_P, B_T).
#   2. Pooled cortex (fixed-effect inverse-variance meta-analysis) and the
#      striatum-minus-pooled-cortex comparison, with SE_cx multiplied by
#      phi = 1 (main) or phi = sqrt(lambda_het) (heterogeneity-adjusted).
#
# Reproduces (bystudy_20261015_151051):
#   lambda_heterogeneity_cingulate_vs_motor.csv  -> 23,730 genes, lambda_het 1.263
#   pooled_cortex_heterogeneity_sensitivity.csv  -> main 3,374 / 1,867;
#                                                   adjusted 1,567 / 1,523
# The script stops if any of these numbers is not reproduced.
#
# Input : <run>/objects/results_by_contrast.rds (04_main_DESeq2_by_study2.R)
#         contrasts A (striatum HD vs CTRL), B_P (cingulate HD vs CTRL),
#         B_T (motor cortex HD vs CTRL), E_P (interaction); unshrunken log2FC
# Output: new time-stamped folder final_runs/heterogeneity_<date_time>/
# =============================================================================

run  <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis/final_runs/bystudy_20261015_151051"
out  <- file.path(dirname(run), paste0("heterogeneity_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out, recursive = TRUE)

res_all <- readRDS(file.path(run, "objects", "results_by_contrast.rds"))
A  <- res_all[["A"]]      # striatum HD vs control (Paryani)
P  <- res_all[["B_P"]]    # cingulate HD vs control (Paryani)
Tm <- res_all[["B_T"]]    # motor cortex HD vs control (Lin / GSE79666)
I  <- res_all[["E_P"]]    # interaction: striatal minus cingulate HD effect

# genes with a reliable estimate and a defined padj in a contrast
est <- function(x) rownames(x)[x$reliable %in% TRUE & !is.na(x$padj)]

# ---- 1. lambda_het over all genes that passed the count filter in both models -
# (= all genes with a finite log2FC and SE in both cortical contrasts; 23,730).
# Stricter gene definitions give lambda_het 1.263-1.266 (checked 2026-11-10).
zhet <- function(gg) (P[gg, "log2FC"] - Tm[gg, "log2FC"]) /
  sqrt(P[gg, "lfcSE"]^2 + Tm[gg, "lfcSE"]^2)
g_cx       <- intersect(rownames(P), rownames(Tm))
z_het      <- zhet(g_cx)
dropped    <- sum(!is.finite(z_het))
lambda_het <- median(z_het[is.finite(z_het)]^2) / qchisq(0.5, df = 1)

# ---- 2. Pooled cortex and striatum vs pooled cortex, SE_cx scaled by phi -----
g <- Reduce(intersect, lapply(list(A, P, Tm, I), est))   # genes estimable in all four

# per-gene heterogeneity test on the 22,371 pooling genes (as in 11_pooled_cortex.R)
padj_het <- p.adjust(2 * pnorm(-abs(zhet(g))), "BH")

lam_tab <- data.frame(genes = length(g_cx), dropped = dropped, lambda_het = lambda_het,
                      pooling_genes = length(g),
                      pooling_genes_heterogeneous_padj05 = sum(padj_het < 0.05))
print(lam_tab)

pool <- function(phi) {
  a    <- (1 / P[g, "lfcSE"]^2) / (1 / P[g, "lfcSE"]^2 + 1 / Tm[g, "lfcSE"]^2)
  bcx  <- a * P[g, "log2FC"] + (1 - a) * Tm[g, "log2FC"]
  secx <- sqrt(a^2 * P[g, "lfcSE"]^2 + (1 - a)^2 * Tm[g, "lfcSE"]^2) * phi
  covs <- (A[g, "lfcSE"]^2 + P[g, "lfcSE"]^2 - I[g, "lfcSE"]^2) / 2   # Cov(b_str, b_cing)
  d    <- A[g, "log2FC"] - bcx
  vd   <- A[g, "lfcSE"]^2 + secx^2 - 2 * a * covs
  data.frame(gene_id = g, weight_cingulate = a,
             pooled_cortex_log2FC = bcx, pooled_cortex_SE = secx,
             pooled_cortex_padj = p.adjust(2 * pnorm(-abs(bcx / secx)), "BH"),
             striatum_minus_pooled_cortex = d, region_difference_SE = sqrt(vd),
             region_difference_padj = p.adjust(2 * pnorm(-abs(d / sqrt(vd))), "BH"))
}
m <- pool(1)
# The published analysis widened SE_cx by sqrt of lambda_het rounded to 3 decimals
# (1.263), as reported in the text; the unrounded value is shown for comparison.
lambda_used <- round(lambda_het, 3)
h <- pool(sqrt(lambda_used))
h_exact <- pool(sqrt(lambda_het))
cat(sprintf("\nUnrounded lambda_het %.5f: %d pooled DEGs, %d region-dependent genes\n",
            lambda_het, sum(h_exact$pooled_cortex_padj < 0.05),
            sum(h_exact$region_difference_padj < 0.05)))

sig_m_pool <- m$gene_id[m$pooled_cortex_padj < 0.05]
sig_h_pool <- h$gene_id[h$pooled_cortex_padj < 0.05]
sig_m_reg  <- m$gene_id[m$region_difference_padj < 0.05]
sig_h_reg  <- h$gene_id[h$region_difference_padj < 0.05]

sens_tab <- data.frame(
  analysis         = c("main", "heterogeneity_adjusted"),
  SE_multiplier    = c(1, sqrt(lambda_used)),
  genes            = length(g),
  pooled_DEGs      = c(length(sig_m_pool), length(sig_h_pool)),
  region_dependent = c(length(sig_m_reg),  length(sig_h_reg)),
  pooled_retained_pct           = c(100, round(100 * mean(sig_m_pool %in% sig_h_pool), 1)),
  region_dependent_retained_pct = c(100, round(100 * mean(sig_m_reg  %in% sig_h_reg),  1)))
print(sens_tab)

# ---- 3. Checks against the published numbers ---------------------------------
stopifnot(
  length(g_cx) == 23730, dropped == 0, abs(lambda_het - 1.263) < 5e-4,
  length(g) == 22371, sum(padj_het < 0.05) == 2,
  length(sig_m_pool) == 3374, length(sig_m_reg) == 1867,
  length(sig_h_pool) == 1567, length(sig_h_reg) == 1523,
  all(sig_h_pool %in% sig_m_pool), all(sig_h_reg %in% sig_m_reg)
)
cat("\nAll published numbers reproduced.\n")

# ---- 4. Save ------------------------------------------------------------------
write.csv(lam_tab,  file.path(out, "lambda_heterogeneity_cingulate_vs_motor.csv"), row.names = FALSE)
write.csv(sens_tab, file.path(out, "pooled_cortex_heterogeneity_sensitivity.csv"), row.names = FALSE)

per_gene <- merge(m, h[, c("gene_id", "pooled_cortex_SE", "pooled_cortex_padj",
                           "region_difference_SE", "region_difference_padj")],
                  by = "gene_id", suffixes = c("_main", "_het_adjusted"))
write.csv(per_gene, file.path(out, "pooled_cortex_per_gene_main_and_het_adjusted.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = sig_h_reg),
          file.path(out, "region_dependent_genes_heterogeneity_adjusted.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Outputs written to:", out, "\n")