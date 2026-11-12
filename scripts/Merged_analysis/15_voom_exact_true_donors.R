# =============================================================================
#  15_voom_exact_true_donors.R   (written 1 Nov 2026)
#  Why: script 14 reproduced the fold changes of the saved donor-aware limma-voom run
#  exactly, but not its numbers of significant genes. Cause: limma's contrasts.fit()
#  computes standard errors only APPROXIMATELY when the design is non-orthogonal and the
#  fit has precision weights (voom) - see ?contrasts.fit. The approximation depends on
#  how the design is parametrised, so the sensitivity run (~ batch + sex + group) and
#  script 14 (~ 0 + group + batch + sex) gave different p-values for the same model.
#  Fix: fit each comparison with a design in which that comparison IS a single
#  coefficient (same model, same column space, same voom weights and donor correlation),
#  which gives exact standard errors - as in the null calibration, where the tested
#  pseudo-group was a single coefficient.
#  Reruns with exact standard errors: donor-aware limma-voom (old donors, for comparison,
#  and true donors), limma-voom without donor blocking (Table 10 row), and everything
#  built on the donor-aware fit: high-confidence lists, the 137 most robust interaction
#  genes, core genes with high-confidence lists, the cross-study check of the
#  high-confidence interaction genes, and the striatum vs pooled cortex robustness.
#  NOT repeated (not affected by limma): composition scores and the without-juvenile
#  analysis from script 14 (truedonor_20261101_151301) stay valid.
#  Output: NEW folder Merged_analysis/final_runs/voomexact_<date>_<time>/
#  Run:    Session > Restart R, Source. About 5-10 minutes.
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(edgeR) })   # edgeR attaches limma
t0 <- Sys.time()
base <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run  <- file.path(base, "final_runs/bystudy_20261015_151051")
sens <- file.path(base, "final_runs/sensitivity_20261017_163558")
in_dir <- file.path(base, "results", "objects")
out_dir <- file.path(base, "final_runs", format(Sys.time(), "voomexact_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }
mins <- function() paste0("[", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min]")
slim <- function(d) { for (a in c("mu", "H", "replaceCounts", "replaceCooks")) assays(d)[[a]] <- NULL; d }
latest_dir <- function(prefix, must) {
  d <- list.dirs(file.path(base, "final_runs"), recursive = FALSE)
  d <- d[grepl(paste0("^", prefix, "_[0-9]"), basename(d)) & file.exists(file.path(d, must))]
  if (!length(d)) stop("No final_runs/", prefix, "_* folder containing ", must); sort(d)[length(d)] }
latest_file <- function(dir, pattern) { f <- list.files(dir, pattern = pattern, full.names = TRUE)
if (!length(f)) stop("No file matching ", pattern, " in ", dir); f[which.max(file.info(f)$mtime)] }

# True donor of every GEO sample (donor table of Paryani et al. 2024, Supplementary Data 1;
# batch 1-3: number in the sample title; batch 4: OF019-OF052 = table rows S43-S76)
key <- data.frame(
  title = c(
    "H4161_Accum", "H4161_Caud", "H4161_Cing", "C4812_Caud", "C4915_Accum", "C4915_Caud",
    "C4915_Cing", "H5109_Accum", "H5109_Caud", "H5109_Cing", "C5227_Accum", "C5227_Caud",
    "C5227_Cing", "H5263_Accum", "H5263_Caud", "H5263_Cing", "H5493_Accum", "H5493_Caud",
    "H5575_Accum", "C5596_Accum", "C5596_Caud", "C5596_Cing", "H5608_Accum", "H5608_Caud",
    "H5608_Cing", "C5700_Accum", "C5700_Caud", "C5700_Cing", "C5832_Caud", "C4812_Cing", "C98_Cing",
    "H3859_Cing", "H4929_Cing", "H5300_Cing", "H5493_Cing", "H5575_Caud", "H5575_Cing", "C180_Cing",
    "C4518_Cing", "C5404_Cing", "H5301_Cing", "C5382_Cing", "OF019", "OF020", "OF021", "OF022",
    "OF023", "OF024", "OF025", "OF026", "OF027", "OF028", "OF029", "OF030", "OF031", "OF032",
    "OF033", "OF034", "OF035", "OF036", "OF037", "OF038", "OF039", "OF040", "OF041", "OF042",
    "OF043", "OF044", "OF045", "OF046", "OF047", "OF048", "OF049", "OF050", "OF051", "OF052"),
  gsm = c(
    "GSM7754685", "GSM7754686", "GSM7754687", "GSM7754688", "GSM7754689", "GSM7754690",
    "GSM7754691", "GSM7754692", "GSM7754693", "GSM7754694", "GSM7754695", "GSM7754696",
    "GSM7754697", "GSM7754698", "GSM7754699", "GSM7754700", "GSM7754701", "GSM7754702",
    "GSM7754703", "GSM7754704", "GSM7754705", "GSM7754706", "GSM7754707", "GSM7754708",
    "GSM7754709", "GSM7754710", "GSM7754711", "GSM7754712", "GSM7754713", "GSM7754714",
    "GSM7754715", "GSM7754716", "GSM7754717", "GSM7754718", "GSM7754719", "GSM7754720",
    "GSM7754721", "GSM7754722", "GSM7754723", "GSM7754724", "GSM7754725", "GSM7754726",
    "GSM7754727", "GSM7754728", "GSM7754729", "GSM7754730", "GSM7754731", "GSM7754732",
    "GSM7754733", "GSM7754734", "GSM7754735", "GSM7754736", "GSM7754737", "GSM7754738",
    "GSM7754739", "GSM7754740", "GSM7754741", "GSM7754742", "GSM7754743", "GSM7754744",
    "GSM7754745", "GSM7754746", "GSM7754747", "GSM7754748", "GSM7754749", "GSM7754750",
    "GSM7754751", "GSM7754752", "GSM7754753", "GSM7754754", "GSM7754755", "GSM7754756",
    "GSM7754757", "GSM7754758", "GSM7754759", "GSM7754760"),
  srr = c(
    "SRR25877173", "SRR25877171", "SRR25877170", "SRR25877168", "SRR25877167", "SRR25877165",
    "SRR25877164", "SRR25877162", "SRR25877161", "SRR25877160", "SRR25877159", "SRR25877158",
    "SRR25877157", "SRR25877156", "SRR25877155", "SRR25877154", "SRR25877150", "SRR25877151",
    "SRR25877149", "SRR25877146", "SRR25877142", "SRR25877143", "SRR25877141", "SRR25877140",
    "SRR25877139", "SRR25877138", "SRR25877137", "SRR25877136", "SRR25877135", "SRR25877134",
    "SRR25877133", "SRR25877132", "SRR25877107", "SRR25877105", "SRR25877104", "SRR25877102",
    "SRR25877101", "SRR25877099", "SRR25877098", "SRR25877097", "SRR25877096", "SRR25877095",
    "SRR25877094", "SRR25877093", "SRR25877092", "SRR25877091", "SRR25877090", "SRR25877089",
    "SRR25877085", "SRR25877086", "SRR25877084", "SRR25877081", "SRR25877077", "SRR25877078",
    "SRR25877076", "SRR25877075", "SRR25877074", "SRR25877073", "SRR25877072", "SRR25877071",
    "SRR25877070", "SRR25877069", "SRR25877068", "SRR25877063", "SRR25877131", "SRR25877128",
    "SRR25877127", "SRR25877126", "SRR25877125", "SRR25877124", "SRR25877123", "SRR25877122",
    "SRR25877121", "SRR25877120", "SRR25877119", "SRR25877118"),
  true_donor = c(
    "T4161", "T4161", "T4161", "T4812", "T4915", "T4915", "T4915", "T5109", "T5109", "T5109",
    "T5227", "T5227", "T5227", "T5263", "T5263", "T5263", "T5493", "T5493", "T5575", "T5596",
    "T5596", "T5596", "T5608", "T5608", "T5608", "T5700", "T5700", "T5700", "T5832", "T4812", "T98",
    "T3859", "T4929", "T5300", "T5493", "T5575", "T5575", "T180", "T4518", "T5404", "T5301",
    "T5382", "T178", "T1905", "T3859", "T4410", "T4410", "T4273", "T4273", "T4285", "T4285",
    "T4285", "T4501", "T4501", "T4638", "T4638", "T4638", "T5263", "T5266", "T5266", "T5266",
    "T5286", "T5286", "T5534", "T5534", "T5534", "T5693", "T5693", "T5693", "T5714", "T5714",
    "T5714", "T5775", "T5775", "T5798", "T5798"),
  stringsAsFactors = FALSE)

# ---- 1. Samples, true donors and juvenile flag ---------------------------------------
ddsP <- slim(readRDS(file.path(run, "objects", "dds_ModelP_Paryani.rds"))); invisible(gc())
cd <- as.data.frame(colData(ddsP))
idx <- ifelse(cd$sample_id %in% key$title, match(cd$sample_id, key$title),
              ifelse(cd$sample_id %in% key$gsm, match(cd$sample_id, key$gsm), match(cd$sample_id, key$srr)))
stopifnot(!anyNA(idx))
cd$title <- key$title[idx]; cd$true_donor <- key$true_donor[idx]
old_donor  <- as.character(cd$donor); true_donor <- cd$true_donor
split_true <- names(which(tapply(old_donor, true_donor, function(x) length(unique(x))) > 1))
say("Paryani samples: ", nrow(cd), " | donors old: ", length(unique(old_donor)), ", true: ", length(unique(true_donor)))
say("True donors that were split in the analysis: ", paste(split_true, collapse = ", "),
    " (samples: ", paste(cd$title[true_donor %in% split_true], collapse = ", "), ")")
sh <- cd$group == "Striatum_HD"
say("Striatal HD samples (null calibration): donor assignment identical to true donors: ",
    all(tapply(old_donor[sh], true_donor[sh], function(x) length(unique(x))) == 1) &&
      all(tapply(true_donor[sh], old_donor[sh], function(x) length(unique(x))) == 1))
jv_old <- toupper(as.character(cd$juvenile_onset)) %in% c("TRUE", "Y", "YES", "1", "J") | toupper(as.character(cd$grade)) == "J"
jv_new <- jv_old | true_donor == "T3859"
say("Juvenile-onset samples: old flag ", sum(jv_old), " (expected 16), corrected ", sum(jv_new),
    " (added: ", paste(cd$title[jv_new & !jv_old], collapse = ", "), ")")

# ---- 2. DESeq2 results of the main model (unchanged; same rule as the main run) ------------
rn <- resultsNames(ddsP); ref <- levels(ddsP$group)[1]
g <- function(lv, rnx = rn, refx = ref) { v <- setNames(numeric(length(rnx)), rnx); if (lv == refx) return(v)
v[paste0("group_", lv, "_vs_", refx)] <- 1; v }
tab_P <- function(dds, vec, groups) {
  r <- results(dds, contrast = vec, alpha = 0.05); cnt <- counts(dds)
  nz <- function(lv) rowSums(cnt[, dds$group == lv, drop = FALSE]) > 0
  rel <- Reduce(`&`, lapply(groups, nz)) & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10)
  data.frame(gene_id = rownames(r), b = r$log2FoldChange, se = r$lfcSE, padj = r$padj, estimable = rel,
             tested = rel & !is.na(r$padj), sig = rel & !is.na(r$padj) & r$padj < 0.05)
}
SH <- "Striatum_HD"; SC <- "Striatum_CTRL"; CH <- "Cortex_HD"; CC <- "Cortex_CTRL"
cvec <- function(k, rnx = rn, refx = ref) {
  G <- function(lv) g(lv, rnx, refx)
  switch(k, A = G(SH) - G(SC), B_P = G(CH) - G(CC), C = G(SC) - G(CC), D = G(SH) - G(CH),
         E_P = G(SH) - G(SC) - G(CH) + G(CC)) }
cgrp <- list(A = c(SH, SC), B_P = c(CH, CC), C = c(SC, CC), D = c(SH, CH), E_P = c(SH, SC, CH, CC))
K <- names(cgrp)
lab <- c(A = "Striatum, HD vs control", B_P = "Cingulate cortex, HD vs control", C = "Controls, striatum vs cingulate",
         D = "HD, striatum vs cingulate", E_P = "Region-dependent HD effect (interaction)")
des <- lapply(setNames(K, K), function(k) tab_P(ddsP, cvec(k), cgrp[[k]]))
say("Self-check DESeq2 (main run: 6267, 1000, 5811, 12115, 589): ", paste(sapply(des, function(d) sum(d$sig)), collapse = ", "))
gene_ids <- rownames(ddsP); samples <- colnames(ddsP)
rm(ddsP); invisible(gc())

ddsL <- slim(readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds")))
rL <- results(ddsL, name = "condition_HD_vs_CTRL", alpha = 0.05); cL <- counts(ddsL)
relL <- rowSums(cL[, ddsL$condition == "HD"]) > 0 & rowSums(cL[, ddsL$condition == "CTRL"]) > 0 &
  !(!is.na(rL$log2FoldChange) & abs(rL$log2FoldChange) > 10)
BT <- data.frame(gene_id = rownames(rL), b = rL$log2FoldChange, se = rL$lfcSE, padj = rL$padj,
                 tested = relL & !is.na(rL$padj), sig = relL & !is.na(rL$padj) & rL$padj < 0.05)
rm(ddsL, rL, cL); invisible(gc())
say("Self-check motor cortex (714): ", sum(BT$sig), " ", mins())

# ---- 3. Donor-aware limma-voom, old and true donors ------------------------------------
txi <- readRDS(latest_file(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$"))
stopifnot(identical(txi$countsFromAbundance, "no"), all(samples %in% colnames(txi$counts)), all(gene_ids %in% rownames(txi$counts)))
ab <- txi$abundance[, samples]; ln <- txi$length[, samples]; ct <- txi$counts[, samples]; rm(txi); invisible(gc())
nc  <- ab * rowMeans(ln)                                   # lengthScaledTPM counts (as in the sensitivity run)
lst <- t(t(nc) * (colSums(ct) / colSums(nc)))[gene_ids, ]
rm(ab, ln, ct, nc); invisible(gc())
y <- calcNormFactors(DGEList(counts = lst))                  # TMM
grp <- factor(as.character(cd$group), levels = c(SC, SH, CC, CH))
cov_df <- data.frame(batch_study = factor(as.character(cd$batch_study)), sex = factor(as.character(cd$sex)))
dmat <- model.matrix(~ batch_study + sex + grp, data = cbind(cov_df, grp = grp))   # base design (as in the sensitivity run)
# one design per comparison: same column space, the comparison is a single coefficient
design_k <- function(k) {
  if (k == "E_P") {
    tis  <- factor(ifelse(grepl("^Striatum", grp), "Striatum", "Cortex"), levels = c("Cortex", "Striatum"))
    cond <- factor(ifelse(grepl("_HD$", grp), "HD", "CTRL"), levels = c("CTRL", "HD"))
    mm <- model.matrix(~ batch_study + sex + tis * cond, data = cbind(cov_df, tis = tis, cond = cond))
    return(list(mm = mm, coef = "tisStriatum:condHD")) }
  ref <- c(A = SC, B_P = CC, C = CC, D = CH)[[k]]; tgt <- c(A = SH, B_P = CH, C = SC, D = SH)[[k]]
  g2 <- relevel(grp, ref = ref)
  list(mm = model.matrix(~ batch_study + sex + g2, data = cbind(cov_df, g2 = g2)), coef = paste0("g2", tgt)) }
for (k in K) { dk <- design_k(k); stopifnot(dk$coef %in% colnames(dk$mm), qr(dk$mm)$rank == qr(dmat)$rank, ncol(dk$mm) == ncol(dmat)) }
cm_base <- matrix(0, ncol(dmat), 5, dimnames = list(colnames(dmat), K))       # contrasts exactly as in the sensitivity run
cm_base["grpStriatum_HD", "A"] <- 1; cm_base[c("grpCortex_HD", "grpCortex_CTRL"), "B_P"] <- c(1, -1)
cm_base["grpCortex_CTRL", "C"] <- -1; cm_base[c("grpStriatum_HD", "grpCortex_HD"), "D"] <- c(1, -1)
cm_base[c("grpStriatum_HD", "grpCortex_HD", "grpCortex_CTRL"), "E_P"] <- c(1, -1, 1)
fit_exact <- function(v, block = NULL, correlation = NULL) lapply(setNames(K, K), function(k) {
  dk <- design_k(k)
  fit <- if (is.null(block)) eBayes(lmFit(v, dk$mm)) else eBayes(lmFit(v, dk$mm, block = block, correlation = correlation))
  tt <- topTable(fit, coef = dk$coef, number = Inf, sort.by = "none")
  data.frame(gene_id = rownames(tt), b = tt$logFC, se = fit$stdev.unscaled[, dk$coef] * sqrt(fit$s2.post), padj = tt$adj.P.Val) })
run_voom <- function(block, label) {
  say("limma-voom with ", label, " (", length(unique(block)), " donors)... ", mins())
  v   <- voom(y, dmat)
  dc1 <- duplicateCorrelation(v, dmat, block = block)
  v2  <- voom(y, dmat, block = block, correlation = dc1$consensus)
  dc2 <- duplicateCorrelation(v2, dmat, block = block)
  say("  consensus within-donor correlation: ", round(dc2$consensus, 3), " ", mins())
  fa <- eBayes(contrasts.fit(lmFit(v2, dmat, block = block, correlation = dc2$consensus), cm_base))
  approx <- lapply(setNames(K, K), function(k) data.frame(gene_id = rownames(fa), padj = p.adjust(fa$p.value[, k], "BH")))
  list(res = fit_exact(v2, block, dc2$consensus), approx = approx, rho = dc2$consensus)
}
say("limma-voom without donor blocking (exact)... ", mins())
V_nb   <- list(res = fit_exact(voom(y, dmat)))
V_old  <- run_voom(old_donor,  "the donors used so far (exact SE)")
V_true <- run_voom(true_donor, "the true donors (exact SE)")
rm(y, lst); invisible(gc())
for (k in K) write.csv(V_true$res[[k]], file.path(out_dir, paste0("voom_donor_TRUE_exact_", k, ".csv")), row.names = FALSE)

# comparison with the saved sensitivity run (approximate SE, old donors)
vd <- read.csv(file.path(sens, "tables", "SENS_voom_donor_results.csv"))
rel <- function(k) des[[k]]$estimable       # DESeq2 estimability rule, as the 'reliable' flag of the sensitivity run
rep_chk <- do.call(rbind, lapply(K, function(k) { v <- V_old$res[[k]][match(des[[k]]$gene_id, V_old$res[[k]]$gene_id), ]
s <- vd[match(des[[k]]$gene_id, vd$gene_id), ]
data.frame(contrast = lab[[k]], rho_log2FC_vs_saved = round(cor(v$b, s[[paste0("log2FC_", k)]], method = "spearman", use = "complete.obs"), 4),
           median_SE_ratio_exact_vs_saved = round(median(v$se / s[[paste0("lfcSE_", k)]], na.rm = TRUE), 3),
           DEGs_saved_approx = sum(s[[paste0("sig_", k)]] %in% TRUE),
           DEGs_recomputed_approx = sum(rel(k) & (V_old$approx[[k]]$padj[match(des[[k]]$gene_id, V_old$approx[[k]]$gene_id)] < 0.05) %in% TRUE),
           DEGs_exact_old_donors = sum(rel(k) & (v$padj < 0.05) %in% TRUE)) }))
show("Donor-aware limma-voom, old donors: saved run vs recomputed with the same approximate SE vs exact SE (DEGs among estimable genes)", rep_chk)
say("If DEGs_saved_approx = DEGs_recomputed_approx, the difference is fully explained by the contrasts.fit approximation.")

# ---- 4. High-confidence lists --------------------------------------------------------
hc_of <- function(V, k) { d <- des[[k]]; v <- V$res[[k]][match(d$gene_id, V$res[[k]]$gene_id), ]
d$sig & (v$padj < 0.05) %in% TRUE & (sign(v$b) == sign(d$b)) %in% TRUE }
hc_saved <- function(k) { d <- des[[k]]; s <- vd[match(d$gene_id, vd$gene_id), ]
d$sig & s[[paste0("sig_", k)]] %in% TRUE & (sign(s[[paste0("log2FC_", k)]]) == sign(d$b)) %in% TRUE }
HC <- lapply(setNames(K, K), function(k) hc_of(V_true, k))
rho_k <- function(V, k) { d <- des[[k]]; v <- V$res[[k]][match(d$gene_id, V$res[[k]]$gene_id), ]
round(cor(d$b[d$tested], v$b[d$tested], method = "spearman", use = "complete.obs"), 2) }
hc_tab <- do.call(rbind, lapply(K, function(k) data.frame(contrast = lab[[k]], DESeq2_DEGs = sum(des[[k]]$sig),
                                                          HC_saved_run_approx = sum(hc_saved(k)), HC_exact_old_donors = sum(hc_of(V_old, k)), HC_exact_true_donors = sum(HC[[k]]),
                                                          HC_true_higher = sum(HC[[k]] & des[[k]]$b > 0), HC_true_lower = sum(HC[[k]] & des[[k]]$b < 0),
                                                          pct_of_DEGs_true = round(100 * sum(HC[[k]]) / sum(des[[k]]$sig), 1),
                                                          rho_DESeq2_vs_voom_true = rho_k(V_true, k),
                                                          in_both_old_and_true = sum(HC[[k]] & hc_of(V_old, k)))))
show("High-confidence genes (DESeq2 and donor-aware limma-voom, same direction)", hc_tab)
write.csv(hc_tab, file.path(out_dir, "HighConfidence_counts_exact_true_donors.csv"), row.names = FALSE)
for (k in K) write.csv(data.frame(gene_id = des[[k]]$gene_id, log2FC_DESeq2 = round(des[[k]]$b, 3),
                                  padj_DESeq2 = signif(des[[k]]$padj, 3))[HC[[k]], ],
                       file.path(out_dir, paste0("HighConfidence_genes_exact_true_donors_", k, ".csv")), row.names = FALSE)

# Table 10 rows: limma-voom without donor blocking and donor-aware (true donors), exact SE
ret_row <- function(V, label) do.call(cbind, c(list(data.frame(analysis = label)), lapply(c("A", "B_P", "E_P"), function(k) {
  d <- des[[k]]; v <- V$res[[k]][match(d$gene_id, V$res[[k]]$gene_id), ]
  x <- data.frame(paste0(round(100 * sum(d$sig & (v$padj < 0.05) %in% TRUE & (sign(v$b) == sign(d$b)) %in% TRUE) / sum(d$sig), 1),
                         " (", rho_k(V, k), ")")); names(x) <- k; x })))
t10_tab <- rbind(ret_row(V_nb, "limma-voom, donors independent (thesis: 85.9 (0.96), 70.1 (0.97), 72.3 (0.93))"),
                 ret_row(V_true, "limma-voom, donor-aware, true donors (thesis: 72.2 (0.95), 43.0 (0.96), 89.6 (0.92))"))
show("Table 10 rows with exact standard errors: % of main DEGs retained (Spearman rho)", t10_tab)
write.csv(t10_tab, file.path(out_dir, "Table10_voom_rows_exact.csv"), row.names = FALSE)

# ---- 5. Most robust interaction genes (high-confidence AND without batch 4) ---------------
hcdir <- latest_dir("highconf", "Interaction_genes_robustness_flags.csv")
fl <- read.csv(file.path(hcdir, "Interaction_genes_robustness_flags.csv"))
wo4_ids <- fl$gene_id[fl$significant_without_batch4 %in% TRUE]
EP <- des$E_P
say("\nInteraction genes significant without batch 4 (from ", basename(hcdir), "): ", length(wo4_ids), " (thesis 141)")
say("Interaction genes high-confidence AND without batch 4: old ", sum(fl$high_confidence & fl$significant_without_batch4),
    " (thesis 137), true donors ", sum(HC$E_P & EP$gene_id %in% wo4_ids))

# ---- 6. Core genes with high-confidence lists ----------------------------------------------
core <- function(sA, sB) {
  ids <- Reduce(intersect, list(des$A$gene_id[sA], des$B_P$gene_id[sB], BT$gene_id[BT$sig]))
  a <- des$A[match(ids, des$A$gene_id), ]; b <- des$B_P[match(ids, des$B_P$gene_id), ]; t <- BT[match(ids, BT$gene_id), ]
  ids[sign(a$b) == sign(b$b) & sign(a$b) == sign(t$b)] }
core_main <- core(des$A$sig, des$B_P$sig); core_hc <- core(HC$A, HC$B_P)
say("Core genes: main ", length(core_main), " (thesis 40); with high-confidence lists, true donors ", length(core_hc),
    " (thesis 22), all among the 40: ", all(core_hc %in% core_main))

# ---- 7. Cross-study check of the stricter interaction sets ---------------------------------
A <- des$A; BP <- des$B_P
ids <- intersect(A$gene_id[A$tested], BT$gene_id[BT$tested])
a <- A[match(ids, A$gene_id), ]; t <- BT[match(ids, BT$gene_id), ]; cc <- BP[match(ids, BP$gene_id), ]; e <- EP[match(ids, EP$gene_id), ]
d_ST <- a$b - t$b; sig_ST <- p.adjust(2 * pnorm(-abs(d_ST / sqrt(a$se^2 + t$se^2))), "BH") < 0.05
closer <- abs(t$b - cc$b) < abs(t$b - a$b)
hcE_ids <- EP$gene_id[HC$E_P]
check_set <- function(name, inset) {
  inset <- inset & !is.na(inset)
  ft <- fisher.test(table(factor(inset, c(FALSE, TRUE)), factor(sig_ST, c(FALSE, TRUE))))
  data.frame(gene_set = name, genes_tested_in_both = sum(inset), pct_differ_from_motor_cortex = round(100 * sum(inset & sig_ST) / sum(inset), 1),
             other_genes_pct = round(100 * sum(!inset & sig_ST) / sum(!inset), 1), odds_ratio = round(unname(ft$estimate), 1),
             motor_cortex_closer_to_cingulate_pct = round(100 * sum(inset & closer) / sum(inset), 1)) }
reg <- rbind(check_set("All interaction genes (thesis 584; 64.9%)", e$sig),
             check_set("High-confidence, true donors (thesis 525; 67.0%)", ids %in% hcE_ids),
             check_set("High-confidence + without batch 4, true donors (thesis 136; 77.2%)", ids %in% intersect(hcE_ids, wo4_ids)))
show("Cross-study regional check (striatal vs motor-cortex HD effect)", reg)
write.csv(reg, file.path(out_dir, "Regional_crossstudy_exact_true_donors.csv"), row.names = FALSE)

# ---- 8. Striatum vs pooled cortex: robustness with the exact true-donor voom -------------------------
pcdir <- latest_dir("pooled_cortex", "PooledCortex_all_genes.rds")
prdir <- latest_dir("pooled_robustness", "PooledRobustness_region_dependent_genes_flags.csv")
Gp <- readRDS(file.path(pcdir, "PooledCortex_all_genes.rds"))
pf <- read.csv(file.path(prdir, "PooledRobustness_region_dependent_genes_flags.csv"))
pid <- Gp$gene_id; main_sig <- Gp$region_difference_padj < 0.05; d_main <- Gp$striatum_minus_pooled_cortex
tb <- BT$b[match(pid, BT$gene_id)]; tse <- BT$se[match(pid, BT$gene_id)]
pool_diff <- function(V) {
  s <- V$res$A[match(pid, V$res$A$gene_id), ]; c <- V$res$B_P[match(pid, V$res$B_P$gene_id), ]; e <- V$res$E_P[match(pid, V$res$E_P$gene_id), ]
  a <- (1 / c$se^2) / (1 / c$se^2 + 1 / tse^2)
  b_cx <- a * c$b + (1 - a) * tb; se_cx <- sqrt(a^2 * c$se^2 + (1 - a)^2 * tse^2)
  cov_SC <- (s$se^2 + c$se^2 - e$se^2) / 2; d <- s$b - b_cx; v <- s$se^2 + se_cx^2 - 2 * a * cov_SC
  ok <- is.finite(d) & is.finite(v) & v > 0; p <- rep(NA_real_, length(d)); p[ok] <- 2 * pnorm(-abs(d[ok] / sqrt(v[ok])))
  list(d = d, padj = p.adjust(p, "BH")) }
robust <- function(x) main_sig & (x$padj < 0.05) %in% TRUE & (sign(x$d) == sign(d_main)) %in% TRUE
wo4p <- pid %in% pf$gene_id[pf$significant_without_batch4 %in% TRUE]
hc_old <- robust(pool_diff(V_old)); hc_new <- robust(pool_diff(V_true))
neg <- d_main < 0
prow <- function(lab, x) data.frame(gene_set = lab, genes = sum(x), more_negative_in_striatum = sum(x & neg),
                                    more_positive_in_striatum = sum(x & !neg), pct_of_1867 = round(100 * sum(x) / sum(main_sig), 1))
ptab <- rbind(prow("All region-dependent genes", main_sig),
              prow("High-confidence, old donors, exact (thesis 1,469)", hc_old), prow("High-confidence, true donors, exact", hc_new),
              prow("Significant without batch 4 (unchanged)", main_sig & wo4p),
              prow("High-confidence + without batch 4, old donors, exact (thesis 973)", hc_old & wo4p),
              prow("High-confidence + without batch 4, true donors, exact", hc_new & wo4p))
show("Striatum vs pooled cortex genes", ptab)
write.csv(ptab, file.path(out_dir, "PooledRobustness_exact_true_donors.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = pid, gene_symbol = Gp$gene_symbol, striatum_minus_pooled_cortex = d_main,
                     high_confidence_true_donors = hc_new, significant_without_batch4 = wo4p)[main_sig, ],
          file.path(out_dir, "PooledRobustness_flags_exact_true_donors.csv"), row.names = FALSE)
say(mins())

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone ", mins(), ". Send Claude the console output (or run_log.txt). Outputs: ", out_dir)