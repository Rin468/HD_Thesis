# =============================================================================
#  09_high_confidence_DEGs.R   (written 1 Nov 2026)
#  Why: the null calibration (nullcal_20261031_225350) showed that methods which
#  treat samples from the same donor as independent give false positives
#  (DESeq2: lambda 1.45, false DEGs in every donor-level split), whereas
#  donor-aware limma-voom was calibrated (lambda 1.02).
#  What: "high-confidence" genes = significant (padj < 0.05) in BOTH the main DESeq2
#  model and the donor-aware limma-voom model, with the same direction.
#    1. Counts for every Paryani contrast (striatum, cingulate, tissue differences,
#       interaction).
#    2. Interaction genes that also survive without the HD-only batch 4.
#    3. Core genes (all three regions) recomputed with the high-confidence lists.
#    4. Cross-study regional check (as in 08) on the high-confidence interaction genes.
#  Nothing is refitted except, if needed, the model without batch 4 (a few minutes).
#
#  Output: NEW folder Merged_analysis/final_runs/highconf_<date>_<time>/
#  Run:    Session > Restart R, open this file, Source.
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2) })

base     <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run      <- file.path(base, "final_runs/bystudy_20261015_151051")
sens     <- file.path(base, "final_runs/sensitivity_20261017_163558")
out_dir  <- file.path(base, "final_runs", format(Sys.time(), "highconf_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }

# ---- 1. Main models and DESeq2 results (same rule as the main run) --------------
ddsP <- readRDS(file.path(run, "objects", "dds_ModelP_Paryani.rds"))
ddsL <- readRDS(file.path(run, "objects", "dds_ModelL_GSE79666.rds"))
rn <- resultsNames(ddsP); ref <- levels(ddsP$group)[1]
g <- function(lv) { v <- setNames(numeric(length(rn)), rn); if (lv == ref) return(v)
nm <- paste0("group_", lv, "_vs_", ref); stopifnot(nm %in% rn); v[nm] <- 1; v }
tab_P <- function(dds, vec, groups) {
  r <- results(dds, contrast = vec, alpha = 0.05); cnt <- counts(dds)
  nz <- function(lv) rowSums(cnt[, dds$group == lv, drop = FALSE]) > 0
  rel <- Reduce(`&`, lapply(groups, nz)) & !(!is.na(r$log2FoldChange) & abs(r$log2FoldChange) > 10)
  data.frame(gene_id = rownames(r), b = r$log2FoldChange, se = r$lfcSE, padj = r$padj,
             tested = rel & !is.na(r$padj), sig = rel & !is.na(r$padj) & r$padj < 0.05)
}
SH <- "Striatum_HD"; SC <- "Striatum_CTRL"; CH <- "Cortex_HD"; CC <- "Cortex_CTRL"
contr <- list(A   = list(g(SH) - g(SC), c(SH, SC)),
              B_P = list(g(CH) - g(CC), c(CH, CC)),
              C   = list(g(SC) - g(CC), c(SC, CC)),
              D   = list(g(SH) - g(CH), c(SH, CH)),
              E_P = list(g(SH) - g(SC) - g(CH) + g(CC), c(SH, SC, CH, CC)))
lab <- c(A = "Striatum, HD vs control", B_P = "Cingulate cortex, HD vs control",
         C = "Controls, striatum vs cingulate", D = "HD, striatum vs cingulate",
         E_P = "Region-dependent HD effect (interaction)")
des <- lapply(contr, function(z) tab_P(ddsP, z[[1]], z[[2]]))

rL <- results(ddsL, name = "condition_HD_vs_CTRL", alpha = 0.05); cL <- counts(ddsL)
relL <- rowSums(cL[, ddsL$condition == "HD"]) > 0 & rowSums(cL[, ddsL$condition == "CTRL"]) > 0 &
  !(!is.na(rL$log2FoldChange) & abs(rL$log2FoldChange) > 10)
BT <- data.frame(gene_id = rownames(rL), b = rL$log2FoldChange, se = rL$lfcSE, padj = rL$padj,
                 tested = relL & !is.na(rL$padj), sig = relL & !is.na(rL$padj) & rL$padj < 0.05)
say("Self-check (main run: A 6267, B_P 1000, C 5811, D 12115, E_P 589, motor cortex 714): ",
    paste(names(des), sapply(des, function(d) sum(d$sig)), sep = " ", collapse = ", "),
    ", motor cortex ", sum(BT$sig))

# ---- 2. Donor-aware limma-voom results (saved per gene by the sensitivity run) ----
vd <- read.csv(file.path(sens, "tables", "SENS_voom_donor_results.csv"))
sym <- setNames(vd$gene_symbol, vd$gene_id)

hc_tab <- list(); counts_tab <- list()
for (k in names(contr)) {
  d <- des[[k]]
  v <- vd[match(d$gene_id, vd$gene_id), ]
  v_sig <- v[[paste0("sig_", k)]] %in% TRUE
  same  <- sign(d$b) == sign(v[[paste0("log2FC_", k)]])
  hc    <- d$sig & v_sig & same %in% TRUE
  hc_tab[[k]] <- data.frame(gene_id = d$gene_id, gene_symbol = unname(sym[d$gene_id]),
                            log2FC_DESeq2 = round(d$b, 3), padj_DESeq2 = signif(d$padj, 3),
                            log2FC_voom_donor = round(v[[paste0("log2FC_", k)]], 3),
                            padj_voom_donor = signif(v[[paste0("padj_", k)]], 3),
                            high_confidence = hc)[hc, ]
  counts_tab[[k]] <- data.frame(comparison = lab[[k]],
                                DEGs_DESeq2 = sum(d$sig), DEGs_voom_donor = sum(v_sig),
                                high_confidence = sum(hc),
                                higher_in_first_group = sum(hc & d$b > 0), lower_in_first_group = sum(hc & d$b < 0),
                                pct_of_DESeq2_DEGs = round(100 * sum(hc) / sum(d$sig), 1))
  write.csv(hc_tab[[k]], file.path(out_dir, paste0("HighConfidence_genes_", k, ".csv")), row.names = FALSE)
}
counts_tab <- do.call(rbind, counts_tab)
show("High-confidence genes (significant in DESeq2 and donor-aware limma-voom, same direction)", counts_tab)
write.csv(counts_tab, file.path(out_dir, "HighConfidence_counts.csv"), row.names = FALSE)

# ---- 3. Interaction without the HD-only batch 4 -----------------------------------
EP <- des$E_P
s1_file <- file.path(sens, "objects", "results_S1_without_batch4.rds")
s1 <- tryCatch(readRDS(s1_file), error = function(e) NULL)
s1_tab <- NULL
if (is.data.frame(s1) && all(c("gene_id", "log2FC_E_P", "padj_E_P") %in% names(s1))) {
  s1_tab <- data.frame(gene_id = s1$gene_id, b = s1$log2FC_E_P, padj = s1$padj_E_P,
                       sig = if ("sig_E_P" %in% names(s1)) s1$sig_E_P %in% TRUE else s1$padj_E_P < 0.05)
} else if (is.list(s1) && "E_P" %in% names(s1) && is.data.frame(s1$E_P)) {
  x <- s1$E_P; bcol <- intersect(c("log2FC", "log2FoldChange", "b"), names(x))[1]
  idc <- if ("gene_id" %in% names(x)) x$gene_id else rownames(x)
  s1_tab <- data.frame(gene_id = idc, b = x[[bcol]], padj = x$padj,
                       sig = if ("sig" %in% names(x)) x$sig %in% TRUE else x$padj < 0.05)
}
if (is.null(s1_tab)) {
  say("\nSaved S1 result not in a known format (structure below); refitting the model without batch 4.")
  if (!is.null(s1)) cat(capture.output(str(s1, max.level = 1)), sep = "\n", file = log_file, append = TRUE)
  hd_only <- names(which(tapply(grepl("HD", ddsP$group), ddsP$batch_study, all)))
  say("HD-only batch removed: ", paste(hd_only, collapse = ", "))
  d1 <- ddsP[, !(ddsP$batch_study %in% hd_only)]
  colData(d1) <- droplevels(colData(d1))
  d1 <- DESeq(d1, quiet = TRUE)
  rn1 <- resultsNames(d1); ref1 <- levels(d1$group)[1]
  g1 <- function(lv) { v <- setNames(numeric(length(rn1)), rn1); if (lv == ref1) return(v)
  v[paste0("group_", lv, "_vs_", ref1)] <- 1; v }
  t1 <- tab_P(d1, g1(SH) - g1(SC) - g1(CH) + g1(CC), c(SH, SC, CH, CC))
  s1_tab <- data.frame(gene_id = t1$gene_id, b = t1$b, padj = t1$padj, sig = t1$sig)
}
m1 <- s1_tab[match(EP$gene_id, s1_tab$gene_id), ]
wo4 <- EP$sig & m1$sig %in% TRUE & sign(EP$b) == sign(m1$b)
say("\nInteraction genes still significant without batch 4, same direction: ", sum(wo4, na.rm = TRUE),
    " of ", sum(EP$sig), " (", round(100 * sum(wo4, na.rm = TRUE) / sum(EP$sig), 1),
    "%; sensitivity run reported 23.9%)")
hcE <- EP$gene_id %in% hc_tab$E_P$gene_id
robust <- hcE & wo4 %in% TRUE
say("Interaction genes high-confidence AND without batch 4: ", sum(robust))
rob_tab <- data.frame(gene_id = EP$gene_id, gene_symbol = unname(sym[EP$gene_id]),
                      interaction_log2FC = round(EP$b, 3), interaction_padj = signif(EP$padj, 3),
                      high_confidence = hcE, significant_without_batch4 = wo4 %in% TRUE)[EP$sig, ]
write.csv(rob_tab, file.path(out_dir, "Interaction_genes_robustness_flags.csv"), row.names = FALSE)

# ---- 4. Core genes (all three regions) --------------------------------------------
core <- function(sA, sB) {
  ids <- Reduce(intersect, list(des$A$gene_id[sA], des$B_P$gene_id[sB], BT$gene_id[BT$sig]))
  a <- des$A[match(ids, des$A$gene_id), ]; b <- des$B_P[match(ids, des$B_P$gene_id), ]
  t <- BT[match(ids, BT$gene_id), ]
  ids[sign(a$b) == sign(b$b) & sign(a$b) == sign(t$b)]
}
core_main <- core(des$A$sig, des$B_P$sig)
core_hc   <- core(des$A$gene_id %in% hc_tab$A$gene_id, des$B_P$gene_id %in% hc_tab$B_P$gene_id)
say("\nCore genes, main lists: ", length(core_main), " (thesis: 40)")
say("Core genes, high-confidence striatal and cingulate lists: ", length(core_hc),
    " (", sum(core_hc %in% core_main), " of them among the original core genes)")
a <- des$A[match(core_hc, des$A$gene_id), ]
write.csv(data.frame(gene_id = core_hc, gene_symbol = unname(sym[core_hc]), direction = ifelse(a$b > 0, "higher in HD", "lower in HD"),
                     in_original_core = core_hc %in% core_main),
          file.path(out_dir, "Core_genes_high_confidence.csv"), row.names = FALSE)

# ---- 5. Cross-study regional check on the stricter interaction sets ---------------
A <- des$A; BP <- des$B_P
ids <- intersect(A$gene_id[A$tested], BT$gene_id[BT$tested])
a <- A[match(ids, A$gene_id), ]; t <- BT[match(ids, BT$gene_id), ]; cc <- BP[match(ids, BP$gene_id), ]
e <- EP[match(ids, EP$gene_id), ]
d_ST <- a$b - t$b; padj_z <- p.adjust(2 * pnorm(-abs(d_ST / sqrt(a$se^2 + t$se^2))), "BH")
sig_ST <- padj_z < 0.05; closer <- abs(t$b - cc$b) < abs(t$b - a$b)
check_set <- function(name, inset) {
  inset <- inset & !is.na(inset)
  ft <- fisher.test(table(factor(inset, c(FALSE, TRUE)), factor(sig_ST, c(FALSE, TRUE))))
  data.frame(gene_set = name, genes_tested_in_both = sum(inset),
             also_differ_from_motor_cortex = sum(inset & sig_ST),
             pct = round(100 * sum(inset & sig_ST) / sum(inset), 1),
             same_direction_as_interaction = sum(inset & sig_ST & sign(e$b) == sign(d_ST)),
             motor_cortex_closer_to_cingulate_pct = round(100 * sum(inset & closer) / sum(inset), 1),
             other_genes_pct = round(100 * sum(!inset & sig_ST) / sum(!inset), 1),
             odds_ratio = round(unname(ft$estimate), 1), fisher_p = signif(ft$p.value, 3))
}
reg <- rbind(check_set("All interaction genes", e$sig),
             check_set("High-confidence interaction genes", ids %in% hc_tab$E_P$gene_id),
             check_set("High-confidence and significant without batch 4", ids %in% rob_tab$gene_id[rob_tab$high_confidence & rob_tab$significant_without_batch4]))
show("Cross-study regional check (striatal HD effect vs motor-cortex HD effect)", reg)
write.csv(reg, file.path(out_dir, "Regional_crossstudy_by_gene_set.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude the console output (or run_log.txt). Outputs: ", out_dir)