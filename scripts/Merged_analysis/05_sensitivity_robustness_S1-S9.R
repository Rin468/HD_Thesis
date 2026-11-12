# =============================================================================
#  05_sensitivity_Paryani.R
#  SENSITIVITY ANALYSES of the main per-study DESeq2 results
#  (main run: final_runs/bystudy_20261008_230208, script 04_main_DESeq2_by_study.R)
#
#  PART 1  Donor-aware check (Paryani, 76 samples, same design as Model P):
#            limma-voom on lengthScaledTPM counts, design ~ batch_study + sex + group
#              (a) voom_noblock : no donor term        -> shows the effect of the METHOD alone
#              (b) voom_donor   : duplicateCorrelation, block = donor
#                                                       -> adds the repeated-donor correction
#            Reading: DESeq2 vs (a) = method difference; (a) vs (b) = donor effect.
#  PART 2  DESeq2 reruns (same model, subset of samples / other model structure):
#            S1 without batch 4 (HD-only batch)
#            S2 without the 10 tissue-corrected / flagged samples
#            S3 without C5596_Accum (possible outlier)
#            S4 striatum-only model   (~ batch_study + sex + condition)  -> contrast A
#            S5 cortex-only model     (~ batch_study + sex + condition)  -> contrast B_P
#            S6 without juvenile-onset samples
#            S7 GSE79666 without CTRL-3 (~ sex + condition)               -> contrast B_T
#            S8 % mitochondrial reads as an extra covariate
#               (~ batch_study + sex + mt_pct_z + group; batch 4 has a much higher share)
#            S9 without the 4 samples whose sex was corrected from expression
#               (H5608_Accum, H5608_Caud, H5608_Cing, C98_Cing)
#  CHANGES 2026-10-15: compares against the NEWEST bystudy_* run (printed in the log);
#  samples with in_model = FALSE in that run (OF045) are left out everywhere; S8 added.
#  CHANGES 2026-10-17: S9 added (the run stops S9 unless exactly 4 samples are removed);
#  same_dir_pct_of_tested fixed - the numerator now counts only genes tested in the
#  new analysis (before, genes with padj = NA but a log2FC could push it above 100%).
#  retained_pct and rho_log2FC were not affected by the bug.
#  Every analysis is compared with the main result: DEGs retained (significant in
#  both, same direction), same-direction %, Spearman rho of log2FC, DEGs only in the
#  new analysis, and the 3-region "core" genes (HD_common_all_three).
#
#  NOTE: subsets have fewer samples, so fewer DEGs are expected. Judge robustness by
#  the same-direction % and rho first, by retained DEGs second.
#
#  Input   main run objects (colData_final, results_by_contrast, MASTER table) and
#          Merged_analysis/results/objects/merged_txi_*, gene_annotation_*
#  Output  NEW folder Merged_analysis/final_runs/sensitivity_<date>_<time>/
#  Run     Session > Restart R, open this file, Source. Runtime ~15-30 min.
#  Save as HD_Thesis/scripts/Merged_analysis/05_sensitivity_Paryani.R
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
for (p in c("limma", "edgeR")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is missing: BiocManager::install('", p, "')")
}
suppressPackageStartupMessages(library(edgeR))   # also attaches limma

# ---- 0. Paths and helpers ----------------------------------------------------
base_dir   <- "C:/!Virag.ORI/Bioinfo/szakdoga"
merged_dir <- file.path(base_dir, "Merged_analysis")
in_dir     <- file.path(merged_dir, "results", "objects")
runs_avail <- sort(list.files(file.path(merged_dir, "final_runs"), pattern = "^bystudy_[0-9]{8}_[0-9]{6}$"))
if (length(runs_avail) == 0) stop("No bystudy_* main run found - run 04_main_DESeq2_by_study.R first")
main_run   <- file.path(merged_dir, "final_runs", tail(runs_avail, 1))   # newest main run
run_dir    <- file.path(merged_dir, "final_runs", format(Sys.time(), "sensitivity_%Y%m%d_%H%M%S"))
for (s in c("tables", "objects", "figures"))
  dir.create(file.path(run_dir, s), recursive = TRUE, showWarnings = FALSE)
out_path <- function(sub, file) file.path(run_dir, sub, file)
log_file <- file.path(run_dir, "run_log.txt")

say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = ""); cat(txt, "\n", sep = "", file = log_file, append = TRUE) }
show_tab <- function(title, x) {
  out <- c("", paste0("-- ", title, " --"), capture.output(print(x)))
  cat(out, sep = "\n"); cat(out, sep = "\n", file = log_file, append = TRUE)
}
run_step <- function(name, expr) {
  say("\n==== ", name, " ====  [", format(Sys.time(), "%H:%M:%S"), "]")
  tryCatch(expr, error = function(e) { say("ERROR in '", name, "': ", conditionMessage(e)); stop(e) })
  invisible(NULL)
}
latest <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) stop("No file matching ", pattern, " in ", dir)
  f[which.max(file.info(f)$mtime)]
}
save_png <- function(plot, file, w = 8.5, h = 6.5)
  tryCatch(ggsave(out_path("figures", file), plot, width = w, height = h, dpi = 300),
           error = function(e) say("Figure ", file, " could not be saved: ", conditionMessage(e)))
scor <- function(x, y) { ok <- complete.cases(x, y); if (sum(ok) < 3) NA_real_ else round(cor(x[ok], y[ok], method = "spearman"), 3) }
isT  <- function(x) x %in% TRUE
theme_thesis <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        panel.border = element_rect(colour = "#c3c2b7"),
        axis.text = element_text(colour = "#52514e"), plot.subtitle = element_text(colour = "#52514e"))

say("Run folder: ", run_dir)
say("Main run compared against: ", main_run)
say("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | 05_sensitivity_Paryani.R")


# ---- 1. Load inputs ----------------------------------------------------------
run_step("1. Load inputs", {
  txi      <- readRDS(latest(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$"))
  ann      <- as.data.frame(readRDS(latest(in_dir, "^gene_annotation_GENCODEv44(_v[0-9]+)?\\.rds$")))
  cd_all   <- as.data.frame(readRDS(file.path(main_run, "objects", "colData_final.rds")))
  res_main <- readRDS(file.path(main_run, "objects", "results_by_contrast.rds"))
  mt       <- readRDS(file.path(main_run, "objects", "MASTER_DESeq2_by_study.rds"))
  stopifnot(identical(colnames(txi$counts), rownames(cd_all)), nrow(cd_all) == 90,
            identical(txi$countsFromAbundance, "no"),
            all(c("A", "B_P", "B_T", "C", "D", "E_P") %in% names(res_main)))
  cd_all$condition   <- factor(as.character(cd_all$condition), levels = c("CTRL", "HD"))
  cd_all$tissue      <- factor(as.character(cd_all$tissue), levels = c("Striatum", "Cortex"))
  cd_all$group       <- factor(as.character(cd_all$group),
                               levels = c("Striatum_CTRL", "Striatum_HD", "Cortex_CTRL", "Cortex_HD"))
  cd_all$sex         <- factor(as.character(cd_all$sex))
  cd_all$batch_study <- factor(as.character(cd_all$batch_study))
  cd_all$donor       <- as.character(cd_all$donor)
  stopifnot(nrow(ann) == nrow(txi$counts) || all(rownames(txi$counts) %in% ann$gene_id))
  if (!"in_model" %in% names(cd_all)) cd_all$in_model <- TRUE
  isP <- cd_all$study == "Paryani" & cd_all$in_model
  isL <- cd_all$study == "Tseng" & cd_all$in_model
  say("Excluded as in the main run: ",
      if (any(!cd_all$in_model)) paste(cd_all$sample_label[!cd_all$in_model], collapse = ", ") else "none")
  # % of counts from the 37 mitochondrial (MT-) genes; z-scored within Paryani for S8
  mt_rows <- ann$gene_id[grepl("^MT-", ann$gene_symbol)]
  cd_all$mt_pct <- round(100 * colSums(txi$counts[intersect(mt_rows, rownames(txi$counts)), ]) / colSums(txi$counts), 2)
  cd_all$mt_pct_z <- NA_real_
  cd_all$mt_pct_z[isP] <- as.vector(scale(cd_all$mt_pct[isP]))
  say("Mitochondrial genes found: ", length(mt_rows))
  show_tab("% mitochondrial counts by batch (Paryani, median [min-max])",
           tapply(cd_all$mt_pct[isP], droplevels(cd_all$batch_study[isP]),
                  function(x) sprintf("%.1f [%.1f-%.1f]", median(x), min(x), max(x))))
  write.csv(cd_all[, c("sample_id", "sample_label", "study", "batch_study", "group", "mt_pct")],
            out_path("tables", "QC_mitochondrial_percent_per_sample.csv"), row.names = FALSE)
  core_ids <- mt$gene_id[isT(mt$HD_common_all_three)]
  say("Samples: ", sum(isP), " Paryani + ", sum(isL), " GSE79666 | core genes (significant in all three regions, same direction): ",
      length(core_ids))
  show_tab("Paryani: corrected group x batch", table(cd_all$group[isP], droplevels(cd_all$batch_study[isP])))
  show_tab("Paryani: donors per sample (number of donors with n samples)", table(table(cd_all$donor[isP])))
  assign("txi", txi, envir = .GlobalEnv); assign("ann", ann, envir = .GlobalEnv)
  assign("cd_all", cd_all, envir = .GlobalEnv); assign("res_main", res_main, envir = .GlobalEnv)
  assign("core_ids", core_ids, envir = .GlobalEnv)
  assign("isP", isP, envir = .GlobalEnv); assign("isL", isL, envir = .GlobalEnv)
})


# ---- Comparison helpers ------------------------------------------------------
# alt = named list (names = main contrast names) of data.frames with log2FC, padj, reliable, sig
compare_to_main <- function(alt, analysis) {
  do.call(rbind, lapply(names(alt), function(k) {
    r <- res_main[[k]]
    a <- alt[[k]][match(rownames(r), rownames(alt[[k]])), , drop = FALSE]
    r_sig <- isT(r$sig); a_sig <- isT(a$sig)
    same  <- isT(sign(a$log2FC) == sign(r$log2FC))
    tested <- !is.na(a$padj)
    ok <- isT(r$reliable) & isT(a$reliable) & is.finite(r$log2FC) & is.finite(a$log2FC)
    data.frame(analysis = analysis, contrast = k,
               main_DEGs = sum(r_sig), alt_DEGs = sum(a_sig),
               main_DEGs_tested_in_alt = sum(r_sig & tested),
               retained_sig_same_dir = sum(r_sig & a_sig & same),
               retained_pct = round(100 * sum(r_sig & a_sig & same) / max(1, sum(r_sig)), 1),
               same_dir_pct_of_tested = round(100 * sum(r_sig & tested & same) / max(1, sum(r_sig & tested)), 1),
               alt_only_DEGs = sum(a_sig & !r_sig),
               rho_log2FC = scor(r$log2FC[ok], a$log2FC[ok]),
               stringsAsFactors = FALSE)
  }))
}
core_check <- function(alt, analysis) {
  do.call(rbind, lapply(intersect(names(alt), c("A", "B_P")), function(k) {
    r <- res_main[[k]]; a <- alt[[k]]
    ia <- match(core_ids, rownames(a)); ir <- match(core_ids, rownames(r))
    same <- isT(sign(a$log2FC[ia]) == sign(r$log2FC[ir]))
    data.frame(analysis = analysis, contrast = k, core_genes = length(core_ids),
               tested = sum(!is.na(a$padj[ia])), same_direction = sum(same),
               significant_same_direction = sum(same & isT(a$sig[ia])), stringsAsFactors = FALSE)
  }))
}
comp_list <- list(); core_list <- list(); alt_results <- list()


# ---- 2. PART 1: limma-voom, with and without donor blocking -------------------
run_step("2. Part 1: limma-voom (no block vs donor block), Paryani", {
  txiP <- list(abundance = txi$abundance[, isP], counts = txi$counts[, isP], length = txi$length[, isP])
  cdP  <- droplevels(cd_all[isP, ])
  # lengthScaledTPM counts exactly as tximport(countsFromAbundance = "lengthScaledTPM") computes them
  nc  <- txiP$abundance * rowMeans(txiP$length)
  lst <- t(t(nc) * (colSums(txiP$counts) / colSums(nc)))
  genes <- rownames(res_main$A)                       # same genes as Model P
  stopifnot(all(genes %in% rownames(lst)))
  say("lengthScaledTPM vs original counts, Spearman rho over all genes x samples (sanity): ",
      scor(as.vector(log1p(lst[genes, ])), as.vector(log1p(txiP$counts[genes, ]))))
  y <- calcNormFactors(DGEList(counts = lst[genes, ]))   # TMM
  dmat <- model.matrix(~ batch_study + sex + group, data = cdP)
  if (qr(dmat)$rank < ncol(dmat)) stop("limma design not full rank")
  say("limma design columns: ", paste(colnames(dmat), collapse = ", "))
  need <- c("groupStriatum_HD", "groupCortex_CTRL", "groupCortex_HD")
  if (!all(need %in% colnames(dmat))) stop("Unexpected design column names")
  cm <- matrix(0, ncol(dmat), 5, dimnames = list(colnames(dmat), c("A", "B_P", "C", "D", "E_P")))
  cm["groupStriatum_HD", "A"] <- 1
  cm[c("groupCortex_HD", "groupCortex_CTRL"), "B_P"] <- c(1, -1)
  cm["groupCortex_CTRL", "C"] <- -1
  cm[c("groupStriatum_HD", "groupCortex_HD"), "D"] <- c(1, -1)
  cm[c("groupStriatum_HD", "groupCortex_HD", "groupCortex_CTRL"), "E_P"] <- c(1, -1, 1)

  extract <- function(fit2) {
    setNames(lapply(colnames(cm), function(k) {
      p <- fit2$p.value[, k]
      data.frame(log2FC = fit2$coefficients[, k],
                 lfcSE  = sqrt(fit2$s2.post) * fit2$stdev.unscaled[, k],
                 pvalue = p, padj = p.adjust(p, "BH"),
                 reliable = res_main[[k]]$reliable[match(genes, rownames(res_main[[k]]))],
                 row.names = genes)
    }), colnames(cm))
  }
  add_sig <- function(l) lapply(l, function(d) { d$sig <- isT(d$reliable) & isT(d$padj < 0.05); d })

  # (a) no donor term
  v1   <- voom(y, dmat)
  fitA <- eBayes(contrasts.fit(lmFit(v1, dmat), cm))
  alt_results$voom_noblock <- add_sig(extract(fitA))

  # (b) donor as a random-like block via duplicateCorrelation (two voom/dupCor rounds)
  blk <- cdP$donor
  dc1 <- duplicateCorrelation(v1, dmat, block = blk)
  v2  <- voom(y, dmat, block = blk, correlation = dc1$consensus)
  dc2 <- duplicateCorrelation(v2, dmat, block = blk)
  fitB <- eBayes(contrasts.fit(lmFit(v2, dmat, block = blk, correlation = dc2$consensus), cm))
  alt_results$voom_donor <- add_sig(extract(fitB))
  dc_tab <- data.frame(round = c(1, 2), consensus_correlation = round(c(dc1$consensus.correlation, dc2$consensus.correlation), 4),
                       donors = length(unique(blk)), samples = length(blk))
  show_tab("duplicateCorrelation: consensus within-donor correlation", dc_tab)
  write.csv(dc_tab, out_path("tables", "SENS_limma_duplicateCorrelation.csv"), row.names = FALSE)
  if (dc2$consensus.correlation < 0.05)
    say("NOTE: consensus correlation is close to 0 - donor blocking has little effect.")
  if (dc2$consensus.correlation > 0.9)
    say("WARNING: consensus correlation close to 1 - check donor IDs (possible confounding of donor with group).")

  for (nm in c("voom_noblock", "voom_donor")) {
    comp_list[[nm]] <- compare_to_main(alt_results[[nm]], nm)
    core_list[[nm]] <- core_check(alt_results[[nm]], nm)
    wide <- do.call(cbind, lapply(names(alt_results[[nm]]), function(k) {
      d <- alt_results[[nm]][[k]][, c("log2FC", "lfcSE", "padj", "sig")]
      setNames(d, paste0(names(d), "_", k))
    }))
    write.csv(data.frame(gene_id = genes, gene_symbol = ann$gene_symbol[match(genes, ann$gene_id)], wide),
              out_path("tables", paste0("SENS_", nm, "_results.csv")), row.names = FALSE)
  }
  show_tab("limma-voom vs main DESeq2 (Model P)", rbind(comp_list$voom_noblock, comp_list$voom_donor))

  # Scatter: DESeq2 main vs voom_donor log2FC for A, B_P, E_P
  pd <- do.call(rbind, lapply(c("A", "B_P", "E_P"), function(k) {
    r <- res_main[[k]]; a <- alt_results$voom_donor[[k]]
    ia <- match(rownames(r), rownames(a))
    ok <- isT(r$reliable) & is.finite(r$log2FC) & is.finite(a$log2FC[ia])
    st <- ifelse(isT(r$sig) & isT(a$sig[ia]), "significant in both",
                 ifelse(isT(r$sig), "DESeq2 main only", ifelse(isT(a$sig[ia]), "voom (donor) only", "neither")))
    data.frame(contrast = k, deseq2 = r$log2FC[ok], voom = a$log2FC[ia][ok], status = st[ok])
  }))
  pd$status <- factor(pd$status, levels = c("neither", "DESeq2 main only", "voom (donor) only", "significant in both"))
  pd <- pd[order(pd$status), ]
  p_sc <- ggplot(pd, aes(deseq2, voom, colour = status)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#898781") +
    geom_point(size = 0.6, alpha = 0.6) + facet_wrap(~ contrast, scales = "free") +
    scale_colour_manual(values = c("neither" = "#d9d8d2", "DESeq2 main only" = "#eb6834",
                                   "voom (donor) only" = "#2a78d6", "significant in both" = "#1baf7a")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "log2FC, DESeq2 main (Model P)", y = "log2FC, limma-voom with donor blocking",
         colour = "padj < 0.05", title = "Main DESeq2 vs donor-aware limma-voom (Paryani)") +
    theme_thesis
  save_png(p_sc, "SENS_DESeq2_vs_voom_donor_log2FC.png", w = 11, h = 4.5)
  rm(txiP, nc, lst, y, v1, v2, fitA, fitB); invisible(gc())
})


# ---- 3. PART 2: DESeq2 reruns -------------------------------------------------
g_SHD <- "group_Striatum_HD_vs_Striatum_CTRL"
g_CC  <- "group_Cortex_CTRL_vs_Striatum_CTRL"
g_CHD <- "group_Cortex_HD_vs_Striatum_CTRL"
c_HD  <- "condition_HD_vs_CTRL"
full_specs <- list(
  A   = list(w = setNames(1, g_SHD), sets = c("group == 'Striatum_HD'", "group == 'Striatum_CTRL'")),
  B_P = list(w = setNames(c(1, -1), c(g_CHD, g_CC)), sets = c("group == 'Cortex_HD'", "group == 'Cortex_CTRL'")),
  C   = list(w = setNames(-1, g_CC), sets = c("group == 'Striatum_CTRL'", "group == 'Cortex_CTRL'")),
  D   = list(w = setNames(c(1, -1), c(g_SHD, g_CHD)), sets = c("group == 'Striatum_HD'", "group == 'Cortex_HD'")),
  E_P = list(w = setNames(c(1, -1, 1), c(g_SHD, g_CHD, g_CC)),
             sets = c("group == 'Striatum_HD'", "group == 'Striatum_CTRL'", "group == 'Cortex_HD'", "group == 'Cortex_CTRL'")))
cond_spec <- function(name) setNames(list(list(w = setNames(1, c_HD), sets = c("condition == 'HD'", "condition == 'CTRL'"))), name)

run_scenario <- function(tag, sel, design, cap, group_var, specs) {
  cdx <- droplevels(cd_all[sel, ])
  mm  <- model.matrix(design, data = cdx)
  if (qr(mm)$rank < ncol(mm)) stop("design not full rank (", ncol(mm), " columns)")
  min_s <- min(cap, min(table(cdx[[group_var]])))
  sub <- list(abundance = txi$abundance[, sel], counts = txi$counts[, sel], length = txi$length[, sel],
              countsFromAbundance = "no")
  dds <- DESeqDataSetFromTximport(sub, colData = cdx, design = design)
  keep <- rowSums(counts(dds) >= 10) >= min_s
  dds <- dds[keep, ]
  say(tag, ": ", ncol(dds), " samples | ", paste(deparse(design), collapse = ""), " | genes kept: ", sum(keep),
      " (>= 10 counts in >= ", min_s, " samples)")
  say("  samples per ", group_var, ": ", paste(names(table(cdx[[group_var]])), as.vector(table(cdx[[group_var]])), sep = "=", collapse = ", "))
  dds <- tryCatch(DESeq(dds), error = function(e) {
    say("  DESeq() failed (", conditionMessage(e), ") - retrying with fitType = 'local'")
    DESeq(dds, fitType = "local")
  })
  bad <- !(mcols(dds)$betaConv %in% TRUE)
  say("  non-converged genes removed: ", sum(bad))
  dds <- dds[!bad, ]
  rn <- resultsNames(dds)
  cnt <- counts(dds)
  out <- list()
  for (k in names(specs)) {
    s <- specs[[k]]
    if (!all(names(s$w) %in% rn)) stop("coefficient missing for ", k, ": ", paste(setdiff(names(s$w), rn), collapse = ", "))
    cv <- setNames(numeric(length(rn)), rn); cv[names(s$w)] <- s$w
    res <- results(dds, contrast = cv, alpha = 0.05)
    zero <- Reduce(`|`, lapply(s$sets, function(e) rowSums(cnt[, eval(parse(text = e), cdx), drop = FALSE]) == 0))
    rel  <- !zero & !(!is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 10)
    out[[k]] <- data.frame(log2FC = res$log2FoldChange, lfcSE = res$lfcSE, pvalue = res$pvalue, padj = res$padj,
                           reliable = rel, sig = rel & isT(res$padj < 0.05), row.names = rownames(res))
    say(sprintf("  %-4s significant: %5d (up %5d, down %5d)", k, sum(out[[k]]$sig),
                sum(out[[k]]$sig & isT(out[[k]]$log2FC > 0)), sum(out[[k]]$sig & isT(out[[k]]$log2FC < 0))))
  }
  attr(out, "n_samples") <- ncol(dds); attr(out, "n_genes") <- nrow(dds)
  out
}

scenarios <- list(
  S1_without_batch4 = list(sel = quote(isP & cd_all$batch_study != "Paryani_b4"),
                           design = ~ batch_study + sex + group, cap = 10, gv = "group", specs = full_specs),
  S2_without_corrected_flagged = list(sel = quote(isP & cd_all$region_check == "consistent"),
                                      design = ~ batch_study + sex + group, cap = 10, gv = "group", specs = full_specs),
  S3_without_C5596_Accum = list(sel = quote(isP & cd_all$sample_label != "C5596_Accum"),
                                design = ~ batch_study + sex + group, cap = 10, gv = "group", specs = full_specs),
  S4_striatum_only = list(sel = quote(isP & cd_all$tissue == "Striatum"),
                          design = ~ batch_study + sex + condition, cap = 10, gv = "condition", specs = cond_spec("A")),
  S5_cortex_only = list(sel = quote(isP & cd_all$tissue == "Cortex"),
                        design = ~ batch_study + sex + condition, cap = 10, gv = "condition", specs = cond_spec("B_P")),
  S6_without_juvenile_onset = list(sel = quote(isP & !(cd_all$juvenile_onset %in% TRUE)),
                                   design = ~ batch_study + sex + group, cap = 10, gv = "group", specs = full_specs),
  S7_GSE79666_without_CTRL3 = list(sel = quote(isL & cd_all$sample_label != "CTRL-3"),
                                   design = ~ sex + condition, cap = 7, gv = "condition", specs = cond_spec("B_T")),
  S8_mitochondrial_covariate = list(sel = quote(isP),
                                    design = ~ batch_study + sex + mt_pct_z + group, cap = 10, gv = "group",
                                    specs = full_specs),
  S9_without_sex_corrected = list(sel = quote(isP & !(grepl("^H5608[_-]", cd_all$sample_label) |
                                                        grepl("^C98[_-]Cing$", cd_all$sample_label))),
                                  design = ~ batch_study + sex + group, cap = 10, gv = "group", specs = full_specs))

run_step("3. Part 2: DESeq2 reruns", {
  scen_info <- list()
  for (nm in names(scenarios)) {
    s <- scenarios[[nm]]
    sel <- eval(s$sel)
    n_excl <- sum(isP | isL) - sum(sel)
    say("\n--- ", nm, " ---")
    if (sum(sel) == sum(isP) || sum(sel) == sum(isL)) {
      if (!grepl("only|covariate", nm)) { say("  nothing was excluded (sample label or flag not found) - scenario skipped"); next }
    }
    if (nm == "S9_without_sex_corrected") {
      excl <- cd_all$sample_label[isP & !sel]
      say("  excluded: ", if (length(excl)) paste(excl, collapse = ", ") else "none")
      if (length(excl) != 4) {
        say("  expected 4 excluded samples, found ", length(excl), " - check sample_label spelling; scenario skipped")
        next
      }
    }
    if (nm == "S7_GSE79666_without_CTRL3") {
      show_tab("GSE79666 sex x condition (all 14)", table(cd_all$sex[isL], cd_all$condition[isL]))
      show_tab("GSE79666 sample labels", cd_all[isL, c("sample_label", "sex", "condition")])
    }
    res <- tryCatch(run_scenario(nm, sel, s$design, s$cap, s$gv, s$specs), error = function(e) {
      say("  SCENARIO FAILED: ", conditionMessage(e)); NULL })
    if (is.null(res)) next
    alt_results[[nm]] <- res
    comp_list[[nm]] <- compare_to_main(res, nm)
    core_list[[nm]] <- core_check(res, nm)
    scen_info[[nm]] <- data.frame(scenario = nm, samples_used = attr(res, "n_samples"),
                                  genes_tested = attr(res, "n_genes"),
                                  design = paste(deparse(s$design), collapse = ""), stringsAsFactors = FALSE)
    saveRDS(res, out_path("objects", paste0("results_", nm, ".rds")))
    invisible(gc())
  }
  if (length(scen_info)) {
    sc <- do.call(rbind, scen_info); rownames(sc) <- NULL
    show_tab("Scenario overview", sc)
    write.csv(sc, out_path("tables", "SENS_scenarios_overview.csv"), row.names = FALSE)
  }
})


# ---- 4. Summary tables and figure --------------------------------------------
run_step("4. Summary tables and figure", {
  comp <- do.call(rbind, comp_list); rownames(comp) <- NULL
  core <- do.call(rbind, core_list); rownames(core) <- NULL
  show_tab("ALL ANALYSES vs MAIN DESeq2 (padj < 0.05, reliable genes)", comp)
  show_tab("3-region core genes: behaviour in each analysis", core)
  write.csv(comp, out_path("tables", "SENS_comparison_vs_main.csv"), row.names = FALSE)
  write.csv(core, out_path("tables", "SENS_core_genes_check.csv"), row.names = FALSE)
  saveRDS(alt_results, out_path("objects", "sens_all_results.rds"))
  p <- ggplot(comp, aes(contrast, retained_pct, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    scale_fill_viridis_d(end = 0.9) +
    labs(x = "Main contrast", y = "Main-run DEGs retained (%), significant in both with the same direction",
         fill = "Analysis", title = "Robustness of the main DEG lists",
         subtitle = "Subsets have fewer samples, so lower retention is expected; see same-direction % and rho in the table") +
    theme_thesis
  save_png(p, "SENS_retention_of_main_DEGs.png", w = 10, h = 6)
  writeLines(capture.output(sessionInfo()), out_path("objects", "sessionInfo.txt"))
  say("\nFINISHED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  say("All outputs: ", run_dir)
  say("Please send: run_log.txt (the two summary tables at the end are the key part) + ",
      "SENS_DESeq2_vs_voom_donor_log2FC.png + SENS_retention_of_main_DEGs.png")
})
