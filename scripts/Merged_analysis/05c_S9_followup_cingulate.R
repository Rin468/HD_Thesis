# =============================================================================
#  05c_S9_followup_cingulate.R
#  FOLLOW-UP to S9 of 05_sensitivity_robustness_S1-S9.R (run sensitivity_20261017_163558)
#
#  In S9 (without the 4 sex-corrected samples) the cingulate HD-vs-CTRL list (B_P)
#  dropped from 1,000 to 287 DEGs, while effect sizes were preserved (rho 0.975).
#  3 of the 4 removed samples are cortex. This script separates the candidates:
#    S9a  without H5608_Accum               (relabelled to cortex; 2nd cortex sample of donor H5608)
#    S9b  without C98_Cing                  (sex-corrected cortex control; Cortex_CTRL 10 -> 9)
#    S9c  without H5608_Accum + H5608_Cing  (both cortex samples of donor H5608)
#  Same model, filter and DEG rule as the main run: ~ batch_study + sex + group,
#  >= 10 counts in >= min(10, smallest group) samples, padj < 0.05 on reliable genes.
#
#  Reading:
#    S9a or S9c reproduce the B_P drop -> the main B_P count is inflated by donor H5608
#                                         contributing two near-identical cortex samples
#                                         (pseudo-replication)
#    S9b reproduces the drop            -> losing one of only 10 cortex controls drives it
#    none reproduces it                 -> the drop needs the combination of samples
#
#  Input   newest final_runs/bystudy_* run (objects) + Merged_analysis/results/objects/merged_txi_*
#  Output  NEW folder Merged_analysis/final_runs/S9followup_<date>_<time>/ (run_log.txt, tables/)
#  Run     Session > Restart R, open this file, Source. Runtime ~3-5 min.
#  Save as HD_Thesis/scripts/Merged_analysis/05c_S9_followup_cingulate.R
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))

# ---- 0. Paths and helpers ----------------------------------------------------
base_dir   <- "C:/!Virag.ORI/Bioinfo/szakdoga"
merged_dir <- file.path(base_dir, "Merged_analysis")
in_dir     <- file.path(merged_dir, "results", "objects")
runs_avail <- sort(list.files(file.path(merged_dir, "final_runs"), pattern = "^bystudy_[0-9]{8}_[0-9]{6}$"))
if (length(runs_avail) == 0) stop("No bystudy_* main run found - run 04_main_DESeq2_by_study.R first")
main_run   <- file.path(merged_dir, "final_runs", tail(runs_avail, 1))   # newest main run
run_dir    <- file.path(merged_dir, "final_runs", format(Sys.time(), "S9followup_%Y%m%d_%H%M%S"))
for (s in c("tables", "objects"))
  dir.create(file.path(run_dir, s), recursive = TRUE, showWarnings = FALSE)
out_path <- function(sub, file) file.path(run_dir, sub, file)
log_file <- file.path(run_dir, "run_log.txt")

say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = ""); cat(txt, "\n", sep = "", file = log_file, append = TRUE) }
show_tab <- function(title, x) {
  out <- c("", paste0("-- ", title, " --"), capture.output(print(x)))
  cat(out, sep = "\n"); cat(out, sep = "\n", file = log_file, append = TRUE)
}
latest <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) stop("No file matching ", pattern, " in ", dir)
  f[which.max(file.info(f)$mtime)]
}
scor <- function(x, y) { ok <- complete.cases(x, y); if (sum(ok) < 3) NA_real_ else round(cor(x[ok], y[ok], method = "spearman"), 3) }
isT  <- function(x) x %in% TRUE

say("Run folder: ", run_dir)
say("Main run compared against: ", main_run)
say("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | 05c_S9_followup_cingulate.R")


# ---- 1. Load inputs ----------------------------------------------------------
say("\n==== 1. Load inputs ====  [", format(Sys.time(), "%H:%M:%S"), "]")
txi_file <- latest(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$")
say("Counts: ", txi_file)
txi      <- readRDS(txi_file)
cd_all   <- as.data.frame(readRDS(file.path(main_run, "objects", "colData_final.rds")))
res_main <- readRDS(file.path(main_run, "objects", "results_by_contrast.rds"))
stopifnot(identical(colnames(txi$counts), rownames(cd_all)), nrow(cd_all) == 90,
          identical(txi$countsFromAbundance, "no"),
          all(c("A", "B_P", "C", "D", "E_P") %in% names(res_main)))
cd_all$condition   <- factor(as.character(cd_all$condition), levels = c("CTRL", "HD"))
cd_all$tissue      <- factor(as.character(cd_all$tissue), levels = c("Striatum", "Cortex"))
cd_all$group       <- factor(as.character(cd_all$group),
                             levels = c("Striatum_CTRL", "Striatum_HD", "Cortex_CTRL", "Cortex_HD"))
cd_all$sex         <- factor(as.character(cd_all$sex))
cd_all$batch_study <- factor(as.character(cd_all$batch_study))
if (!"in_model" %in% names(cd_all)) cd_all$in_model <- TRUE
isP <- cd_all$study == "Paryani" & cd_all$in_model
say("Paryani samples in the main model: ", sum(isP),
    " | excluded as in the main run: ", paste(cd_all$sample_label[!cd_all$in_model], collapse = ", "))

candidates <- c("H5608_Accum", "H5608_Caud", "H5608_Cing", "C98_Cing")
show_tab("The 4 sex-corrected samples as used in the main model",
         cd_all[cd_all$sample_label %in% candidates,
                intersect(c("sample_label", "donor", "group", "sex", "batch_study"), names(cd_all))])


# ---- 2. Helpers copied from 05_sensitivity_robustness_S1-S9.R ------------------
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

g_SHD <- "group_Striatum_HD_vs_Striatum_CTRL"
g_CC  <- "group_Cortex_CTRL_vs_Striatum_CTRL"
g_CHD <- "group_Cortex_HD_vs_Striatum_CTRL"
full_specs <- list(
  A   = list(w = setNames(1, g_SHD), sets = c("group == 'Striatum_HD'", "group == 'Striatum_CTRL'")),
  B_P = list(w = setNames(c(1, -1), c(g_CHD, g_CC)), sets = c("group == 'Cortex_HD'", "group == 'Cortex_CTRL'")),
  C   = list(w = setNames(-1, g_CC), sets = c("group == 'Striatum_CTRL'", "group == 'Cortex_CTRL'")),
  D   = list(w = setNames(c(1, -1), c(g_SHD, g_CHD)), sets = c("group == 'Striatum_HD'", "group == 'Cortex_HD'")),
  E_P = list(w = setNames(c(1, -1, 1), c(g_SHD, g_CHD, g_CC)),
             sets = c("group == 'Striatum_HD'", "group == 'Striatum_CTRL'", "group == 'Cortex_HD'", "group == 'Cortex_CTRL'")))

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


# ---- 3. Follow-up scenarios ----------------------------------------------------
scenarios <- list(
  S9a_without_H5608_Accum      = c("H5608_Accum"),
  S9b_without_C98_Cing         = c("C98_Cing"),
  S9c_without_H5608_both_cortex = c("H5608_Accum", "H5608_Cing"))

say("\n==== 2. Follow-up DESeq2 reruns ====  [", format(Sys.time(), "%H:%M:%S"), "]")
comp_list <- list(); alt_results <- list()
for (nm in names(scenarios)) {
  drop <- scenarios[[nm]]
  sel  <- isP & !(cd_all$sample_label %in% drop)
  excl <- cd_all$sample_label[isP & !sel]
  say("\n--- ", nm, " ---")
  say("  excluded: ", if (length(excl)) paste(excl, collapse = ", ") else "none")
  if (!setequal(excl, drop)) {
    say("  expected to exclude ", paste(drop, collapse = ", "), " - check sample_label spelling; scenario skipped")
    next
  }
  res <- tryCatch(run_scenario(nm, sel, ~ batch_study + sex + group, 10, "group", full_specs),
                  error = function(e) { say("  SCENARIO FAILED: ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  alt_results[[nm]] <- res
  comp_list[[nm]]   <- compare_to_main(res, nm)
  invisible(gc())
}


# ---- 4. Summary ----------------------------------------------------------------
say("\n==== 3. Summary ====  [", format(Sys.time(), "%H:%M:%S"), "]")
if (length(comp_list) == 0) stop("No scenario finished - see the log above")
comp <- do.call(rbind, comp_list); rownames(comp) <- NULL
show_tab("ALL FOLLOW-UP ANALYSES vs MAIN DESeq2 (padj < 0.05, reliable genes)", comp)
show_tab("Cingulate HD vs CTRL (B_P) only", comp[comp$contrast == "B_P", ])
say("\nFor reference (sensitivity_20261017_163558): main B_P 1,000 DEGs; ",
    "S9 without all 4 sex-corrected samples: 287 DEGs, retained 28.7%, rho 0.975.")
write.csv(comp, out_path("tables", "S9followup_comparison_vs_main.csv"), row.names = FALSE)
saveRDS(alt_results, out_path("objects", "S9followup_results.rds"))
writeLines(capture.output(sessionInfo()), out_path("objects", "sessionInfo.txt"))
say("\nFINISHED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("All outputs: ", run_dir)
