# =============================================================================
#  06_null_calibration_Paryani.R
#  NULL CALIBRATION: how often do the methods find "differences" where none exist?
#
#  Samples  Paryani striatum HD samples only (28 samples, repeated donors),
#           corrected tissue labels from the main run. Inside this set there is
#           NO disease contrast: samples are split at random into two
#           pseudo-groups (g1 / g2), and "g2 vs g1" is tested.
#  Two kinds of random split
#    donor_level   all samples of a donor stay in the same pseudo-group
#                  (this is where a method that ignores donors can go wrong)
#    sample_level  samples are split at random, ignoring donors (control)
#  Three analyses per split, all with the main design ~ batch_study + sex + pseudo
#    DESeq2           (same as the main analysis, tximport length offsets)
#    voom_noblock     limma-voom, donors treated as independent
#    voom_donor       limma-voom + duplicateCorrelation(block = donor)
#  Measured per split: genes with padj < 0.05, fraction of genes with p < 0.05
#  (expected 0.05 under a correct null), and the inflation factor lambda
#  (median chi-square from the p-values / 0.455; expected 1).
#
#  Interpretation guide
#    A calibrated method gives lambda ~ 1, ~5% of p-values < 0.05 and (almost)
#    no genes with padj < 0.05 in a random split.
#    If a method is clearly inflated in donor_level splits but not in
#    sample_level splits, donor dependence is the cause.
#    NOTE: sample_level splits put correlated samples of one donor into both
#    groups; methods that ignore donors are then expected to be a little
#    CONSERVATIVE (lambda < 1) - this is not an error.
#  Limits: 28 samples of one group; heterogeneity inside HD (onset, batch) adds
#    noise that is random with respect to the labels but not zero; this tests
#    the donor structure of striatum HD, not every real contrast.
#
#  HOW TO RUN
#   1) Leave test_mode <- TRUE: 2 splits per type, a few minutes. Check the log.
#   2) Set test_mode <- FALSE for the full run (50 voom / 10 DESeq2 splits per
#      type by default). The log prints the time of every split.
#  Output  NEW folder Merged_analysis/final_runs/nullcal_<date>_<time>/
#  Run     Session > Restart R, open this file, Source.
# =============================================================================

suppressPackageStartupMessages({ library(DESeq2); library(edgeR); library(ggplot2) })  # edgeR attaches limma

# ---- Settings -------------------------------------------------------------------
# 31 Oct 2026: final thesis run. 20 splits per type for all three methods
# (DESeq2 is the method being judged, so it gets as many splits as voom).
test_mode      <- FALSE
n_splits_voom  <- 20      # splits per type for the two limma-voom analyses
n_splits_deseq <- 20      # splits per type for DESeq2 (the same splits as voom)
if (test_mode) { n_splits_voom <- 2; n_splits_deseq <- 2 }
set.seed(20261012)

base_dir   <- "C:/!Virag.ORI/Bioinfo/szakdoga"
merged_dir <- file.path(base_dir, "Merged_analysis")
in_dir     <- file.path(merged_dir, "results", "objects")
# 31 Oct 2026: point to the analysis of record (final sample decisions, OF045 excluded).
# The 2-split run of 12 Oct used the superseded provisional run bystudy_20261008_230208.
main_run   <- file.path(merged_dir, "final_runs", "bystudy_20261015_151051")
if (!file.exists(file.path(main_run, "objects", "colData_final.rds"))) {
  print(list.files(file.path(main_run, "objects")))
  stop("colData_final.rds not found in the final main run - send Claude the file list printed above.")
}
run_dir    <- file.path(merged_dir, "final_runs", format(Sys.time(), "nullcal_%Y%m%d_%H%M%S"))
for (s in c("tables", "objects", "figures")) dir.create(file.path(run_dir, s), recursive = TRUE, showWarnings = FALSE)
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
say("Run folder: ", run_dir)
say("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | test_mode = ", test_mode,
    " | splits per type: voom ", n_splits_voom, ", DESeq2 ", n_splits_deseq)

# ---- 1. Samples and data ----------------------------------------------------------
txi <- readRDS(latest(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$"))
cd  <- as.data.frame(readRDS(file.path(main_run, "objects", "colData_final.rds")))
stopifnot(identical(txi$countsFromAbundance, "no"), all(rownames(cd) %in% colnames(txi$counts)))
# The final run excludes OF045, so keep only the samples of the final colData, in its order
for (k in c("abundance", "counts", "length")) txi[[k]] <- txi[[k]][, rownames(cd), drop = FALSE]
stopifnot(identical(colnames(txi$counts), rownames(cd)))
sel <- cd$study == "Paryani" & as.character(cd$group) == "Striatum_HD"
cdx <- droplevels(cd[sel, ])
cdx$batch_study <- factor(as.character(cdx$batch_study))
cdx$sex         <- factor(as.character(cdx$sex))
cdx$donor       <- as.character(cdx$donor)
say("Striatum HD samples: ", nrow(cdx), " | donors: ", length(unique(cdx$donor)),
    " | donors with >= 2 samples: ", sum(table(cdx$donor) > 1))
if (!"donor" %in% names(cd)) stop("No 'donor' column in colData_final.rds - send Claude names(cd).")
if (nrow(cdx) != 28) say("NOTE: ", nrow(cdx), " striatum HD samples (28 expected in the 76-sample design; ",
                         "27 is expected if OF045 is excluded in the final run).")
show_tab("Samples per donor (number of donors with n samples)", table(table(cdx$donor)))
show_tab("Batch x sex", table(cdx$batch_study, cdx$sex))

# Gene filter independent of the random labels
keep <- rowSums(txi$counts[, sel] >= 10) >= 10
txi_sub <- list(abundance = txi$abundance[keep, sel], counts = txi$counts[keep, sel],
                length = txi$length[keep, sel], countsFromAbundance = "no")
rm(txi); invisible(gc())
say("Genes kept (>= 10 counts in >= 10 of the 28 samples): ", sum(keep))

dds0 <- DESeqDataSetFromTximport(txi_sub, colData = cdx, design = ~ 1)
nc   <- txi_sub$abundance * rowMeans(txi_sub$length)             # lengthScaledTPM counts for voom
lst  <- t(t(nc) * (colSums(txi_sub$counts) / colSums(nc)))
y    <- calcNormFactors(DGEList(counts = lst))                   # TMM
rm(nc, lst); invisible(gc())

# ---- 2. Random splits (saved, so they can be replayed) ----------------------------
donors <- unique(cdx$donor)
n_per  <- vapply(donors, function(d) sum(cdx$donor == d), numeric(1))
names(n_per) <- donors
design_ok <- function(g2) {
  d <- cdx; d$pseudo <- factor(ifelse(g2, "g2", "g1"), levels = c("g1", "g2"))
  mm <- model.matrix(~ batch_study + sex + pseudo, data = d)
  qr(mm)$rank == ncol(mm)
}
make_donor_split <- function() {
  repeat {
    g2_donors <- donors[runif(length(donors)) < 0.5]
    g2 <- cdx$donor %in% g2_donors
    if (sum(g2) >= 11 && sum(g2) <= 17 && design_ok(g2)) return(g2)
  }
}
make_sample_split <- function() {
  repeat {
    g2 <- rep(FALSE, nrow(cdx)); g2[sample(nrow(cdx), nrow(cdx) / 2)] <- TRUE
    if (design_ok(g2)) return(g2)
  }
}
splits <- list(donor_level  = lapply(seq_len(n_splits_voom), function(i) make_donor_split()),
               sample_level = lapply(seq_len(n_splits_voom), function(i) make_sample_split()))
saveRDS(list(sample_ids = rownames(cdx), splits = splits), out_path("objects", "null_splits.rds"))

# ---- 3. Analyses ---------------------------------------------------------------------
p_stats <- function(p, padj) {
  ok <- is.finite(p)
  data.frame(n_genes = sum(ok), n_sig = sum(padj < 0.05, na.rm = TRUE),
             frac_p05 = mean(p[ok] < 0.05), frac_p001 = mean(p[ok] < 0.001),
             lambda = median(qchisq(p[ok], 1, lower.tail = FALSE)) / qchisq(0.5, 1))
}
voom_run <- function(dmat, blocked) {
  v <- voom(y, dmat)
  cons <- NA_real_
  if (!blocked) {
    fit <- eBayes(lmFit(v, dmat))
  } else {
    dc1 <- duplicateCorrelation(v, dmat, block = cdx$donor)
    v2  <- voom(y, dmat, block = cdx$donor, correlation = dc1$consensus)
    dc2 <- duplicateCorrelation(v2, dmat, block = cdx$donor)
    fit <- eBayes(lmFit(v2, dmat, block = cdx$donor, correlation = dc2$consensus))
    cons <- dc2$consensus.correlation
  }
  p <- fit$p.value[, "pseudog2"]
  cbind(p_stats(p, p.adjust(p, "BH")), consensus_cor = cons)
}
deseq_run <- function(pseudo) {
  d <- dds0
  d$pseudo <- pseudo
  design(d) <- ~ batch_study + sex + pseudo
  d <- DESeq(d, quiet = TRUE)
  r <- results(d, name = "pseudo_g2_vs_g1")
  cbind(p_stats(r$pvalue, r$padj), consensus_cor = NA_real_)
}
safe <- function(label, expr) tryCatch(expr, error = function(e) { say("   ", label, " FAILED: ", conditionMessage(e)); NULL })

rows <- list()
for (type in names(splits)) {
  for (i in seq_along(splits[[type]])) {
    t0  <- Sys.time()
    g2  <- splits[[type]][[i]]
    pseudo <- factor(ifelse(g2, "g2", "g1"), levels = c("g1", "g2"))
    dd  <- cdx; dd$pseudo <- pseudo
    dmat <- model.matrix(~ batch_study + sex + pseudo, data = dd)
    donors_split <- sum(tapply(g2, cdx$donor, function(x) length(unique(x)) > 1))
    base <- data.frame(split_type = type, iter = i, n_g2 = sum(g2), donors_in_both_groups = donors_split)
    res_i <- list()
    res_i[["voom_noblock"]] <- safe("voom_noblock", voom_run(dmat, FALSE))
    res_i[["voom_donor"]]   <- safe("voom_donor",   voom_run(dmat, TRUE))
    if (i <= n_splits_deseq) res_i[["DESeq2"]] <- safe("DESeq2", deseq_run(pseudo))
    res_i <- res_i[!vapply(res_i, is.null, logical(1))]
    if (length(res_i) > 0) {
      rows[[length(rows) + 1]] <- do.call(rbind, lapply(names(res_i), function(m) cbind(base, method = m, res_i[[m]])))
      write.csv(do.call(rbind, rows), out_path("tables", "NULL_calibration_per_split_partial.csv"), row.names = FALSE)
    }
    say(sprintf("[%s] %s split %d/%d | g2 = %d samples, %d donors in both groups | genes padj<0.05: %s | %.0f s",
                format(Sys.time(), "%H:%M:%S"), type, i, length(splits[[type]]), sum(g2), donors_split,
                paste(names(res_i), vapply(res_i, function(x) x$n_sig, numeric(1)), sep = "=", collapse = ", "),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, out_path("tables", "NULL_calibration_per_split.csv"), row.names = FALSE)

# ---- 4. Summary and figures ----------------------------------------------------------
summ <- do.call(rbind, lapply(split(res, list(res$split_type, res$method), drop = TRUE), function(d) data.frame(
  split_type = d$split_type[1], method = d$method[1], n_splits = nrow(d),
  median_genes_padj_lt_0.05 = median(d$n_sig), max_genes_padj_lt_0.05 = max(d$n_sig),
  pct_splits_with_any_DEG = round(100 * mean(d$n_sig > 0), 1),
  mean_fraction_p_lt_0.05 = round(mean(d$frac_p05), 4),
  median_lambda = round(median(d$lambda), 3), min_lambda = round(min(d$lambda), 3), max_lambda = round(max(d$lambda), 3),
  median_donor_correlation = if (all(is.na(d$consensus_cor))) NA else round(median(d$consensus_cor, na.rm = TRUE), 3))))
rownames(summ) <- NULL
show_tab("NULL CALIBRATION SUMMARY (expected under a correct null: lambda ~ 1, fraction p<0.05 ~ 0.05, no DEGs)", summ)
write.csv(summ, out_path("tables", "NULL_calibration_summary.csv"), row.names = FALSE)

res$method_lab <- factor(res$method, levels = c("DESeq2", "voom_noblock", "voom_donor"),
                         labels = c("DESeq2", "limma-voom\n(donors independent)", "limma-voom\n(donor-aware)"))
res$type_lab <- factor(res$split_type, levels = c("donor_level", "sample_level"),
                       labels = c("Donor-level random split\n(donors kept together)",
                                  "Sample-level random split\n(donors ignored; control)"))
theme_n <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "#eeeeee"))
p1 <- ggplot(res, aes(method_lab, lambda)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#898781") +
  geom_boxplot(outlier.shape = NA, fill = "#e1e0d9", colour = "#52514e") +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7, colour = "#2a78d6") +
  facet_wrap(~ type_lab) +
  labs(x = NULL, y = "Inflation factor lambda (1 = well calibrated)",
       title = "Null calibration: random pseudo-groups inside striatum HD samples",
       subtitle = "Each point is one random split; no real disease contrast exists in these comparisons") + theme_n
p2 <- ggplot(res, aes(method_lab, n_sig)) +
  geom_boxplot(outlier.shape = NA, fill = "#e1e0d9", colour = "#52514e") +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7, colour = "#eb6834") +
  scale_y_continuous(transform = scales::pseudo_log_transform(base = 10), breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ type_lab) +
  labs(x = NULL, y = "Genes with padj < 0.05 (false positives)",
       title = "False-positive genes in random pseudo-group comparisons",
       subtitle = "A correct method finds (almost) none") + theme_n
for (nm in c("p1", "p2")) {
  f <- if (nm == "p1") "NULL_lambda_by_method.png" else "NULL_false_positive_genes_by_method.png"
  tryCatch(ggsave(out_path("figures", f), get(nm), width = 9, height = 5.5, dpi = 300),
           error = function(e) say("Figure ", f, " could not be saved: ", conditionMessage(e)))
}
writeLines(capture.output(sessionInfo()), out_path("objects", "sessionInfo.txt"))
say("\nFINISHED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | test_mode = ", test_mode)
say("All outputs: ", run_dir)
if (test_mode) say("This was a dry run. Set test_mode <- FALSE at the top and Source again for the full run.")
say("Please send: the NULL CALIBRATION SUMMARY table from run_log.txt + the two PNG figures.")
