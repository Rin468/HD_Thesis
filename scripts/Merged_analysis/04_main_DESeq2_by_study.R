# =============================================================================
#  04_main_DESeq2_by_study.R
#  MAIN ANALYSIS: one DESeq2 model per study, combined into ONE master table
#
#  Data   Paryani et al. 2024, GSE242195 (76 samples): striatum (caudate +
#         nucleus accumbens) and cingulate cortex, HD and control, 4 batches
#         Lin et al. 2016, GSE79666 (14 samples; labelled "Tseng" in the data
#         objects): motor cortex BA4, 7 HD + 7 control
#
#  Design (fixed before looking at any result of this script)
#    Model P  Paryani only   ~ batch_study + sex + group
#             group = tissue x condition: Striatum_CTRL, Striatum_HD,
#             Cortex_CTRL, Cortex_HD (one combined factor, as in the DESeq2
#             vignette section "Interactions")
#    Model L  GSE79666 only  ~ sex + condition  (~ condition if sex cannot be
#             estimated)
#    Why: DESeq2 estimates ONE dispersion (variability) per gene per model. The
#    two studies differ in laboratory, donor cohort, library preparation and
#    sequencing depth, so a joint model forces one variability estimate on two
#    differently behaving datasets. With separate models the striatum-vs-cortex
#    comparison stays inside one study (same donors, lab and batches) and BA4
#    becomes an independent dataset. The earlier joint 90-sample model
#    (04_final_deseq2_master_table_v2.R) becomes a sensitivity analysis.
#
#  Contrasts   A    striatum: HD vs CTRL                  (Model P)
#              B_P  cingulate cortex: HD vs CTRL          (Model P)
#              B_T  motor cortex BA4: HD vs CTRL          (Model L)
#              C    controls: striatum vs cingulate       (Model P)
#              D    HD: striatum vs cingulate             (Model P)
#              E_P  interaction A - B_P                   (Model P)
#  Between-dataset comparisons (BA4 vs cingulate, striatum vs BA4): difference
#  of two independent estimates, z = (b1 - b2) / sqrt(se1^2 + se2^2),
#  Benjamini-Hochberg adjusted (Altman & Bland 2003, BMJ 326:219).
#  NOTE: these differences mix region AND study (lab, protocol, brain bank).
#
#  Input   (Merged_analysis/results/objects)
#            merged_txi_Paryani76_Tseng14*.rds
#            merged_colData_Paryani76_Tseng14_donorSexRevised*.rds
#            gene_annotation_GENCODEv44*.rds
#  Steps    1. load and check inputs
#           2. tissue verification by marker genes (all 90) + documented corrections
#           3. per-study DESeqDataSets, filtering, sample QC (donors, PCA, variability)
#           4. DESeq2 fits (Model P, Model L) + dispersion comparison
#           5. contrasts, reliability flags, ashr shrinkage
#           6. master table + comparisons between regions/datasets
#           7. summaries, gene lists, figures
#           8. Excel workbook, session info
#  Output   NEW folder Merged_analysis/final_runs/bystudy_<date>_<time>/
#  Run      Session > Restart R (no old workspace loaded), open this file, Source.
# =============================================================================

#  CHANGES 2026-10-15 (final run)
#    - OF045 is EXCLUDED from Model P (mixed tissue: cortical-neuron genes at
#      cortex level and PENK above every cortex sample). It stays in colData_final
#      with region_check = "excluded" and in_model = FALSE, so the QC plots still
#      show it. Set exclude_from_models <- character(0) to keep it.
#    - Model P groups after exclusion: 10 / 27 / 10 / 28 (checked).
#    - The E_P = A - B_P algebra check is now reported on reliable genes only.
# =============================================================================

exclude_from_models <- c("OF045")

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
for (p in c("ashr", "limma")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is missing. Install it first: ",
         if (p == "ashr") "install.packages('ashr')" else "BiocManager::install('limma')")
}

# ---- 0. Paths, run folder, helpers ------------------------------------------
base_dir   <- "C:/!Virag.ORI/Bioinfo/szakdoga"
merged_dir <- file.path(base_dir, "Merged_analysis")
in_dir     <- file.path(merged_dir, "results", "objects")
run_dir    <- file.path(merged_dir, "final_runs", format(Sys.time(), "bystudy_%Y%m%d_%H%M%S"))
for (s in c("tables", "objects", "figures", "gene_lists"))
  dir.create(file.path(run_dir, s), recursive = TRUE, showWarnings = FALSE)
out_path <- function(sub, file) file.path(run_dir, sub, file)
log_file <- file.path(run_dir, "run_log.txt")

say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_file, append = TRUE)
}
show_tab <- function(title, x) {
  out <- c("", paste0("-- ", title, " --"), capture.output(print(x)))
  cat(out, sep = "\n")
  cat(out, sep = "\n", file = log_file, append = TRUE)
}
mem_mb <- function() round(sum(gc()[, 2]))   # memory used by R (Mb) after garbage collection
run_step <- function(name, expr) {           # runs a block; logs the error and stops if it fails
  say("\n==== ", name, " ====  [", format(Sys.time(), "%H:%M:%S"), ", R memory in use: ", mem_mb(), " Mb]")
  tryCatch(expr, error = function(e) {
    say("ERROR in '", name, "': ", conditionMessage(e))
    stop(e)
  })
  invisible(NULL)
}
latest <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) stop("No file matching ", pattern, " in ", dir)
  f[which.max(file.info(f)$mtime)]
}
save_png <- function(plot, file, w = 8.5, h = 6.5) {
  tryCatch(ggsave(out_path("figures", file), plot, width = w, height = h, dpi = 300),
           error = function(e) say("Figure ", file, " could not be saved: ", conditionMessage(e)))
}
scor <- function(x, y) {                     # Spearman correlation, NA if too few pairs
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  round(cor(x[ok], y[ok], method = "spearman"), 3)
}
maxdiff <- function(x) { x <- abs(x[is.finite(x)]); if (length(x)) signif(max(x), 3) else NA }
isT <- function(x) x %in% TRUE               # TRUE stays TRUE, FALSE and NA become FALSE
strip_version <- function(x) unique(sub("\\..*$", "", x))
theme_thesis <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        panel.border     = element_rect(colour = "#c3c2b7"),
        axis.text        = element_text(colour = "#52514e"),
        plot.subtitle    = element_text(colour = "#52514e"))
col_tissue    <- c(Striatum = "#eb6834", Cortex = "#2a78d6")
col_condition <- c(CTRL = "#898781", HD = "#2a78d6")

# PCA on the 500 most variable genes: R^2 table, coordinates, plot(s), robust outlier check
pca_qc <- function(mat, meta, tag, colour_by, colour_values, shape_by, r2_vars,
                   label_rows = NULL, facet_by = NULL, group_by = NULL) {
  top <- order(apply(mat, 1, var), decreasing = TRUE)[seq_len(min(500, nrow(mat)))]
  pc  <- prcomp(t(mat[top, ]))
  pv  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  npc <- min(6, ncol(pc$x))
  df  <- data.frame(meta, pc$x[, seq_len(npc), drop = FALSE], check.names = FALSE)
  r2_vars <- r2_vars[vapply(r2_vars, function(v) length(unique(df[[v]])) > 1, logical(1))]
  r2 <- vapply(r2_vars, function(v) vapply(seq_len(npc), function(k)
    summary(lm(df[[paste0("PC", k)]] ~ df[[v]]))$r.squared, numeric(1)), numeric(npc))
  r2 <- matrix(r2, nrow = npc,
               dimnames = list(sprintf("PC%d (%.1f%%)", seq_len(npc), pv[seq_len(npc)]), r2_vars))
  show_tab(paste0("R^2 of PC1-", npc, " with covariates (", tag, ")"), round(r2, 2))
  write.csv(data.frame(PC = rownames(r2), round(r2, 3), check.names = FALSE),
            out_path("tables", paste0("QC_PCA_R2_", tag, ".csv")), row.names = FALSE)
  write.csv(df, out_path("tables", paste0("QC_PCA_coordinates_", tag, ".csv")), row.names = FALSE)
  
  p <- ggplot(df, aes(PC1, PC2, colour = .data[[colour_by]], shape = .data[[shape_by]])) +
    geom_point(size = 2.6, alpha = 0.9) +
    scale_colour_manual(values = colour_values) +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
         colour = colour_by, shape = shape_by,
         title = paste0("PCA, 500 most variable genes (", tag, ")")) +
    theme_thesis
  if (!is.null(label_rows) && any(label_rows)) {
    p <- p + geom_text(data = df[label_rows, , drop = FALSE], aes(label = sample_label),
                       size = 2.3, vjust = -1, colour = "#52514e", show.legend = FALSE)
  }
  save_png(p, paste0("QC_PCA_", tag, ".png"))
  if (!is.null(facet_by)) {
    save_png(p + facet_wrap(stats::as.formula(paste("~", facet_by))),
             paste0("QC_PCA_", tag, "_by_", facet_by, ".png"), w = 10, h = 7)
  }
  if (!is.null(group_by)) {                  # robust z (median/MAD) within group on PC1 and PC2
    rz <- function(x) {
      s <- mad(x)
      if (!is.finite(s) || s == 0) return(rep(0, length(x)))
      (x - median(x)) / s
    }
    z1  <- ave(df$PC1, df[[group_by]], FUN = rz)
    z2  <- ave(df$PC2, df[[group_by]], FUN = rz)
    out <- df[abs(z1) > 3.5 | abs(z2) > 3.5, c("sample_label", group_by, "PC1", "PC2"), drop = FALSE]
    if (nrow(out) == 0) {
      say("Possible outliers in ", tag, " (robust z > 3.5 within ", group_by, " on PC1 or PC2): none")
    } else {
      show_tab(paste0("Possible outliers in ", tag, " (robust z > 3.5 within ", group_by,
                      " on PC1 or PC2) - inspect, do not remove automatically"), out)
    }
  }
  invisible(list(df = df, pv = pv))
}

# Build a DESeqDataSet for a subset of samples (tximport offsets kept) and filter genes
build_dds <- function(sel, design, min_samples, tag) {
  cdx <- droplevels(cd[sel, ])
  mm  <- model.matrix(design, data = cdx)
  if (qr(mm)$rank < ncol(mm)) stop("Design for ", tag, " is not full rank")
  txi_sub <- list(abundance = txi$abundance[, sel], counts = txi$counts[, sel],
                  length = txi$length[, sel], countsFromAbundance = "no")
  dds  <- DESeqDataSetFromTximport(txi_sub, colData = cdx, design = design)
  keep <- rowSums(counts(dds) >= 10) >= min_samples
  say(tag, ": ", ncol(dds), " samples | design ", paste(deparse(design), collapse = ""),
      " (", ncol(mm), " coefficients, full rank) | genes kept: ", sum(keep), " of ", length(keep),
      " (>= 10 counts in >= ", min_samples, " samples)")
  dds[keep, ]
}

# Fit DESeq2, remove non-converged genes, save the object and the dispersion plot
fit_deseq <- function(dds, tag, file_tag) {
  dds <- DESeq(dds)
  if ("replace" %in% names(mcols(dds)))
    say(tag, ": genes with an outlier count replaced (Cook's distance): ", sum(mcols(dds)$replace, na.rm = TRUE))
  bad <- !(mcols(dds)$betaConv %in% TRUE)
  say(tag, ": non-converged genes removed: ", sum(bad),
      if (any(bad) && sum(bad) <= 20) paste0(" (", paste(rownames(dds)[bad], collapse = ", "), ")") else "")
  dds <- dds[!bad, ]
  say(tag, ": resultsNames = ", paste(resultsNames(dds), collapse = ", "))
  saveRDS(dds, out_path("objects", paste0("dds_", file_tag, ".rds")))
  png(out_path("figures", paste0("QC_dispersion_", file_tag, ".png")), width = 1800, height = 1500, res = 250)
  tryCatch(plotDispEsts(dds, main = paste0("Dispersion estimates: ", tag)),
           error = function(e) say("Dispersion plot failed: ", conditionMessage(e)),
           finally = dev.off())
  dds
}

# One contrast: numeric contrast vector, reliability flags, ashr shrinkage
run_contrast <- function(k) {
  s   <- specs[[k]]
  dds <- if (s$model == "P") dds_P else dds_L
  rn  <- resultsNames(dds)
  if (!all(names(s$w) %in% rn))
    stop("Contrast ", k, ": coefficient not found: ", paste(setdiff(names(s$w), rn), collapse = ", "))
  cv <- setNames(numeric(length(rn)), rn)
  cv[names(s$w)] <- s$w
  n_genes <- nrow(dds)
  res <- tryCatch(results(dds, contrast = cv, alpha = 0.05),
                  error = function(e) { say("  ", k, ": results() FAILED - ", conditionMessage(e)); NULL })
  if (is.null(res)) {
    return(data.frame(log2FC = rep(NA_real_, n_genes), log2FC_shr = NA_real_, lfcSE = NA_real_,
                      pvalue = NA_real_, padj = NA_real_, reliable = FALSE, sig = FALSE,
                      row.names = rownames(dds)))
  }
  cnt      <- counts(dds)
  zero     <- Reduce(`|`, lapply(s$sets, function(set) rowSums(cnt[, set, drop = FALSE]) == 0))
  extreme  <- !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 10
  reliable <- !zero & !extreme
  
  # ashr shrinkage on estimable genes only (the call lfcShrink(type = "ashr") makes)
  shr <- rep(NA_real_, n_genes)
  ok  <- reliable & is.finite(res$log2FoldChange) & is.finite(res$lfcSE) & isT(res$lfcSE > 0)
  fit <- tryCatch(ashr::ash(res$log2FoldChange[ok], res$lfcSE[ok], mixcompdist = "normal", method = "shrink"),
                  error = function(e) { say("  ", k, ": ashr FAILED - ", conditionMessage(e)); NULL })
  if (!is.null(fit)) shr[ok] <- fit$result$PosteriorMean
  
  out <- data.frame(log2FC = res$log2FoldChange, log2FC_shr = shr, lfcSE = res$lfcSE,
                    pvalue = res$pvalue, padj = res$padj, reliable = reliable,
                    sig = reliable & isT(res$padj < 0.05), row.names = rownames(res))
  say(sprintf("%-4s %-58s significant: %5d (up %5d, down %5d) | median SE: %.3f | shrunk genes: %5d",
              k, s$label, sum(out$sig), sum(out$sig & isT(out$log2FC > 0)), sum(out$sig & isT(out$log2FC < 0)),
              median(res$lfcSE[ok]), sum(ok)))
  out
}

say("Run folder: ", run_dir)
say("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | 04_main_DESeq2_by_study.R")
start_mb <- mem_mb()
say("R memory in use before the analysis: ", start_mb, " Mb")
if (start_mb > 500)
  say("WARNING: ", start_mb, " Mb already in use - an old workspace is probably loaded. ",
      "Recommended: Tools > Global Options > General: untick 'Restore .RData', then Session > Restart R.")


# ---- 1. Load and check the inputs -------------------------------------------
run_step("1. Load inputs", {
  f_txi <- latest(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$")
  f_cd  <- latest(in_dir, "^merged_colData_Paryani76_Tseng14_donorSexRevised(_v[0-9]+)?\\.rds$")
  f_ann <- latest(in_dir, "^gene_annotation_GENCODEv44(_v[0-9]+)?\\.rds$")
  say("Merged txi : ", f_txi)
  say("colData    : ", f_cd)
  say("Annotation : ", f_ann)
  txi <- readRDS(f_txi)
  cd  <- as.data.frame(readRDS(f_cd))
  ann <- as.data.frame(readRDS(f_ann))
  stopifnot(identical(colnames(txi$counts), rownames(cd)),
            identical(rownames(txi$counts), ann$gene_id),
            identical(txi$countsFromAbundance, "no"),
            nrow(cd) == 90, !anyNA(cd$sex), !anyNA(cd$group), !anyNA(cd$batch_study),
            all(c("sample_id", "sample_label", "donor", "study", "region", "tissue", "condition")
                %in% names(cd)))
  cd$condition   <- factor(as.character(cd$condition), levels = c("CTRL", "HD"))
  cd$tissue      <- factor(as.character(cd$tissue), levels = c("Striatum", "Cortex"))
  cd$sex         <- factor(as.character(cd$sex))
  cd$batch_study <- factor(as.character(cd$batch_study))
  if (anyNA(cd$condition) || anyNA(cd$tissue)) stop("condition must be CTRL/HD and tissue Striatum/Cortex")
  if (sum(cd$study == "Paryani") != 76 || sum(cd$study == "Tseng") != 14)
    stop("Expected 76 Paryani and 14 GSE79666 ('Tseng') samples in the column 'study'")
  show_tab("Samples per group x study (GEO labels)", table(cd$group, cd$study))
})


# ---- 2. Tissue verification with marker genes + documented corrections ------
run_step("2. Tissue verification (all 90 samples) and corrections", {
  dds0   <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ 1)
  keep90 <- rowSums(counts(dds0) >= 10) >= 10
  dds0   <- dds0[keep90, ]
  say("Genes kept for the joint QC transformation: ", sum(keep90), " of ", length(keep90))
  vsd90 <- vst(dds0, blind = TRUE)            # blind: QC must not use the design
  rm(dds0); invisible(gc())
  sym <- ann$gene_symbol[match(rownames(vsd90), ann$gene_id)]
  
  striatal_markers <- c("PPP1R1B", "PENK", "DRD1", "DRD2", "ADORA2A", "GPR88", "RGS9", "TAC1")
  cortical_markers <- c("SLC17A7", "SATB2", "NEUROD6", "TBR1")   # absent from striatum with or without HD
  if (!all(c(striatal_markers, cortical_markers) %in% sym)) stop("Marker genes missing after filtering")
  cd$striatal_score <- round(colMeans(assay(vsd90)[sym %in% striatal_markers, , drop = FALSE]), 2)
  cd$cortical_score <- round(colMeans(assay(vsd90)[sym %in% cortical_markers, , drop = FALSE]), 2)
  
  # Reference ranges from batches 1-3 + GSE79666 (their labels behaved as expected in the PCA)
  ref   <- cd$batch_study != "Paryani_b4"
  thr   <- mean(c(median(cd$cortical_score[ref & cd$tissue == "Striatum"]),
                  median(cd$cortical_score[ref & cd$tissue == "Cortex"])))
  agree <- ifelse(cd$cortical_score > thr, "Cortex", "Striatum") == as.character(cd$tissue)
  rng_str <- range(cd$cortical_score[ref & agree & cd$tissue == "Striatum"])
  rng_ctx <- range(cd$cortical_score[ref & agree & cd$tissue == "Cortex"])
  say("Cortical-marker score: reference striatum ", rng_str[1], " to ", rng_str[2],
      " | reference cortex ", rng_ctx[1], " to ", rng_ctx[2], " | midpoint ", round(thr, 2))
  
  # Documented correction table (decided from the marker check, script 03b)
  fix <- data.frame(
    sample_label = c("OF021", "OF024", "OF038", "OF049", "H5608_Accum", "OF025", "OF050",
                     "H5263_Cing", "OF045", "OF033"),
    tissue_new   = c("Cortex", "Cortex", "Cortex", "Cortex", "Cortex", "Striatum", "Striatum",
                     "Cortex", "Striatum", "Striatum"),
    action       = c(rep("relabelled", 7), "kept_flagged",
                     if ("OF045" %in% exclude_from_models) "excluded" else "kept_flagged", "kept_flagged"),
    reason       = c(rep("cortical-neuron marker score within the cortex reference range", 5),
                     rep("cortical-neuron marker score within the striatum reference range", 2),
                     "intermediate cortical score with low striatal score: consistent with cortex",
                     "high striatal and cortical scores: mixed tissue, identity cannot be assigned",
                     "intermediate cortical score with high striatal score: possible cortical admixture"),
    stringsAsFactors = FALSE)
  if (!all(exclude_from_models %in% fix$sample_label))
    stop("exclude_from_models contains a sample that is not in the correction table")
  i <- match(fix$sample_label, cd$sample_label)
  if (anyNA(i)) stop("Correction table contains unknown sample labels")
  fix$sample_id      <- cd$sample_id[i]
  fix$condition      <- as.character(cd$condition[i])
  fix$region_GEO     <- as.character(cd$region[i])
  fix$tissue_GEO     <- as.character(cd$tissue[i])
  fix$striatal_score <- cd$striatal_score[i]
  fix$cortical_score <- cd$cortical_score[i]
  
  rel    <- fix$action == "relabelled"
  ok_ctx <- fix$cortical_score[rel & fix$tissue_new == "Cortex"]   >= rng_ctx[1] - 0.5
  ok_str <- fix$cortical_score[rel & fix$tissue_new == "Striatum"] <= rng_str[2] + 0.5
  if (!all(ok_ctx) || !all(ok_str)) stop("Marker scores no longer support a planned relabelling - send run_log.txt")
  if (!all(fix$tissue_new[rel] != fix$tissue_GEO[rel]) || !all(fix$tissue_new[!rel] == fix$tissue_GEO[!rel]))
    stop("Correction table is inconsistent with the GEO labels")
  unexpected <- setdiff(cd$sample_label[!agree], fix$sample_label)
  if (length(unexpected) > 0)
    stop("Marker/label mismatch not covered by the correction table: ", paste(unexpected, collapse = ", "))
  
  # Apply: tissue/group corrected; GEO labels kept in separate columns for transparency
  cd$tissue_GEO <- cd$tissue
  cd$group_GEO  <- cd$group
  tis <- as.character(cd$tissue)
  tis[i[rel]] <- fix$tissue_new[rel]
  cd$tissue <- factor(tis, levels = c("Striatum", "Cortex"))
  cd$group  <- factor(paste(cd$tissue, cd$condition, sep = "_"),
                      levels = c("Striatum_CTRL", "Striatum_HD", "Cortex_CTRL", "Cortex_HD"))
  cd$region_final <- as.character(cd$region)
  cd$region_final[i[rel]] <- paste0(fix$tissue_new[rel], "_unspecified")
  cd$region_check <- "consistent"
  cd$region_check[i] <- fix$action
  cd$exclude_in_sensitivity <- cd$region_check != "consistent"
  cd$in_model <- !(cd$sample_label %in% exclude_from_models)
  say("Samples excluded from the DESeq2 models: ",
      if (any(!cd$in_model)) paste(cd$sample_label[!cd$in_model], collapse = ", ") else "none")

  show_tab("Tissue corrections (90 samples in colData; excluded samples are not modelled)",
           fix[, c("sample_label", "condition", "region_GEO", "tissue_GEO", "tissue_new", "action",
                   "striatal_score", "cortical_score")])
  show_tab("Group: GEO label (rows) vs corrected (columns)", table(cd$group_GEO, cd$group))
  show_tab("Corrected group x study", table(cd$group, cd$study))
  show_tab("Sex (expression-based) x corrected group x study", ftable(table(cd$study, cd$group, cd$sex)))
  if (!all(as.vector(table(cd$group)) == c(10, 28, 17, 35))) stop("Unexpected group sizes after correction")
  
  saveRDS(cd, out_path("objects", "colData_final.rds"))
  write.csv(cd, out_path("tables", "colData_final.csv"), row.names = FALSE)
  write.csv(fix, out_path("tables", "region_corrections.csv"), row.names = FALSE)
  
  mk <- data.frame(sample_label = cd$sample_label, tissue = cd$tissue, check = cd$region_check,
                   striatal = cd$striatal_score, cortical = cd$cortical_score)
  p_mk <- ggplot(mk, aes(cortical, striatal, colour = tissue, shape = check)) +
    geom_vline(xintercept = c(rng_str[2], rng_ctx[1]), linetype = "dashed", colour = "#898781") +
    geom_point(size = 2.6, alpha = 0.9) +
    geom_text(data = mk[mk$check != "consistent", ], aes(label = sample_label),
              size = 2.4, vjust = -1, colour = "#52514e", show.legend = FALSE) +
    scale_colour_manual(values = col_tissue) +
    labs(x = "Cortical neuron markers, mean VST (SLC17A7, SATB2, NEUROD6, TBR1)",
         y = "Striatal MSN markers, mean VST",
         colour = "Tissue (corrected)", shape = "Tissue check",
         title = "Tissue identity of every sample by marker genes",
         subtitle = "Dashed lines: upper limit of reference striatum / lower limit of reference cortex") +
    theme_thesis
  save_png(p_mk, "QC_tissue_markers_corrected.png")
  
  # Overview of all 90 samples (shows the study effect that motivates separate models)
  pca_qc(assay(vsd90), cd, "all90_uncorrected", colour_by = "tissue", colour_values = col_tissue,
         shape_by = "study", r2_vars = c("study", "batch_study", "tissue", "condition", "group", "sex"),
         label_rows = cd$region_check != "consistent", facet_by = "batch_study")
  rm(vsd90); invisible(gc())
})


# ---- 3. Per-study DESeqDataSets and sample QC --------------------------------
run_step("3. Per-study DESeqDataSets and sample QC", {
  selP <- cd$study == "Paryani" & cd$in_model
  selL <- cd$study == "Tseng" & cd$in_model
  if (any(cd$tissue[selL] != "Cortex")) stop("All GSE79666 samples should be cortex (BA4)")
  
  # Model L: adjust for sex only if both sexes are present and sex is not confounded with condition
  cdL0 <- droplevels(cd[selL, ])
  show_tab("GSE79666: sex x condition", table(cdL0$sex, cdL0$condition))
  design_L <- ~ sex + condition
  if (nlevels(cdL0$sex) < 2) {
    design_L <- ~ condition
    say("GSE79666: only one sex present - design reduced to ~ condition")
  } else {
    mmL <- model.matrix(design_L, data = cdL0)
    if (qr(mmL)$rank < ncol(mmL)) {
      design_L <- ~ condition
      say("WARNING: sex is confounded with condition in GSE79666 - design reduced to ~ condition")
    }
  }
  design_P <- ~ batch_study + sex + group
  
  dds_P <- build_dds(selP, design_P, min_samples = 10, tag = "Model P (Paryani)")
  dds_L <- build_dds(selL, design_L, min_samples = 7,  tag = "Model L (GSE79666)")
  rm(txi); invisible(gc())                   # the merged raw object is no longer needed
  cdP <- as.data.frame(colData(dds_P))
  cdL <- as.data.frame(colData(dds_L))
  expected_P <- c(Striatum_CTRL = 10, Striatum_HD = 28, Cortex_CTRL = 10, Cortex_HD = 28)
  expected_P["Striatum_HD"] <- 28 - sum(c("OF045", "OF033") %in% exclude_from_models)
  if (!all(table(cdP$group)[names(expected_P)] == expected_P))
    stop("Unexpected Model P group sizes: ", paste(names(table(cdP$group)), table(cdP$group), collapse = ", "))
  say("Model P group sizes: ", paste(names(expected_P), expected_P, sep = " = ", collapse = ", "))
  show_tab("Model P: samples per group x batch", table(cdP$group, cdP$batch_study))
  show_tab("Model P: samples per group x sex", table(cdP$group, cdP$sex))
  
  # Donor structure in Paryani: several samples from one donor within a tissue are not independent
  rep_tab <- do.call(rbind, lapply(split(cdP, list(cdP$tissue, cdP$condition), drop = TRUE), function(x) {
    rt <- table(x$region_final)
    data.frame(tissue = as.character(x$tissue[1]), condition = as.character(x$condition[1]),
               samples = nrow(x), donors = length(unique(x$donor)),
               donors_with_2plus_samples = sum(table(as.character(x$donor)) > 1),
               regions = paste(names(rt), as.vector(rt), sep = ": ", collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  rownames(rep_tab) <- NULL
  show_tab("Model P: samples and donors per tissue x condition (donor IDs as reconstructed in 01b)", rep_tab)
  write.csv(rep_tab, out_path("tables", "QC_donor_structure_Paryani.csv"), row.names = FALSE)
  
  # Blind VST per study (QC only)
  vst_P <- vst(dds_P, blind = TRUE)
  vst_L <- vst(dds_L, blind = TRUE)
  saveRDS(vst_P, out_path("objects", "vst_blind_ModelP_Paryani.rds"))
  saveRDS(vst_L, out_path("objects", "vst_blind_ModelL_GSE79666.rds"))
  mP <- assay(vst_P)
  mL <- assay(vst_L)
  mP_corr <- limma::removeBatchEffect(mP, batch = cdP$batch_study,
                                      covariates = model.matrix(~ sex, data = cdP)[, -1, drop = FALSE],
                                      design = model.matrix(~ group, data = cdP))
  sex_in_L <- "sex" %in% all.vars(design_L)
  mL_corr <- if (sex_in_L) {
    limma::removeBatchEffect(mL, covariates = model.matrix(~ sex, data = cdL)[, -1, drop = FALSE],
                             design = model.matrix(~ condition, data = cdL))
  } else mL
  
  pca_qc(mP, cdP, "Paryani_uncorrected", colour_by = "tissue", colour_values = col_tissue,
         shape_by = "condition", r2_vars = c("batch_study", "tissue", "condition", "group", "sex"),
         label_rows = cdP$region_check != "consistent", facet_by = "batch_study")
  pcP <- pca_qc(mP_corr, cdP, "Paryani_batch_sex_removed_display_only", colour_by = "tissue",
                colour_values = col_tissue, shape_by = "condition",
                r2_vars = c("batch_study", "tissue", "condition", "group", "sex"),
                label_rows = cdP$region_check != "consistent", group_by = "group")
  if ("grade" %in% names(cdP)) {
    p_gr <- ggplot(pcP$df, aes(PC1, PC2, colour = factor(grade), shape = tissue)) +
      geom_point(size = 2.6, alpha = 0.9) +
      scale_colour_viridis_d(end = 0.9, na.value = "#c3c2b7") +
      labs(x = sprintf("PC1 (%.1f%%)", pcP$pv[1]), y = sprintf("PC2 (%.1f%%)", pcP$pv[2]),
           colour = "Grade", shape = "Tissue",
           title = "PCA, Paryani (batch and sex removed for display) - coloured by grade") +
      theme_thesis
    save_png(p_gr, "QC_PCA_Paryani_by_grade.png")
  }
  pca_qc(mL, cdL, "GSE79666_uncorrected", colour_by = "condition", colour_values = col_condition,
         shape_by = "sex", r2_vars = c("condition", "sex"),
         label_rows = rep(TRUE, nrow(cdL)), group_by = "condition")
  
  # Within-group variability (descriptive): median per-gene SD of VST values in each group,
  # genes tested in both studies and expressed above the median in both
  common_vst <- intersect(rownames(mP), rownames(mL))
  hiP   <- rowMeans(mP)[common_vst] >= median(rowMeans(mP))
  hiL   <- rowMeans(mL)[common_vst] >= median(rowMeans(mL))
  g_use <- common_vst[hiP & hiL]
  wsd <- function(m, grp) vapply(split(seq_along(grp), grp), function(ix)
    median(apply(m[g_use, ix, drop = FALSE], 1, sd)), numeric(1))
  vP <- wsd(mP_corr, cdP$group)
  vL <- wsd(mL_corr, cdL$condition)
  var_tab <- data.frame(
    dataset = c(rep("Paryani (batch + sex removed)", length(vP)),
                rep(if (sex_in_L) "GSE79666 (sex removed)" else "GSE79666", length(vL))),
    group   = c(names(vP), paste0("BA4_", names(vL))),
    samples = c(as.vector(table(cdP$group)[names(vP)]), as.vector(table(cdL$condition)[names(vL)])),
    median_within_group_SD = round(c(vP, vL), 3), stringsAsFactors = FALSE)
  show_tab(paste0("Within-group variability, median per-gene SD of VST values (", length(g_use),
                  " genes expressed in both studies; descriptive)"), var_tab)
  write.csv(var_tab, out_path("tables", "QC_within_group_variability.csv"), row.names = FALSE)
  rm(vst_P, vst_L, mP, mL, mP_corr, mL_corr, pcP); invisible(gc())
})


# ---- 4. DESeq2 fits + dispersion comparison ----------------------------------
run_step("4. DESeq2 fits (Model P and Model L)", {
  dds_P <- fit_deseq(dds_P, "Model P (Paryani)", "ModelP_Paryani")
  dds_L <- fit_deseq(dds_L, "Model L (GSE79666)", "ModelL_GSE79666")
  
  design_tab <- data.frame(
    model   = c("P", "L"),
    dataset = c("Paryani et al. 2024 (GSE242195)", "Lin et al. 2016 (GSE79666)"),
    samples = c(ncol(dds_P), ncol(dds_L)),
    groups  = c(paste(names(table(dds_P$group)), as.vector(table(dds_P$group)), sep = ": ", collapse = "; "),
                paste(names(table(dds_L$condition)), as.vector(table(dds_L$condition)), sep = ": ", collapse = "; ")),
    design  = c(paste(deparse(design_P), collapse = ""), paste(deparse(design_L), collapse = "")),
    genes_tested = c(nrow(dds_P), nrow(dds_L)), stringsAsFactors = FALSE)
  show_tab("Models", design_tab)
  write.csv(design_tab, out_path("tables", "MODELS_design_summary.csv"), row.names = FALSE)
  
  # Do the two studies differ in variability? (the reason for separate models)
  common <- intersect(rownames(dds_P), rownames(dds_L))
  iP <- match(common, rownames(dds_P))
  iL <- match(common, rownames(dds_L))
  gP <- mcols(dds_P)$dispGeneEst[iP]
  gL <- mcols(dds_L)$dispGeneEst[iL]
  okd <- is.finite(gP) & is.finite(gL) & gP > 1e-6 & gL > 1e-6
  asym <- function(d) {
    co <- attr(dispersionFunction(d), "coefficients")
    if (is.null(co)) NA_real_ else signif(unname(co[1]), 3)
  }
  aP <- asym(dds_P); aL <- asym(dds_L)
  disp_tab <- data.frame(
    measure = c("Genes tested in both models",
                "Median gene-wise dispersion, Paryani",
                "Median gene-wise dispersion, GSE79666",
                "Median final dispersion, Paryani",
                "Median final dispersion, GSE79666",
                "Asymptotic dispersion of the fitted trend, Paryani",
                "Asymptotic dispersion of the fitted trend, GSE79666",
                "Biological CV at high counts (sqrt of asymptotic dispersion), Paryani",
                "Biological CV at high counts (sqrt of asymptotic dispersion), GSE79666",
                "Median per-gene ratio of gene-wise dispersions, GSE79666 / Paryani",
                "Genes with higher gene-wise dispersion in GSE79666 (%)"),
    value = c(length(common),
              signif(median(gP[okd]), 3), signif(median(gL[okd]), 3),
              signif(median(dispersions(dds_P)[iP], na.rm = TRUE), 3),
              signif(median(dispersions(dds_L)[iL], na.rm = TRUE), 3),
              aP, aL, round(sqrt(aP), 3), round(sqrt(aL), 3),
              signif(median(gL[okd] / gP[okd]), 3),
              round(100 * mean(gL[okd] > gP[okd]), 1)),
    stringsAsFactors = FALSE)
  show_tab("Dispersion (variability) of the two studies, genes tested in both", disp_tab)
  write.csv(disp_tab, out_path("tables", "QC_dispersion_comparison.csv"), row.names = FALSE)
  
  dd <- data.frame(P = gP[okd], L = gL[okd])
  p_d <- ggplot(dd, aes(P, L)) +
    geom_point(size = 0.4, alpha = 0.25, colour = "#52514e") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#eb6834") +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Gene-wise dispersion, Paryani (Model P)", y = "Gene-wise dispersion, GSE79666 (Model L)",
         title = "Gene-wise dispersion estimates of the two studies",
         subtitle = paste0(nrow(dd), " genes tested in both; dashed line = equal dispersion; median ratio L/P = ",
                           signif(median(dd$L / dd$P), 3))) +
    theme_thesis
  save_png(p_d, "QC_dispersion_comparison.png")
})


# ---- 5. Contrasts, reliability flags, shrinkage -------------------------------
run_step("5. Contrasts, flags and shrinkage", {
  cdP <- as.data.frame(colData(dds_P))
  cdL <- as.data.frame(colData(dds_L))
  P_SC <- cdP$group == "Striatum_CTRL"; P_SH <- cdP$group == "Striatum_HD"
  P_CC <- cdP$group == "Cortex_CTRL";   P_CH <- cdP$group == "Cortex_HD"
  L_C  <- cdL$condition == "CTRL";      L_H  <- cdL$condition == "HD"
  
  g_SHD <- "group_Striatum_HD_vs_Striatum_CTRL"
  g_CC  <- "group_Cortex_CTRL_vs_Striatum_CTRL"
  g_CHD <- "group_Cortex_HD_vs_Striatum_CTRL"
  c_HD  <- "condition_HD_vs_CTRL"
  
  # w = weights on the model coefficients; sets = groups that must not be all-zero
  specs <- list(
    A   = list(model = "P", w = setNames(1, g_SHD), sets = list(P_SH, P_SC),
               label = "Striatum (Paryani): HD vs CTRL"),
    B_P = list(model = "P", w = setNames(c(1, -1), c(g_CHD, g_CC)), sets = list(P_CH, P_CC),
               label = "Cingulate cortex (Paryani): HD vs CTRL"),
    B_T = list(model = "L", w = setNames(1, c_HD), sets = list(L_H, L_C),
               label = "Motor cortex BA4 (GSE79666): HD vs CTRL"),
    C   = list(model = "P", w = setNames(-1, g_CC), sets = list(P_SC, P_CC),
               label = "Controls (Paryani): striatum vs cingulate"),
    D   = list(model = "P", w = setNames(c(1, -1), c(g_SHD, g_CHD)), sets = list(P_SH, P_CH),
               label = "HD (Paryani): striatum vs cingulate"),
    E_P = list(model = "P", w = setNames(c(1, -1, 1), c(g_SHD, g_CHD, g_CC)),
               sets = list(P_SH, P_SC, P_CH, P_CC),
               label = "Interaction (Paryani): striatum minus cingulate HD effect"))
  
  res_list <- list()
  for (k in names(specs)) {
    res_list[[k]] <- run_contrast(k)
    invisible(gc())
  }
  
  # Algebra checks on unshrunken fold changes (both should be ~0)
  L <- function(k) res_list[[k]]$log2FC
  relE <- res_list$A$reliable & res_list$B_P$reliable & res_list$E_P$reliable
  say("Check E_P = A - B_P (reliable genes, n = ", sum(relE), "), max |diff|: ",
      maxdiff((L("E_P") - (L("A") - L("B_P")))[relE]), "   <- should be ~0")
  say("Check E_P = A - B_P (all genes, incl. non-estimable), max |diff|: ",
      maxdiff(L("E_P") - (L("A") - L("B_P"))), "   (non-zero is expected only from flagged genes)")
  say("Check E_P = D - C,    max |diff|: ", maxdiff(L("E_P") - (L("D") - L("C"))))
  saveRDS(res_list, out_path("objects", "results_by_contrast.rds"))
  saveRDS(specs, out_path("objects", "contrast_specifications.rds"))
})


# ---- 6. MASTER TABLE -----------------------------------------------------------
run_step("6. Master table", {
  genes_P   <- rownames(dds_P)
  genes_L   <- rownames(dds_L)
  all_genes <- ann$gene_id[ann$gene_id %in% c(genes_P, genes_L)]
  mt <- data.frame(gene_id     = all_genes,
                   gene_symbol = ann$gene_symbol[match(all_genes, ann$gene_id)],
                   stringsAsFactors = FALSE)
  if ("gene_type" %in% names(ann)) mt$gene_type <- ann$gene_type[match(all_genes, ann$gene_id)]
  mt$tested_Paryani    <- mt$gene_id %in% genes_P
  mt$tested_GSE79666   <- mt$gene_id %in% genes_L
  mt$baseMean_Paryani  <- round(mcols(dds_P)$baseMean[match(mt$gene_id, genes_P)], 2)
  mt$baseMean_GSE79666 <- round(mcols(dds_L)$baseMean[match(mt$gene_id, genes_L)], 2)
  
  getv <- function(k, col) res_list[[k]][[col]][match(mt$gene_id, rownames(res_list[[k]]))]
  for (k in names(res_list)) {
    mt[[paste0("log2FC_", k)]]     <- round(getv(k, "log2FC"), 4)
    mt[[paste0("log2FC_shr_", k)]] <- round(getv(k, "log2FC_shr"), 4)
    mt[[paste0("lfcSE_", k)]]      <- round(getv(k, "lfcSE"), 4)
    mt[[paste0("pvalue_", k)]]     <- getv(k, "pvalue")
    mt[[paste0("padj_", k)]]       <- getv(k, "padj")
    mt[[paste0("reliable_", k)]]   <- isT(getv(k, "reliable"))
    mt[[paste0("sig_", k)]]        <- isT(getv(k, "sig"))
  }
  
  # Difference of two independent estimates (different studies, different donors)
  diff_test <- function(k1, k2) {
    l1 <- getv(k1, "log2FC"); s1 <- getv(k1, "lfcSE")
    l2 <- getv(k2, "log2FC"); s2 <- getv(k2, "lfcSE")
    ok <- isT(getv(k1, "reliable")) & isT(getv(k2, "reliable")) &
      is.finite(l1) & is.finite(l2) & is.finite(s1) & is.finite(s2) & isT(s1 > 0) & isT(s2 > 0)
    z <- rep(NA_real_, length(l1)); p <- z; q <- z
    z[ok] <- (l1[ok] - l2[ok]) / sqrt(s1[ok]^2 + s2[ok]^2)
    p[ok] <- 2 * pnorm(-abs(z[ok]))
    q[ok] <- p.adjust(p[ok], method = "BH")
    list(diff = l1 - l2, z = z, p = p, q = q, sig = ok & isT(q < 0.05), n = sum(ok))
  }
  dC <- diff_test("B_T", "B_P")              # BA4 minus cingulate
  dS <- diff_test("A", "B_T")                # striatum minus BA4
  mt$diffLFC_BA4_minus_cingulate  <- round(dC$diff, 4)
  mt$diff_z_BA4_vs_cingulate      <- round(dC$z, 3)
  mt$diff_pvalue_BA4_vs_cingulate <- dC$p
  mt$diff_padj_BA4_vs_cingulate   <- dC$q
  mt$cortex_datasets_differ       <- dC$sig
  mt$diffLFC_striatum_minus_BA4   <- round(dS$diff, 4)
  mt$diff_z_striatum_vs_BA4       <- round(dS$z, 3)
  mt$diff_pvalue_striatum_vs_BA4  <- dS$p
  mt$diff_padj_striatum_vs_BA4    <- dS$q
  mt$striatum_BA4_differ          <- dS$sig
  say("Difference tests of independent estimates: BA4 vs cingulate ", dC$n, " genes, striatum vs BA4 ",
      dS$n, " genes")
  
  # HD response categories for two regions + verdict that requires the formal test
  categorise <- function(sA, sB, lA, lB, nA, nB) {
    upA <- lA > 0; upB <- lB > 0
    ifelse(sA & sB & upA == upB, ifelse(upA, "shared_up", "shared_down"),
           ifelse(sA & sB, "opposite_direction",
                  ifelse(sA, ifelse(upA, paste0(nA, "_only_up"), paste0(nA, "_only_down")),
                         ifelse(sB, ifelse(upB, paste0(nB, "_only_up"), paste0(nB, "_only_down")), "not_significant"))))
  }
  verdict <- function(category, supported, nA, nB, test) {
    v   <- rep(NA_character_, length(category))
    has <- !is.na(category)
    ct  <- category[has]
    sup <- isT(supported[has])
    v[has] <- ifelse(ct == "not_significant", "not significant",
                     ifelse(ct == "opposite_direction", "opposite direction",
                            ifelse(grepl("^shared", ct),
                                   ifelse(sup, paste0("shared, effect size differs (", test, ")"), "shared"),
                                   ifelse(grepl(paste0("^", nA, "_only"), ct),
                                          ifelse(sup, paste0(nA, "-specific (", test, " supported)"), paste0(nA, " only (", test, " n.s.)")),
                                          ifelse(sup, paste0(nB, "-specific (", test, " supported)"), paste0(nB, " only (", test, " n.s.)"))))))
    v
  }
  both <- mt$tested_Paryani & mt$tested_GSE79666
  
  # (1) striatum vs cingulate INSIDE Paryani: same donors, lab, batches; test = interaction E_P
  mt$HD_category_striatum_vs_cingulate <- categorise(mt$sig_A, mt$sig_B_P, mt$log2FC_A, mt$log2FC_B_P,
                                                     "striatum", "cingulate")
  mt$HD_category_striatum_vs_cingulate[!mt$tested_Paryani] <- NA
  mt$tissue_verdict_striatum_vs_cingulate <- verdict(mt$HD_category_striatum_vs_cingulate, mt$sig_E_P,
                                                     "striatum", "cingulate", "interaction")
  # (2) striatum (Paryani) vs BA4 (GSE79666): region AND study differ; test = difference z-test
  mt$HD_category_striatum_vs_BA4 <- categorise(mt$sig_A, mt$sig_B_T, mt$log2FC_A, mt$log2FC_B_T,
                                               "striatum", "BA4")
  mt$HD_category_striatum_vs_BA4[!both] <- NA
  mt$tissue_verdict_striatum_vs_BA4 <- verdict(mt$HD_category_striatum_vs_BA4, mt$striatum_BA4_differ,
                                               "striatum", "BA4", "difference test")
  # (3) all three regions (genes tested in both studies)
  lab <- cbind(ifelse(mt$sig_A, "striatum", NA), ifelse(mt$sig_B_P, "cingulate", NA), ifelse(mt$sig_B_T, "BA4", NA))
  mt$HD_pattern_3regions <- apply(lab, 1, function(z) if (all(is.na(z))) "none" else paste(z[!is.na(z)], collapse = "+"))
  mt$HD_pattern_3regions[!both] <- NA
  n_up   <- isT(mt$sig_A & mt$log2FC_A > 0) + isT(mt$sig_B_P & mt$log2FC_B_P > 0) + isT(mt$sig_B_T & mt$log2FC_B_T > 0)
  n_down <- isT(mt$sig_A & mt$log2FC_A < 0) + isT(mt$sig_B_P & mt$log2FC_B_P < 0) + isT(mt$sig_B_T & mt$log2FC_B_T < 0)
  mt$HD_pattern_direction <- ifelse(!both, NA,
                                    ifelse(n_up + n_down == 0, "none",
                                           ifelse(n_down == 0, "up", ifelse(n_up == 0, "down", "mixed"))))
  mt$HD_common_all_three <- both & isT(mt$HD_pattern_3regions == "striatum+cingulate+BA4") &
    mt$HD_pattern_direction %in% c("up", "down")
  # (4) does BA4 reproduce the cingulate result?
  sigBP_both <- mt$sig_B_P & both
  mt$BA4_same_direction       <- ifelse(sigBP_both & mt$reliable_B_T, sign(mt$log2FC_B_T) == sign(mt$log2FC_B_P), NA)
  mt$BA4_replicates_cingulate <- ifelse(sigBP_both, mt$sig_B_T & isT(sign(mt$log2FC_B_T) == sign(mt$log2FC_B_P)), NA)
  
  mt <- mt[order(pmin(mt$padj_A, mt$padj_B_P, mt$padj_B_T, na.rm = TRUE), na.last = TRUE), ]
  rownames(mt) <- NULL
  
  compact_cols <- c("gene_id", "gene_symbol", if ("gene_type" %in% names(mt)) "gene_type",
                    "tested_Paryani", "tested_GSE79666",
                    "HD_pattern_3regions", "HD_pattern_direction", "HD_common_all_three",
                    "HD_category_striatum_vs_cingulate", "tissue_verdict_striatum_vs_cingulate",
                    "HD_category_striatum_vs_BA4", "tissue_verdict_striatum_vs_BA4",
                    "log2FC_shr_A", "padj_A", "log2FC_shr_B_P", "padj_B_P", "log2FC_shr_B_T", "padj_B_T",
                    "log2FC_shr_E_P", "padj_E_P",
                    "BA4_same_direction", "BA4_replicates_cingulate",
                    "cortex_datasets_differ", "diff_padj_BA4_vs_cingulate",
                    "striatum_BA4_differ", "diff_padj_striatum_vs_BA4",
                    "baseMean_Paryani", "baseMean_GSE79666")
  compact <- mt[, compact_cols]
  write.csv(mt, out_path("tables", "MASTER_DESeq2_by_study_full.csv"), row.names = FALSE)
  write.csv(compact, out_path("tables", "MASTER_DESeq2_by_study_compact.csv"), row.names = FALSE)
  saveRDS(mt, out_path("objects", "MASTER_DESeq2_by_study.rds"))
  say("Master table: ", nrow(mt), " genes (tested in Paryani: ", sum(mt$tested_Paryani),
      ", in GSE79666: ", sum(mt$tested_GSE79666), ", in both: ", sum(both), ") x ", ncol(mt),
      " columns; compact version: ", ncol(compact), " columns")
})


# ---- 7. Summaries, gene lists, figures ----------------------------------------
run_step("7. Summaries, gene lists, figures", {
  summ <- do.call(rbind, lapply(names(specs), function(k) {
    r <- res_list[[k]]
    data.frame(contrast = k, model = specs[[k]]$model, description = specs[[k]]$label,
               genes_in_model = nrow(r),
               genes_tested = sum(!is.na(r$padj)),
               flagged_not_estimable = sum(!r$reliable),
               significant = sum(r$sig),
               up = sum(r$sig & isT(r$log2FC > 0)),
               down = sum(r$sig & isT(r$log2FC < 0)),
               significant_and_abs_shrunk_log2FC_gt1 = sum(r$sig & isT(abs(r$log2FC_shr) > 1)),
               median_lfcSE_reliable = round(median(r$lfcSE[r$reliable], na.rm = TRUE), 3),
               stringsAsFactors = FALSE)
  }))
  show_tab("DEG summary (padj < 0.05, non-estimable genes excluded)", summ)
  write.csv(summ, out_path("tables", "SUMMARY_DEG_counts_per_contrast.csv"), row.names = FALSE)
  
  both <- mt$tested_Paryani & mt$tested_GSE79666
  lev2 <- function(nA, nB) c("shared_up", "shared_down", "opposite_direction",
                             paste0(nA, "_only_up"), paste0(nA, "_only_down"),
                             paste0(nB, "_only_up"), paste0(nB, "_only_down"), "not_significant")
  cat_table <- function(category, supported, levels_cat, test) {
    s_lev <- c(paste0(test, "_padj_lt_0.05"), paste0(test, "_ns"))
    ct <- table(factor(category, levels = levels_cat),
                factor(ifelse(isT(supported), s_lev[1], s_lev[2]), levels = s_lev))
    data.frame(category = rownames(ct), as.data.frame.matrix(ct), total = as.vector(rowSums(ct)),
               check.names = FALSE)
  }
  tab_SC <- cat_table(mt$HD_category_striatum_vs_cingulate[mt$tested_Paryani], mt$sig_E_P[mt$tested_Paryani],
                      lev2("striatum", "cingulate"), "interaction")
  tab_SB <- cat_table(mt$HD_category_striatum_vs_BA4[both], mt$striatum_BA4_differ[both],
                      lev2("striatum", "BA4"), "difference")
  show_tab("HD categories: striatum (A) vs cingulate (B_P) within Paryani x interaction E_P", tab_SC)
  show_tab("HD categories: striatum (A, Paryani) vs BA4 (B_T, GSE79666) x difference test", tab_SB)
  write.csv(tab_SC, out_path("tables", "SUMMARY_categories_striatum_vs_cingulate.csv"), row.names = FALSE)
  write.csv(tab_SB, out_path("tables", "SUMMARY_categories_striatum_vs_BA4.csv"), row.names = FALSE)
  
  pat_levels <- c("striatum+cingulate+BA4", "striatum+cingulate", "striatum+BA4", "cingulate+BA4",
                  "striatum", "cingulate", "BA4", "none")
  pt <- table(factor(mt$HD_pattern_3regions[both], levels = pat_levels),
              factor(mt$HD_pattern_direction[both], levels = c("up", "down", "mixed", "none")))
  pat_tab <- data.frame(pattern = rownames(pt), as.data.frame.matrix(pt), total = as.vector(rowSums(pt)),
                        check.names = FALSE)
  show_tab("Significant HD effect in which regions (genes tested in both studies)", pat_tab)
  write.csv(pat_tab, out_path("tables", "SUMMARY_HD_pattern_3regions.csv"), row.names = FALSE)
  
  rel_PT <- both & mt$reliable_B_P & mt$reliable_B_T
  rel_AP <- mt$reliable_A & mt$reliable_B_P
  rel_AT <- both & mt$reliable_A & mt$reliable_B_T
  sigP   <- mt$sig_B_P & both
  agree_tab <- data.frame(
    measure = c("Genes tested in both studies",
                "Spearman rho, cingulate vs BA4 HD effect (shrunken log2FC, all reliable genes)",
                "Spearman rho, cingulate vs BA4 HD effect (genes significant in cingulate)",
                "Genes significant in cingulate (B_P) and tested in BA4",
                "  of these: same direction in BA4 (%)",
                "  of these: significant in BA4 with the same direction",
                "Genes significant in BA4 (B_T)",
                "Genes with a significant BA4 vs cingulate difference",
                "Spearman rho, striatum (A) vs cingulate (B_P) HD effect",
                "Spearman rho, striatum (A) vs BA4 (B_T) HD effect",
                "Genes with a significant striatum vs BA4 difference",
                "Genes significant in all three regions with the same direction"),
    value = c(sum(both),
              scor(mt$log2FC_shr_B_P[rel_PT], mt$log2FC_shr_B_T[rel_PT]),
              scor(mt$log2FC_shr_B_P[rel_PT & sigP], mt$log2FC_shr_B_T[rel_PT & sigP]),
              sum(sigP),
              round(100 * mean(mt$BA4_same_direction[sigP], na.rm = TRUE), 1),
              sum(mt$BA4_replicates_cingulate[sigP], na.rm = TRUE),
              sum(mt$sig_B_T),
              sum(mt$cortex_datasets_differ),
              scor(mt$log2FC_shr_A[rel_AP], mt$log2FC_shr_B_P[rel_AP]),
              scor(mt$log2FC_shr_A[rel_AT], mt$log2FC_shr_B_T[rel_AT]),
              sum(mt$striatum_BA4_differ),
              sum(mt$HD_common_all_three)),
    stringsAsFactors = FALSE)
  show_tab("Agreement between regions and datasets", agree_tab)
  write.csv(agree_tab, out_path("tables", "SUMMARY_agreement_regions_datasets.csv"), row.names = FALSE)
  
  # Gene lists for enrichment (unversioned Ensembl IDs); universe = genes with a padj in that contrast
  write_list <- function(name, ids) writeLines(strip_version(ids), out_path("gene_lists", paste0(name, ".txt")))
  for (k in c("A", "B_P", "B_T", "E_P")) {
    write_list(paste0("universe_", k), mt$gene_id[!is.na(mt[[paste0("padj_", k)]])])
    write_list(paste0("contrast_", k, "_up"),   mt$gene_id[mt[[paste0("sig_", k)]] & isT(mt[[paste0("log2FC_", k)]] > 0)])
    write_list(paste0("contrast_", k, "_down"), mt$gene_id[mt[[paste0("sig_", k)]] & isT(mt[[paste0("log2FC_", k)]] < 0)])
  }
  write_list("universe_both_studies", mt$gene_id[both])
  for (d in c("up", "down")) {
    write_list(paste0("striatum_specific_vs_cingulate_", d, "_interaction_supported"),
               mt$gene_id[isT(mt$HD_category_striatum_vs_cingulate == paste0("striatum_only_", d)) & mt$sig_E_P])
    write_list(paste0("cingulate_specific_vs_striatum_", d, "_interaction_supported"),
               mt$gene_id[isT(mt$HD_category_striatum_vs_cingulate == paste0("cingulate_only_", d)) & mt$sig_E_P])
    write_list(paste0("shared_striatum_cingulate_", d),
               mt$gene_id[isT(mt$HD_category_striatum_vs_cingulate == paste0("shared_", d))])
    write_list(paste0("all_three_regions_", d),
               mt$gene_id[mt$HD_common_all_three & isT(mt$HD_pattern_direction == d)])
    write_list(paste0("cingulate_replicated_in_BA4_", d),
               mt$gene_id[isT(mt$BA4_replicates_cingulate) & isT(if (d == "up") mt$log2FC_B_P > 0 else mt$log2FC_B_P < 0)])
  }
  for (pp in setdiff(pat_levels, "none")) {
    for (d in c("up", "down", "mixed")) {
      sel <- both & isT(mt$HD_pattern_3regions == pp) & isT(mt$HD_pattern_direction == d)
      if (any(sel)) write_list(paste0("pattern3_", gsub("+", "_and_", pp, fixed = TRUE), "_", d), mt$gene_id[sel])
    }
  }
  
  # Figures: HD effect of two regions against each other
  scatter_cat <- function(x, y, category, nA, nB, xlab, ylab, title, file) {
    lab <- ifelse(grepl("^shared", category), "shared",
                  ifelse(grepl(paste0("^", nA, "_only"), category), paste(nA, "only"),
                         ifelse(grepl(paste0("^", nB, "_only"), category), paste(nB, "only"),
                                ifelse(category == "opposite_direction", "opposite direction", "not significant"))))
    lev <- c("not significant", paste(nA, "only"), paste(nB, "only"), "shared", "opposite direction")
    pd  <- data.frame(x = x, y = y, cat = factor(lab, levels = lev))
    pd  <- pd[complete.cases(pd), ]
    pd  <- pd[order(pd$cat), ]
    p <- ggplot(pd, aes(x, y, colour = cat)) +
      geom_hline(yintercept = 0, colour = "#c3c2b7") + geom_vline(xintercept = 0, colour = "#c3c2b7") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#898781") +
      geom_point(size = 0.7, alpha = 0.7) +
      scale_colour_manual(values = setNames(c("#d9d8d2", "#eb6834", "#2a78d6", "#1baf7a", "#52514e"), lev),
                          drop = FALSE) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(x = xlab, y = ylab, colour = "padj < 0.05 in", title = title,
           subtitle = paste0("Spearman rho = ", scor(pd$x, pd$y), "; dashed line = equal effect")) +
      theme_thesis
    save_png(p, file)
  }
  scatter_cat(mt$log2FC_shr_A[rel_AP], mt$log2FC_shr_B_P[rel_AP], mt$HD_category_striatum_vs_cingulate[rel_AP],
              "striatum", "cingulate",
              "Striatum: HD vs CTRL (shrunken log2FC)", "Cingulate cortex: HD vs CTRL (shrunken log2FC)",
              "HD effect: striatum vs cingulate cortex (Paryani, same donors)",
              "RESULT_HD_striatum_vs_cingulate.png")
  scatter_cat(mt$log2FC_shr_A[rel_AT], mt$log2FC_shr_B_T[rel_AT], mt$HD_category_striatum_vs_BA4[rel_AT],
              "striatum", "BA4",
              "Striatum (Paryani): HD vs CTRL (shrunken log2FC)", "Motor cortex BA4 (GSE79666): HD vs CTRL",
              "HD effect: striatum vs motor cortex BA4 (different studies)",
              "RESULT_HD_striatum_vs_BA4.png")
  
  pd2 <- mt[rel_PT, c("log2FC_shr_B_P", "log2FC_shr_B_T", "log2FC_B_P", "log2FC_B_T", "sig_B_P", "sig_B_T")]
  same <- sign(pd2$log2FC_B_P) == sign(pd2$log2FC_B_T)
  lev3 <- c("not significant", "significant in one dataset", "significant in both, same direction",
            "significant in both, opposite direction")
  pd2$col <- factor(ifelse(pd2$sig_B_P & pd2$sig_B_T & same, lev3[3],
                           ifelse(pd2$sig_B_P & pd2$sig_B_T, lev3[4],
                                  ifelse(pd2$sig_B_P | pd2$sig_B_T, lev3[2], lev3[1]))), levels = lev3)
  pd2 <- pd2[complete.cases(pd2[, c("log2FC_shr_B_P", "log2FC_shr_B_T")]), ]
  pd2 <- pd2[order(pd2$col), ]
  p2 <- ggplot(pd2, aes(log2FC_shr_B_P, log2FC_shr_B_T, colour = col)) +
    geom_hline(yintercept = 0, colour = "#c3c2b7") + geom_vline(xintercept = 0, colour = "#c3c2b7") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#898781") +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = setNames(c("#d9d8d2", "#2a78d6", "#1baf7a", "#52514e"), lev3), drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "Cingulate cortex (Paryani): HD vs CTRL, shrunken log2FC",
         y = "Motor cortex BA4 (GSE79666): HD vs CTRL, shrunken log2FC", colour = NULL,
         title = "HD effect in the two cortical datasets",
         subtitle = paste0("Spearman rho = ", scor(pd2$log2FC_shr_B_P, pd2$log2FC_shr_B_T),
                           "; genes with a significant dataset difference: ", sum(mt$cortex_datasets_differ))) +
    theme_thesis
  save_png(p2, "RESULT_HD_cingulate_vs_BA4.png")
})


# ---- 8. Excel workbook, session info, end ---------------------------------------
run_step("8. Excel workbook and session info", {
  writeLines(capture.output(sessionInfo()), out_path("objects", "sessionInfo.txt"))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    readme <- data.frame(
      item = c(names(specs),
               "log2FC_<k>", "log2FC_shr_<k>", "lfcSE_<k>, pvalue_<k>, padj_<k>", "reliable_<k>", "sig_<k>",
               "tested_Paryani / tested_GSE79666",
               "HD_pattern_3regions, HD_pattern_direction, HD_common_all_three",
               "HD_category_striatum_vs_cingulate, tissue_verdict_striatum_vs_cingulate",
               "HD_category_striatum_vs_BA4, tissue_verdict_striatum_vs_BA4",
               "BA4_same_direction, BA4_replicates_cingulate",
               "diff*_BA4_vs_cingulate, cortex_datasets_differ",
               "diff*_striatum_vs_BA4, striatum_BA4_differ",
               "Models"),
      description = c(unname(vapply(names(specs), function(k) specs[[k]]$label, character(1))),
                      "unshrunken DESeq2 log2 fold change (use for statistics)",
                      "ashr-shrunken log2 fold change, estimable genes only (use for ranking, plots, effect-size cut-offs)",
                      "standard error, Wald p-value, Benjamini-Hochberg adjusted p-value",
                      "FALSE if a compared group had zero counts in every sample or |log2FC| > 10 (not estimable)",
                      "padj < 0.05 and reliable",
                      "gene passed the expression filter of that study (otherwise its columns for that study are empty)",
                      "which of striatum (A), cingulate (B_P), BA4 (B_T) are significant; direction; TRUE = all three, same direction",
                      "striatum vs cingulate within Paryani (same donors, lab, batches); tissue-specific only if the interaction E_P is significant",
                      "striatum (Paryani) vs BA4 (GSE79666); difference test of independent estimates; region and study are confounded",
                      "for genes significant in cingulate: same sign in BA4 / significant in BA4 with the same sign",
                      "BA4 minus cingulate log2FC, z, p, BH padj; TRUE = significant difference (region AND study differ)",
                      "striatum minus BA4 log2FC, z, p, BH padj; TRUE = significant difference (region AND study differ)",
                      paste0("Model P (Paryani): ", paste(deparse(design_P), collapse = ""),
                             "; Model L (GSE79666): ", paste(deparse(design_L), collapse = ""),
                             "; tissue corrected for 7 samples; excluded from the models: ",
                             if (length(exclude_from_models)) paste(exclude_from_models, collapse = ", ") else "none",
                             " (sheet region_corrections)")),
      stringsAsFactors = FALSE)
    wb <- openxlsx::createWorkbook()
    sheets <- list(README = readme, models = design_tab, DEG_summary = summ,
                   dispersion = disp_tab, variability = var_tab, donors_Paryani = rep_tab,
                   striatum_vs_cingulate = tab_SC, striatum_vs_BA4 = tab_SB, pattern_3regions = pat_tab,
                   agreement = agree_tab, master_compact = compact, region_corrections = fix)
    for (nm in names(sheets)) {
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, sheets[[nm]])
      openxlsx::freezePane(wb, nm, firstRow = TRUE)
      openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = 18)
    }
    openxlsx::saveWorkbook(wb, out_path("tables", "MASTER_DESeq2_by_study.xlsx"), overwrite = FALSE)
    say("Excel workbook: tables/MASTER_DESeq2_by_study.xlsx")
  } else {
    say("openxlsx not installed - Excel workbook skipped (all CSV files are complete). ",
        "Optional: install.packages('openxlsx')")
  }
  say("\nFINISHED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | R memory in use: ", mem_mb(), " Mb")
  say("All outputs: ", run_dir)
  say("Please send: run_log.txt + figures QC_PCA_Paryani_uncorrected.png, QC_PCA_GSE79666_uncorrected.png, ",
      "QC_dispersion_comparison.png, RESULT_HD_striatum_vs_cingulate.png, RESULT_HD_cingulate_vs_BA4.png")
})