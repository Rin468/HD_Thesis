# =============================================================================
#  07_thesis_requests_20261017.R
#  RECOVERED on 2026-10-31 from HD_Thesis/.Rhistory (lines 42-384), where it was
#  run from the console on 17 Oct 2026 but never saved as a file.
#  Produced: thesis_inputs_20261017/, thesis_figures/, final_runs/GO_20261017_223957/,
#            Salmon_run_summary_Paryani.csv, send_to_claude_20261017/
#  Code is verbatim; only indentation was lost in the history.
#  NOTE: REQUEST 3 creates a NEW time-stamped GO folder each time it is run.
# =============================================================================

# =============================================================================
#  thesis_requests_20261017.R
#  Five data/figure requests that complete the thesis draft (ranked).
#  Run the SETUP block first, then any request block on its own.
#  Nothing here overwrites existing results: every request writes to a NEW folder.
#
#    REQUEST 1  (essential)     Collect the final-run tables, figures and scripts the thesis cites
#    REQUEST 2  (essential)     Thesis figures with plain-language labels (volcano plots,
#                               DEG counts, robustness)
#    REQUEST 3  (essential)     GO biological-process enrichment on the FINAL DEG lists
#    REQUEST 4  (essential)     Pre-processing record for Methods (Salmon versions, settings,
#                               mapping rates; trimming tools)
#    REQUEST 5  (nice-to-have)  Cell-type marker panel (helps interpret neuronal loss vs regulation)
#
#  Save as HD_Thesis/scripts/Merged_analysis/07_thesis_requests.R
# =============================================================================
# ---- SETUP (run first) --------------------------------------------------------
base_dir <- "C:/!Virag.ORI/Bioinfo/szakdoga"
ma       <- file.path(base_dir, "Merged_analysis")
main_run <- file.path(ma, "final_runs", "bystudy_20261015_151051")      # analysis of record
sens_run <- file.path(ma, "final_runs", "sensitivity_20261017_163558")  # final sensitivity run
s9_run   <- file.path(ma, "final_runs", "S9followup_20261017_171514")   # S9 follow-up
stopifnot(dir.exists(main_run), dir.exists(sens_run), dir.exists(s9_run))
res_main  <- readRDS(file.path(main_run, "objects", "results_by_contrast.rds"))
ann_files <- list.files(file.path(ma, "results", "objects"),
pattern = "^gene_annotation_GENCODEv44(_v[0-9]+)?\\.rds$", full.names = TRUE)
ann       <- as.data.frame(readRDS(ann_files[which.max(file.info(ann_files)$mtime)]))
stopifnot(all(c("A", "B_P", "B_T", "E_P") %in% names(res_main)),
all(c("gene_id", "gene_symbol") %in% names(ann)))
# Plain-language names used in every figure and table (no bare codes)
comparisons <- c(A   = "Striatum (Paryani et al.)",
B_P = "Cingulate cortex (Paryani et al.)",
B_T = "Motor cortex BA4 (Lin et al.)")
cat("Columns stored per comparison:", paste(names(res_main$A), collapse = ", "), "\n")
# ==== REQUEST 1 (essential): collect the final-run outputs the thesis cites ====
# Why: every number in the thesis must trace to a file. The final-run numbers are
# currently only in the project notes; these files confirm them and provide the
# QC and comparison figures.
out_dir <- file.path(ma, "thesis_inputs_20261017")
dir.create(out_dir, showWarnings = FALSE)
wanted <- c(
file.path(main_run, "run_log.txt"),
file.path(main_run, "objects", "sessionInfo.txt"),
file.path(main_run, "tables", c(
"SUMMARY_DEG_counts_per_contrast.csv", "MODELS_design_summary.csv",
"colData_final.csv", "region_corrections.csv", "QC_donor_structure_Paryani.csv",
"QC_dispersion_comparison.csv", "QC_PCA_R2_Paryani_uncorrected.csv",
"QC_PCA_R2_GSE79666_uncorrected.csv", "SUMMARY_agreement_regions_datasets.csv",
"SUMMARY_categories_striatum_vs_cingulate.csv", "SUMMARY_categories_striatum_vs_BA4.csv",
"SUMMARY_HD_pattern_3regions.csv", "MASTER_DESeq2_by_study_compact.csv")),
file.path(main_run, "figures", c(
"QC_PCA_Paryani_uncorrected_by_batch_study.png", "QC_PCA_GSE79666_uncorrected.png",
"QC_tissue_markers_corrected.png", "QC_dispersion_comparison.png",
"RESULT_HD_striatum_vs_cingulate.png", "RESULT_HD_cingulate_vs_BA4.png",
"RESULT_HD_striatum_vs_BA4.png")),
file.path(main_run, "gene_lists", c("all_three_regions_up.txt", "all_three_regions_down.txt")),
file.path(sens_run, "tables", c("SENS_comparison_vs_main.csv", "SENS_core_genes_check.csv")),
file.path(ma, "figures", c("stage1_sex_inference_XIST_vs_Ymarkers_v3.png",
"stage3b_region_marker_check_v2.png")),
file.path(ma, "results", "tables", "sex_markers_log2CPM_inferred_sex_v3.csv"),
file.path(ma, "results", "tables", "QC", "region_marker_scores_per_sample_v2.csv"),
list.files(file.path(base_dir, "Git", "HD_Thesis", "scripts", "Merged_analysis"),
pattern = "^0[1-6].*\\.R$", full.names = TRUE))
found  <- file.exists(wanted)
copied <- file.copy(wanted[found], out_dir, overwrite = FALSE)
print(data.frame(file = basename(wanted), found = found))
cat("Copied", sum(copied), "of", length(wanted), "files to", out_dir, "\n")
cat("Size of MASTER_DESeq2_by_study_compact.csv (MB):",
round(file.size(file.path(main_run, "tables", "MASTER_DESeq2_by_study_compact.csv")) / 1e6, 1), "\n")
# Save the adjusted-p distribution you printed in the console, so it is file-backed
padj_dist <- do.call(rbind, lapply(names(comparisons), function(k) {
d <- res_main[[k]][res_main[[k]]$sig %in% TRUE, ]
data.frame(comparison = comparisons[[k]], DEGs = nrow(d),
padj_up_to_0.001     = sum(d$padj <= 0.001),
padj_0.001_to_0.01   = sum(d$padj > 0.001 & d$padj <= 0.01),
padj_0.01_to_0.05    = sum(d$padj > 0.01 & d$padj <= 0.05),
median_abs_log2FC    = round(median(abs(d$log2FC)), 2))
}))
print(padj_dist)
write.csv(padj_dist, file.path(out_dir, "DEG_padj_distribution_by_tissue.csv"), row.names = FALSE)
# THEN: in File Explorer right-click the folder thesis_inputs_20261017 >
#       Send to > Compressed (zipped) folder, and attach the zip.
# ==== REQUEST 2 (essential): thesis figures with plain-language labels =========
# Why: the assignment's core comparison (HD vs control within each tissue) has no
# final-run figure yet, and the robustness figures use internal codes.
suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)
fig_dir <- file.path(ma, "thesis_figures")
dir.create(fig_dir, showWarnings = FALSE)
# Use shrunken fold changes for plotting if the main run stored them, else the MLE
lfc_col <- grep("shr", names(res_main$A), value = TRUE, ignore.case = TRUE)[1]
if (is.na(lfc_col)) lfc_col <- "log2FC"
cat("Volcano x-axis uses column:", lfc_col, "\n")
x_lab <- if (lfc_col == "log2FC") "log2 fold change (HD relative to control)" else
"Shrunken log2 fold change (HD relative to control)"
file_stub <- c(A = "striatum", B_P = "cingulate_cortex", B_T = "motor_cortex_BA4")
for (k in names(comparisons)) {
d <- as.data.frame(res_main[[k]])
d$gene_symbol <- ann$gene_symbol[match(rownames(d), ann$gene_id)]
d$x <- d[[lfc_col]]
d <- d[d$reliable %in% TRUE & !is.na(d$padj) & is.finite(d$x), ]
up <- (d$sig %in% TRUE) & ((d$log2FC > 0) %in% TRUE)
dn <- (d$sig %in% TRUE) & ((d$log2FC < 0) %in% TRUE)
lev <- c("Not significant",
paste0("Higher in HD (", format(sum(up), big.mark = ","), " genes)"),
paste0("Lower in HD (", format(sum(dn), big.mark = ","), " genes)"))
d$status <- factor(ifelse(up, lev[2], ifelse(dn, lev[3], lev[1])), levels = lev)
sigd <- d[d$sig %in% TRUE & !is.na(d$gene_symbol) & d$gene_symbol != "", ]
lab  <- head(sigd[order(sigd$padj), ], 10)            # 10 smallest adjusted p-values
d    <- d[order(d$status), ]                          # draw grey points first
p <- ggplot(d, aes(x, -log10(padj), colour = status)) +
geom_point(size = 0.8, alpha = 0.6) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
geom_text_repel(data = lab, aes(label = gene_symbol), colour = "black", size = 3,
max.overlaps = Inf, show.legend = FALSE) +
scale_colour_manual(values = setNames(c("grey80", "#D55E00", "#0072B2"), lev)) +
labs(title = paste0(comparisons[[k]], ": Huntington's disease vs control"),
x = x_lab, y = "-log10 adjusted p-value", colour = NULL) +
theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, paste0("Volcano_", file_stub[[k]], ".png")), p, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, paste0("Volcano_", file_stub[[k]], ".pdf")), p, width = 7, height = 5.5)
write.csv(lab[, c("gene_symbol", "log2FC", "padj")],
file.path(fig_dir, paste0("Volcano_labelled_genes_", file_stub[[k]], ".csv")), row.names = FALSE)
}
# DEG counts per tissue
cnt <- do.call(rbind, lapply(names(comparisons), function(k) {
d <- res_main[[k]]; s <- d$sig %in% TRUE
data.frame(tissue = comparisons[[k]], direction = c("Higher in HD", "Lower in HD"),
genes = c(sum(s & d$log2FC > 0, na.rm = TRUE), sum(s & d$log2FC < 0, na.rm = TRUE)))
}))
cnt$tissue <- factor(cnt$tissue, levels = comparisons)
p_cnt <- ggplot(cnt, aes(tissue, genes, fill = direction)) +
geom_col(position = position_dodge(width = 0.75), width = 0.7) +
geom_text(aes(label = format(genes, big.mark = ",")), position = position_dodge(width = 0.75),
vjust = -0.3, size = 3.2) +
scale_fill_manual(values = c("Higher in HD" = "#D55E00", "Lower in HD" = "#0072B2")) +
labs(x = NULL, y = "Differentially expressed genes (adjusted p < 0.05)", fill = NULL) +
theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "DEG_counts_by_tissue.png"), p_cnt, width = 7, height = 4.5, dpi = 300)
write.csv(cnt, file.path(fig_dir, "DEG_counts_by_tissue.csv"), row.names = FALSE)
# Robustness: share of main DEGs retained + rank correlation of fold changes
comp <- read.csv(file.path(sens_run, "tables", "SENS_comparison_vs_main.csv"))
s9f  <- read.csv(file.path(s9_run, "tables", "S9followup_comparison_vs_main.csv"))
comp <- rbind(comp, s9f[, names(comp)])
lab_an <- c(voom_noblock = "Different method (limma-voom)",
voom_donor = "Different method, accounting for repeated donors",
S1_without_batch4 = "Without the HD-only sequencing batch",
S2_without_corrected_flagged = "Without relabelled or flagged samples",
S3_without_C5596_Accum = "Without the possible outlier control",
S4_striatum_only = "Striatum samples modelled alone",
S5_cortex_only = "Cortex samples modelled alone",
S6_without_juvenile_onset = "Without juvenile-onset cases",
S7_GSE79666_without_CTRL3 = "Motor cortex without the only female control",
S8_mitochondrial_covariate = "Adjusted for mitochondrial read share",
S9_without_sex_corrected = "Without the four sex-corrected samples",
S9a_without_H5608_Accum = "Without one donor's second cortex sample",
S9b_without_C98_Cing = "Without one sex-corrected cortex control",
S9c_without_H5608_both_cortex = "Without both cortex samples of one donor")
lab_ct <- c(A = "Striatum: HD vs control", B_P = "Cingulate cortex: HD vs control",
B_T = "Motor cortex BA4: HD vs control",
E_P = "Difference between striatal and cingulate HD effects")
comp <- comp[comp$contrast %in% names(lab_ct) & comp$analysis %in% names(lab_an), ]
comp$analysis_label   <- factor(lab_an[comp$analysis], levels = rev(unname(lab_an)))
comp$comparison_label <- factor(lab_ct[comp$contrast], levels = unname(lab_ct))
p_rob <- ggplot(comp, aes(retained_pct, analysis_label)) +
geom_col(fill = "#5B8DB8", width = 0.7) +
geom_text(aes(label = paste0("\u03c1 = ", sprintf("%.2f", rho_log2FC))), hjust = -0.1, size = 2.7) +
facet_wrap(~ comparison_label, nrow = 1) +
scale_x_continuous(limits = c(0, 130), breaks = c(0, 25, 50, 75, 100)) +
labs(x = "Main-analysis DEGs significant again with the same direction (%)", y = NULL) +
theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "Robustness_of_DEG_lists.png"), p_rob, width = 13, height = 6, dpi = 300)
write.csv(comp, file.path(fig_dir, "Robustness_of_DEG_lists.csv"), row.names = FALSE)
cat("Figures written to", fig_dir, "\n")
# ==== REQUEST 3 (essential): GO biological-process enrichment, FINAL lists ======
# Why: the aims promise functional characterisation, but the existing GO tables come
# from the retired per-region analysis (GEO labels). Universe = genes actually tested.
for (p in c("clusterProfiler", "org.Hs.eg.db")) {
if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2) })
go_dir <- file.path(ma, "final_runs", format(Sys.time(), "GO_%Y%m%d_%H%M%S"))
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
strip <- function(x) sub("\\..*$", "", x)               # remove Ensembl version suffix
gene_lists <- list(); universes <- list()
for (k in names(comparisons)) {
d <- res_main[[k]]
universes[[k]] <- unique(strip(rownames(d)[d$reliable %in% TRUE & !is.na(d$padj)]))
gene_lists[[paste0(k, "_up")]]   <- unique(strip(rownames(d)[d$sig %in% TRUE & d$log2FC > 0]))
gene_lists[[paste0(k, "_down")]] <- unique(strip(rownames(d)[d$sig %in% TRUE & d$log2FC < 0]))
}
d <- res_main$B_P                                       # cingulate check on the stronger subset
gene_lists$B_P_strong_up   <- unique(strip(rownames(d)[d$sig %in% TRUE & d$padj < 0.01 & d$log2FC > 0]))
gene_lists$B_P_strong_down <- unique(strip(rownames(d)[d$sig %in% TRUE & d$padj < 0.01 & d$log2FC < 0]))
universe_of <- function(nm) universes[[sub("_(strong_)?(up|down)$", "", nm)]]
plain <- c(A_up = "Striatum, higher in HD", A_down = "Striatum, lower in HD",
B_P_up = "Cingulate, higher in HD", B_P_down = "Cingulate, lower in HD",
B_T_up = "Motor cortex BA4, higher in HD", B_T_down = "Motor cortex BA4, lower in HD",
B_P_strong_up = "Cingulate (adjusted p < 0.01), higher in HD",
B_P_strong_down = "Cingulate (adjusted p < 0.01), lower in HD")
mapped <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
go_res <- list(); go_summary <- list()
for (nm in names(gene_lists)) {
g <- gene_lists[[nm]]; u <- universe_of(nm)
ego <- enrichGO(gene = g, universe = u, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
readable = TRUE)
n_terms <- if (is.null(ego)) 0 else nrow(as.data.frame(ego))
ego_s   <- if (n_terms > 0) clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust",
select_fun = min) else NULL
go_res[[nm]] <- list(full = ego, simplified = ego_s)
if (n_terms > 0) {
write.csv(as.data.frame(ego),   file.path(go_dir, paste0("GO_BP_", nm, "_all.csv")), row.names = FALSE)
write.csv(as.data.frame(ego_s), file.path(go_dir, paste0("GO_BP_", nm, "_simplified.csv")), row.names = FALSE)
}
go_summary[[nm]] <- data.frame(gene_list = plain[[nm]], DEGs = length(g),
DEGs_mapped_in_OrgDb = sum(g %in% mapped),
universe = length(u), universe_mapped_in_OrgDb = sum(u %in% mapped),
terms_FDR_0.05 = n_terms,
terms_after_simplify = if (n_terms > 0) nrow(as.data.frame(ego_s)) else 0)
cat(nm, ":", n_terms, "terms\n")
}
go_summary <- do.call(rbind, go_summary); rownames(go_summary) <- NULL
print(go_summary)
write.csv(go_summary, file.path(go_dir, "GO_BP_summary.csv"), row.names = FALSE)
saveRDS(go_res, file.path(go_dir, "GO_BP_results.rds"))
top <- do.call(rbind, lapply(names(go_res), function(nm) {
s <- go_res[[nm]]$simplified
if (is.null(s)) return(NULL)
s <- as.data.frame(s); s <- head(s[order(s$p.adjust), ], 5)
data.frame(gene_list = plain[[nm]], term = s$Description, padj = s$p.adjust, count = s$Count)
}))
if (!is.null(top) && nrow(top) > 0) {
top_main <- top[!grepl("adjusted p < 0.01", top$gene_list), ]
top_main$gene_list <- factor(top_main$gene_list, levels = unname(plain[1:6]))
p_go <- ggplot(top_main, aes(gene_list, term, size = count, colour = -log10(padj))) +
geom_point() + scale_colour_viridis_c(name = "-log10 adjusted p") +
labs(x = NULL, y = NULL, size = "Genes in term") +
theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(go_dir, "GO_BP_top_terms_by_tissue.png"), p_go, width = 10, height = 7, dpi = 300)
write.csv(top, file.path(go_dir, "GO_BP_top5_simplified_terms.csv"), row.names = FALSE)
}
writeLines(c(capture.output(sessionInfo()),
paste("clusterProfiler", packageVersion("clusterProfiler")),
paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db"))),
file.path(go_dir, "sessionInfo.txt"))
cat("GO outputs in", go_dir, "\n")
# ==== REQUEST 4 (essential): pre-processing record for Methods =================
# Why: Methods must state software versions, Salmon settings and mapping rates.
# Part A (R, laptop): Paryani Salmon runs. Salmon 2.x may name JSON fields differently,
# so the code prints the available fields first.
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(jsonlite)
get_field <- function(x, k) if (!is.null(x[[k]])) paste(unlist(x[[k]]), collapse = ",") else NA_character_
summarise_salmon <- function(quant_dir, srr_ids) {
mi_files <- list.files(quant_dir, pattern = "^meta_info\\.json$", recursive = TRUE, full.names = TRUE)
names(mi_files) <- sub("_quant$", "", basename(dirname(dirname(mi_files))))
cat("Fields in meta_info.json:", paste(names(fromJSON(mi_files[1])), collapse = ", "), "\n")
missing <- setdiff(srr_ids, names(mi_files))
if (length(missing)) cat("No meta_info.json for:", paste(missing, collapse = ", "), "\n")
do.call(rbind, lapply(intersect(srr_ids, names(mi_files)), function(s) {
mi <- fromJSON(mi_files[[s]])
ci_file <- file.path(dirname(dirname(mi_files[[s]])), "cmd_info.json")
ci <- if (file.exists(ci_file)) fromJSON(ci_file) else list()
data.frame(srr_id = s,
salmon_version  = get_field(mi, "salmon_version"),
library_type    = get_field(mi, "library_types"),
fragments_processed = get_field(mi, "num_processed"),
percent_mapped  = get_field(mi, "percent_mapped"),
index           = get_field(ci, "index"),
libType_option  = get_field(ci, "libType"),
gcBias          = !is.null(ci$gcBias),
seqBias         = !is.null(ci$seqBias),
stringsAsFactors = FALSE)
}))
}
bundle <- readRDS(file.path(base_dir, "Paryani", "RStudio", "GSE242195_analysis_bundle.rds"))
salmon_paryani <- summarise_salmon(file.path(base_dir, "Paryani", "quant"),
bundle$sample_metadata$srr_id)
salmon_paryani$batch <- bundle$sample_metadata$batch[match(salmon_paryani$srr_id,
bundle$sample_metadata$srr_id)]
dir.create(file.path(ma, "thesis_inputs_20261017"), showWarnings = FALSE)
write.csv(salmon_paryani, file.path(ma, "thesis_inputs_20261017", "Salmon_run_summary_Paryani.csv"),
row.names = FALSE)
print(table(salmon_paryani$salmon_version, useNA = "ifany"))
print(table(salmon_paryani$batch, salmon_paryani$library_type, useNA = "ifany"))
print(summary(as.numeric(salmon_paryani$percent_mapped)))
# On the Mac, for GSE79666 (same function, different paths):
#   tseng_meta <- readRDS("/Users/ronalab/Virag/RStudio/Tseng/Results/Tseng_metadata.rds")
#   salmon_lin <- summarise_salmon("/Users/ronalab/Virag/RStudio/Tseng/Outputs/salmon_quant",
#                                  tseng_meta$run_accession)
#   write.csv(salmon_lin, "Salmon_run_summary_Lin2016.csv", row.names = FALSE)
#
# Part B (bash, on the server): versions and the exact trimming commands used
#   fastp --version; fastqc --version; salmon --version
#   find /RoNAS/mvirag/arcticles -name "*.json" -path "*fastp*" | head -3
#   grep -m1 '"command"' $(find /RoNAS/mvirag/arcticles/Paryani_2024 -name "*.json" -path "*fastp*" | head -1)
#   grep -m1 '"command"' $(find /RoNAS/mvirag/arcticles/Tseng_2016 -name "*.json" -path "*fastp*" | head -1)
#   cat /RoNAS/mvirag/arcticles/Paryani_2024/Analysis/libtype_test_results.txt
# ==== REQUEST 5 (nice-to-have): cell-type marker panel =========================
# Why: in bulk tissue, a lower signal for neuronal genes can mean fewer neurons
# rather than less expression per cell. Marker behaviour shows how strongly this
# applies in each tissue (descriptive only; it cannot separate the two).
markers <- c(PPP1R1B = "Medium spiny neuron", PDE10A = "Medium spiny neuron",
DRD1 = "Medium spiny neuron", DRD2 = "Medium spiny neuron", PENK = "Medium spiny neuron",
ADORA2A = "Medium spiny neuron", GPR88 = "Medium spiny neuron",
SLC17A7 = "Cortical excitatory neuron", RBFOX3 = "Neuron (general)",
SNAP25 = "Neuron (general)", GFAP = "Astrocyte", AQP4 = "Astrocyte",
AIF1 = "Microglia", MBP = "Oligodendrocyte", PLP1 = "Oligodendrocyte")
mk <- do.call(rbind, lapply(names(comparisons), function(k) {
d <- res_main[[k]]
i <- match(names(markers), ann$gene_symbol[match(rownames(d), ann$gene_id)])
data.frame(tissue = comparisons[[k]], gene = names(markers), cell_type = unname(markers),
log2FC = d$log2FC[i], lfcSE = d$lfcSE[i], padj = d$padj[i])
}))
mk$tissue <- factor(mk$tissue, levels = comparisons)
mk$gene   <- factor(mk$gene, levels = rev(names(markers)))
dir.create(file.path(ma, "thesis_figures"), showWarnings = FALSE)
p_mk <- ggplot(mk, aes(log2FC, gene, colour = cell_type)) +
geom_vline(xintercept = 0, colour = "grey50") +
geom_errorbar(aes(xmin = log2FC - 1.96 * lfcSE, xmax = log2FC + 1.96 * lfcSE),
width = 0.3, orientation = "y") +
geom_point(size = 2) + facet_wrap(~ tissue, nrow = 1) +
labs(x = "log2 fold change, HD relative to control (bars: 95% interval)", y = NULL,
colour = "Marker for") +
theme_bw(base_size = 10)
ggsave(file.path(ma, "thesis_figures", "Cell_type_markers_by_tissue.png"), p_mk, width = 11, height = 5, dpi = 300)
write.csv(mk, file.path(ma, "thesis_figures", "Cell_type_markers_by_tissue.csv"), row.names = FALSE)
print(mk[order(mk$tissue, mk$cell_type), ])
ma   <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
send <- file.path(ma, "send_to_claude_20261017")
dir.create(send, showWarnings = FALSE)
src <- c(file.path(ma, "thesis_inputs_20261017"),
file.path(ma, "thesis_figures"),
file.path(ma, "final_runs", "GO_20261017_223957"))
for (s in src) {
f <- list.files(s, full.names = TRUE)
f <- f[!grepl("\\.(rds|pdf)$", f, ignore.case = TRUE)]
dir.create(file.path(send, basename(s)), showWarnings = FALSE)
file.copy(f, file.path(send, basename(s)), overwrite = TRUE)
}
round(sum(file.size(list.files(send, recursive = TRUE, full.names = TRUE))) / 1e6, 1)
