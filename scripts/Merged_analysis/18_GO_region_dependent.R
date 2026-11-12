# =============================================================================
#  18_GO_region_dependent.R   (written 1 Nov 2026)
#  GO biological-process over-representation analysis of the genes whose HD effect
#  differs between the striatum and the pooled cortex (cingulate + motor cortex):
#    - all 1,867 region-dependent genes, by direction (1,316 more negative / 551 more
#      positive HD effect in the striatum)
#    - the 966 robust genes (high-confidence AND significant without batch 4), by
#      direction (704 / 262)
#  Background: the 22,371 genes tested in the pooled comparison.
#  Settings as in the main GO analysis (Methods 3.9): clusterProfiler::enrichGO, BP,
#  BH-adjusted p and q < 0.05, redundant terms removed with simplify (cutoff 0.7).
#  LIGHT: reads two saved result files only; no DESeq2 model is loaded.
#  Output: NEW folder Merged_analysis/final_runs/GO_regional_<date>_<time>/
#  Run:    Session > Restart R, Source.  (2-5 minutes)
# =============================================================================

suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2) })
base <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
out_dir <- file.path(base, "final_runs", format(Sys.time(), "GO_regional_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }
latest_dir <- function(prefix, must) {
  d <- list.dirs(file.path(base, "final_runs"), recursive = FALSE)
  d <- d[grepl(paste0("^", prefix, "_[0-9]"), basename(d)) & file.exists(file.path(d, must))]
  if (!length(d)) stop("No final_runs/", prefix, "_* folder containing ", must); sort(d)[length(d)] }

# ---- 1. Inputs ----------------------------------------------------------------------
pcdir <- latest_dir("pooled_cortex", "PooledCortex_all_genes.rds")
fdir  <- latest_dir("voomexact", "PooledRobustness_flags_exact_true_donors.csv")
G  <- readRDS(file.path(pcdir, "PooledCortex_all_genes.rds"))
fl <- read.csv(file.path(fdir, "PooledRobustness_flags_exact_true_donors.csv"))
sigG <- G$region_difference_padj < 0.05
rob  <- fl$gene_id[fl$high_confidence_true_donors %in% TRUE & fl$significant_without_batch4 %in% TRUE]
say("Inputs: ", basename(pcdir), " and ", basename(fdir))
say("Self-check (22371 genes, 1867 region-dependent, 966 robust): ", nrow(G), ", ", sum(sigG), ", ", length(rob))

to_eg <- function(gid) unique(na.omit(unname(mapIds(org.Hs.eg.db, sub("\\..*$", "", gid), "ENTREZID", "ENSEMBL", multiVals = "first"))))
univ <- to_eg(G$gene_id)
neg <- G$striatum_minus_pooled_cortex < 0
lists <- list(
  all_neg = list(lab = "All region-dependent genes, HD effect more negative in striatum", ids = G$gene_id[sigG & neg]),
  all_pos = list(lab = "All region-dependent genes, HD effect more positive in striatum", ids = G$gene_id[sigG & !neg]),
  rob_neg = list(lab = "Robust genes, HD effect more negative in striatum", ids = intersect(rob, G$gene_id[neg])),
  rob_pos = list(lab = "Robust genes, HD effect more positive in striatum", ids = intersect(rob, G$gene_id[!neg])))

# ---- 2. Over-representation analysis ------------------------------------------------------
res <- list(); summ <- list()
for (k in names(lists)) {
  gg <- to_eg(lists[[k]]$ids)
  x <- enrichGO(gg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = univ, minGSSize = 10, maxGSSize = 500,
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  n_all <- if (is.null(x)) 0 else nrow(as.data.frame(x))
  r <- NULL
  if (n_all > 0) {
    write.csv(as.data.frame(x), file.path(out_dir, paste0("GO_BP_", k, "_all_terms.csv")), row.names = FALSE)
    r <- as.data.frame(simplify(x, cutoff = 0.7, by = "p.adjust", select_fun = min))
    write.csv(r, file.path(out_dir, paste0("GO_BP_", k, "_nonredundant.csv")), row.names = FALSE)
  }
  res[[k]] <- r
  summ[[k]] <- data.frame(gene_list = lists[[k]]$lab, genes = length(lists[[k]]$ids), genes_mapped = length(gg),
                          background_mapped = length(univ), enriched_terms = n_all, non_redundant = if (is.null(r)) 0 else nrow(r))
  if (!is.null(r) && nrow(r)) show(paste0(lists[[k]]$lab, ": top 15 non-redundant terms"),
                                   head(data.frame(term = r$Description, genes = r$Count, gene_ratio = r$GeneRatio, padj = signif(r$p.adjust, 2),
                                                   example_genes = sapply(strsplit(r$geneID, "/"), function(z) paste(head(z, 8), collapse = ", "))), 15))
  else say("\n", lists[[k]]$lab, ": no enriched term")
}
summ <- do.call(rbind, summ)
show("Summary (as Table 9)", summ)
write.csv(summ, file.path(out_dir, "GO_BP_regional_summary.csv"), row.names = FALSE)

# ---- 3. Figure: top non-redundant terms of the robust lists (falls back to all genes) --------
pick <- if (!is.null(res$rob_neg) || !is.null(res$rob_pos)) c("rob_neg", "rob_pos") else c("all_neg", "all_pos")
fd <- do.call(rbind, lapply(pick, function(k) { r <- res[[k]]; if (is.null(r) || !nrow(r)) return(NULL); r <- head(r[order(r$p.adjust), ], 10)
data.frame(list = paste0(lists[[k]]$lab, "\n(", length(lists[[k]]$ids), " genes)"), term = r$Description, genes = r$Count, padj = r$p.adjust) }))
if (!is.null(fd) && nrow(fd)) {
  fd$term <- factor(fd$term, levels = rev(unique(fd$term[order(fd$list, fd$padj)])))
  p <- ggplot(fd, aes(list, term)) + geom_point(aes(size = genes, colour = -log10(padj))) +
    scale_colour_gradient(low = "#9ecae1", high = "#08306b", name = "-log10\nadjusted p") +
    scale_size_continuous(range = c(2, 7), name = "Genes") +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 9) + theme(axis.text.x = element_text(size = 8))
  ggsave(file.path(out_dir, "GO_BP_striatum_vs_pooled_cortex.png"), p, width = 8, height = max(3.5, 0.28 * length(unique(fd$term)) + 1.5), dpi = 300)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
say("\nDone. Send Claude the console output and the PNG. Outputs: ", out_dir)