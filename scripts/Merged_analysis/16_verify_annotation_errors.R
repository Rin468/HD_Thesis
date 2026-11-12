# =============================================================================
#  16_verify_annotation_errors.R   (written 1 Nov 2026)
#  Re-check the sex and tissue corrections of the Paryani samples (Table 5) with
#  evidence that the original check did NOT use:
#   1. the corrections as stored in colData (what was changed and why);
#   2. whole-transcriptome similarity: is each sample more similar to striatal or to
#      cortical reference samples over the 2,000 most variable genes? Reference
#      samples are those whose deposited tissue agreed with the marker check, and
#      samples of the SAME donor are excluded as references;
#   3. all samples of each affected donor, side by side (reveals within-donor swaps);
#   4. sex-chromosome genes per sample for the sex-corrected donors.
#  Nothing in the analysis is changed.
#  Output: NEW folder Merged_analysis/final_runs/annotation_check_<date>_<time>/
#  Run:    Session > Restart R, Source.  (2-4 minutes)
# =============================================================================
suppressPackageStartupMessages(library(DESeq2))
base   <- "C:/!Virag.ORI/Bioinfo/szakdoga/Merged_analysis"
run    <- file.path(base, "final_runs/bystudy_20261015_151051")
in_dir <- file.path(base, "results", "objects")
out_dir <- file.path(base, "final_runs", format(Sys.time(), "annotation_check_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_log.txt")
say  <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }
show <- function(title, x) { o <- c("", paste0("-- ", title, " --"), capture.output(print(x, row.names = FALSE)))
cat(o, sep = "\n"); cat(o, sep = "\n", file = log_file, append = TRUE) }
latest_file <- function(dir, pattern) { f <- list.files(dir, pattern = pattern, full.names = TRUE)
if (!length(f)) stop("No file matching ", pattern, " in ", dir); f[which.max(file.info(f)$mtime)] }

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

# ---- 1. Corrections as stored in colData (all 76 Paryani samples incl. the excluded one) ----
cd <- as.data.frame(readRDS(file.path(run, "objects", "colData_final.rds")))
cd <- cd[cd$study == "Paryani", ]
idx <- ifelse(cd$sample_id %in% key$title, match(cd$sample_id, key$title),
              ifelse(cd$sample_id %in% key$gsm, match(cd$sample_id, key$gsm), match(cd$sample_id, key$srr)))
stopifnot(!anyNA(idx))
cd$title <- key$title[idx]; cd$true_donor <- key$true_donor[idx]
geo_tis <- ifelse(grepl("Caud|Accum", as.character(cd$tissue_GEO), ignore.case = TRUE), "Striatum",
                  ifelse(grepl("Cing|Cortex", as.character(cd$tissue_GEO), ignore.case = TRUE), "Cortex", as.character(cd$tissue_GEO)))
cd$geo_tissue <- geo_tis
cd$final_tissue <- as.character(cd$tissue)
say("Values of region_check: ", paste(names(table(cd$region_check)), table(cd$region_check), sep = " = ", collapse = "; "))
changed <- (cd$sex_discordant_with_GEO %in% TRUE) | (cd$geo_tissue != cd$final_tissue) | !(cd$in_model %in% TRUE) |
  !(as.character(cd$region_check) %in% names(which.max(table(cd$region_check))))
tab1 <- data.frame(sample = cd$title, true_donor = cd$true_donor, batch = cd$batch_study,
                   GEO_tissue = as.character(cd$tissue_GEO), final_tissue = cd$final_tissue, region_check = cd$region_check,
                   in_model = cd$in_model, striatal_score = round(cd$striatal_score, 2), cortical_score = round(cd$cortical_score, 2),
                   sex_GEO = cd$sex_reported, sex_expression = cd$sex_inferred, sex_score = round(cd$sex_score, 2))[changed, ]
show("1. Samples corrected, flagged or excluded (as stored in colData)", tab1)
write.csv(tab1, file.path(out_dir, "1_corrections_from_colData.csv"), row.names = FALSE)

# ---- 2. Whole-transcriptome similarity to striatal vs cortical reference samples ----------
txi <- readRDS(latest_file(in_dir, "^merged_txi_Paryani76_Tseng14(_v[0-9]+)?\\.rds$"))
cnt <- round(txi$counts[, rownames(cd)]); rm(txi); invisible(gc())
cnt <- cnt[rowSums(cnt >= 10) >= 10, ]
vs  <- assay(vst(DESeqDataSetFromMatrix(cnt, data.frame(row.names = colnames(cnt), x = rep(1, ncol(cnt))), ~ 1), blind = TRUE))
top <- order(apply(vs, 1, var), decreasing = TRUE)[1:2000]
r   <- cor(vs[top, ], method = "pearson")
ref <- !changed & cd$geo_tissue == cd$final_tissue        # samples whose deposited tissue agreed with the marker check
sim <- t(sapply(seq_len(nrow(cd)), function(i) {
  ok <- ref & cd$true_donor != cd$true_donor[i]           # never use the same donor as reference
  c(striatum = mean(r[i, ok & cd$final_tissue == "Striatum"]), cortex = mean(r[i, ok & cd$final_tissue == "Cortex"])) }))
tab2 <- data.frame(sample = cd$title, true_donor = cd$true_donor, GEO_tissue = as.character(cd$tissue_GEO), final_tissue = cd$final_tissue,
                   mean_r_to_striatum = round(sim[, "striatum"], 3), mean_r_to_cortex = round(sim[, "cortex"], 3),
                   transcriptome_says = ifelse(sim[, "striatum"] > sim[, "cortex"], "Striatum", "Cortex"))
tab2$agrees_with_final <- tab2$transcriptome_says == tab2$final_tissue
tab2$agrees_with_GEO   <- tab2$transcriptome_says == cd$geo_tissue
show("2. Whole-transcriptome check of the corrected, flagged and excluded samples", tab2[changed, ])
say("All other samples: transcriptome agrees with their tissue in ", sum(tab2$agrees_with_final[!changed]), " of ", sum(!changed))
write.csv(tab2, file.path(out_dir, "2_transcriptome_similarity_all_samples.csv"), row.names = FALSE)

# ---- 3. All samples of each affected donor ----------------------------------------------
dons <- unique(cd$true_donor[changed])
tab3 <- data.frame(true_donor = cd$true_donor, sample = cd$title, batch = cd$batch_study, GEO_tissue = as.character(cd$tissue_GEO),
                   final_tissue = cd$final_tissue, transcriptome_says = tab2$transcriptome_says,
                   sex_GEO = cd$sex_reported, sex_expression = cd$sex_inferred)[cd$true_donor %in% dons, ]
show("3. All samples of the affected donors (look for swaps: two samples of one donor relabelled in opposite directions)",
     tab3[order(tab3$true_donor, tab3$sample), ])
write.csv(tab3, file.path(out_dir, "3_affected_donors_all_samples.csv"), row.names = FALSE)

# ---- 4. Sex-chromosome genes for the sex-corrected donors ---------------------------------
ann <- as.data.frame(readRDS(latest_file(in_dir, "^gene_annotation_GENCODEv44(_v[0-9]+)?\\.rds$")))
sexg <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY")
cpm <- log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)
ids <- ann$gene_id[match(sexg, ann$gene_symbol)]; ids <- ids[ids %in% rownames(cpm)]
sd_dons <- unique(cd$true_donor[cd$sex_discordant_with_GEO %in% TRUE])
sel <- cd$true_donor %in% sd_dons
tab4 <- data.frame(sample = cd$title[sel], true_donor = cd$true_donor[sel], sex_GEO = cd$sex_reported[sel],
                   round(t(cpm[ids, sel, drop = FALSE]), 2), check.names = FALSE)
names(tab4)[4:ncol(tab4)] <- ann$gene_symbol[match(ids, ann$gene_id)]
show("4. Sex-chromosome genes, log2(CPM + 1), all samples of the sex-corrected donors (female: XIST high, Y genes ~0)", tab4)
write.csv(tab4, file.path(out_dir, "4_sex_genes_sex_corrected_donors.csv"), row.names = FALSE)
say("\nAlso check in the donor table (Supplementary Data 1) the Gender column of donors: ", paste(sd_dons, collapse = ", "))
say("Done. Send Claude the console output. Outputs: ", out_dir)