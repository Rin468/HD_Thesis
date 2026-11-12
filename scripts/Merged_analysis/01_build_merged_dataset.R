# =============================================================================
#  01_build_merged_dataset.R
#  Stage 1 - Merge Paryani 2024 (GSE242195, n = 76) and Tseng (GSE79666, n = 14)
#            tximport objects into one dataset for a joint DESeq2 analysis.
#
#  No statistical model is fitted here. The script:
#   1. locates and loads the four input .rds files
#   2. validates both txi objects (gene IDs, matrix structure, count type)
#   3. harmonises the sample metadata into one colData table
#   4. merges the abundance / counts / length matrices
#   5. tabulates the sample structure (groups, donors, batches, sex, grade, age)
#   6. infers sex from XIST and Y-linked genes (required for Tseng)
#   7. sets the final sex covariate and checks that the planned design is estimable
#   8. (optional) reads Salmon run metadata (meta_info.json) if present on disk
#   9. saves all outputs without overwriting earlier files (_v2, _v3, ... added)
#
#  Save as: Merged_analysis/scripts/01_build_merged_dataset.R
#  Run it with the "Source" button (Ctrl+Shift+S): the script then stops at the
#  first failed check. Everything printed is also written to a log file.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- 0. Paths, labels, helpers ----------------------------------------------
base_dir   <- "C:/!Virag.ORI/Bioinfo/szakdoga"
merged_dir <- file.path(base_dir, "Merged_analysis")
obj_dir    <- file.path(merged_dir, "results", "objects")
tab_dir    <- file.path(merged_dir, "results", "tables")
log_dir    <- file.path(merged_dir, "results", "logs")
fig_dir    <- file.path(merged_dir, "figures")
for (d in c(obj_dir, tab_dir, log_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

lab_P <- "Paryani"   # GSE242195
lab_T <- "Tseng"     # GSE79666 (study label used in colData, batch and donor IDs)

# Numbers verified earlier / reported in the papers - used only to raise flags
expected_groups   <- c(Striatum_CTRL = 10, Striatum_HD = 31, Cortex_CTRL = 17, Cortex_HD = 32)
expected_donors_P <- c(CTRL = 11, HD = 24)                       # Paryani et al. 2024
expected_sex_T    <- rbind(CTRL = c(female = 1, male = 6),       # GSE79666 paper, Table 1
                           HD   = c(female = 3, male = 4))

log_lines <- character(0)   # everything printed is collected here -> log file
checks    <- list()         # all check tables -> RDS
flags     <- character(0)   # issues to review, printed at the end

emit <- function(lines) {
  cat(lines, sep = "\n")
  log_lines <<- c(log_lines, lines)
}
section <- function(title) emit(c("", strrep("=", 78), title, strrep("=", 78)))
say <- function(...) emit(paste0(...))
show_check <- 
  function(name, x) {
  checks[[name]] <<- x
  emit(c("", paste0("-- ", name, " --"), capture.output(print(x))))
  invisible(x)
}
assert_values <- function(x, allowed, what) {
  if (anyNA(x)) stop(what, ": contains missing values")
  bad <- setdiff(unique(x), allowed)
  if (length(bad) > 0) stop(what, ": unexpected value(s): ", paste(bad, collapse = ", "))
  invisible(TRUE)
}
# Returns a path that does not overwrite an existing file (adds _v2, _v3, ...)
safe_path <- function(dir, filename) {
  path <- file.path(dir, filename)
  if (!file.exists(path)) return(path)
  stem <- tools::file_path_sans_ext(filename)
  ext  <- tools::file_ext(filename)
  i <- 2
  while (file.exists(file.path(dir, sprintf("%s_v%d.%s", stem, i, ext)))) i <- i + 1
  file.path(dir, sprintf("%s_v%d.%s", stem, i, ext))
}


# ---- 1. Locate and load the input files -------------------------------------
section("1. Locate and load input files")

input_files <- c(paryani_bundle = "GSE242195_analysis_bundle.rds",
                 paryani_meta   = "metadata_76.rds",
                 tseng_bundle   = "GSE79666_analysis_bundle.rds",
                 tseng_meta     = "Tseng_metadata.rds")

scan_files <- list.files(base_dir, pattern = "\\.(rds|json)$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)

input_paths <- vapply(names(input_files), function(k) {
  hits <- scan_files[tolower(basename(scan_files)) == tolower(input_files[[k]])]
  if (length(hits) == 0) stop("File not found under base_dir: ", input_files[[k]])
  if (length(hits) > 1) {
    print(file.info(hits)[, c("size", "mtime")])
    stop("Several copies of ", input_files[[k]], " were found (listed above). ",
         "Tell me which one is correct before continuing.")
  }
  hits
}, character(1))
show_check("input_files", data.frame(input = names(input_paths), path = unname(input_paths)))

par_bundle <- readRDS(input_paths[["paryani_bundle"]])
par_meta   <- as.data.frame(readRDS(input_paths[["paryani_meta"]]), stringsAsFactors = FALSE)
tse_bundle <- readRDS(input_paths[["tseng_bundle"]])
tse_meta   <- as.data.frame(readRDS(input_paths[["tseng_meta"]]), stringsAsFactors = FALSE)
par_txi <- par_bundle$txi
tse_txi <- tse_bundle$txi

say("Paryani bundle elements: ", paste(names(par_bundle), collapse = ", "),
    "  | txi elements: ", paste(names(par_txi), collapse = ", "))
say("Tseng bundle elements:   ", paste(names(tse_bundle), collapse = ", "),
    "  | txi elements: ", paste(names(tse_txi), collapse = ", "))

# Are the stand-alone metadata files identical to the copies stored in the bundles?
same_as_bundle <- function(a, b) {
  if (is.null(b)) return("no copy in bundle")
  b <- as.data.frame(b, stringsAsFactors = FALSE)
  if (identical(dim(a), dim(b)) && isTRUE(all.equal(a, b, check.attributes = FALSE))) {
    "identical content"
  } else {
    paste0("DIFFERENT (file ", paste(dim(a), collapse = " x "),
           ", bundle ", paste(dim(b), collapse = " x "), ")")
  }
}
meta_cmp <- data.frame(dataset = c(lab_P, lab_T),
                       file_vs_bundle_copy = c(same_as_bundle(par_meta, par_bundle$sample_metadata),
                                               same_as_bundle(tse_meta, tse_bundle$metadata)))
show_check("metadata_file_vs_bundle_copy", meta_cmp)
if (any(grepl("DIFFERENT", meta_cmp$file_vs_bundle_copy)))
  flags <- c(flags, "A metadata file differs from the copy stored in its bundle (the stand-alone file is used).")

if (!is.null(tse_bundle$study_information))
  emit(c("", "-- Tseng study_information (as stored in the bundle) --",
         capture.output(str(tse_bundle$study_information, max.level = 2))))


# ---- 2. Validate the tximport objects ---------------------------------------
section("2. Validate the tximport objects")

validate_txi <- function(txi, label) {
  mats <- c("abundance", "counts", "length")
  miss <- setdiff(mats, names(txi))
  if (length(miss) > 0) stop(label, ": txi lacks ", paste(miss, collapse = ", "))
  for (m in mats) {
    if (!is.matrix(txi[[m]]) || !is.numeric(txi[[m]])) stop(label, ": txi$", m, " is not a numeric matrix")
    if (!identical(dimnames(txi[[m]]), dimnames(txi$counts))) stop(label, ": txi$", m, " has different dimnames from txi$counts")
    if (anyNA(txi[[m]])) stop(label, ": txi$", m, " contains NA/NaN")
  }
  if (any(txi$counts < 0) || any(txi$abundance < 0)) stop(label, ": negative counts or abundance")
  if (any(txi$length <= 0)) stop(label, ": non-positive values in txi$length")
  if (anyDuplicated(rownames(txi$counts)) > 0) stop(label, ": duplicated gene IDs")
  if (anyDuplicated(colnames(txi$counts)) > 0) stop(label, ": duplicated sample IDs")
  tpm_sum <- colSums(txi$abundance)
  data.frame(dataset             = label,
             genes               = nrow(txi$counts),
             samples             = ncol(txi$counts),
             countsFromAbundance = if (is.null(txi$countsFromAbundance)) "<not stored>" else as.character(txi$countsFromAbundance),
             counts_all_integer  = all(txi$counts == round(txi$counts)),
             TPM_colSum_min      = round(min(tpm_sum)),
             TPM_colSum_max      = round(max(tpm_sum)),
             versioned_ENSG_ids  = all(grepl("^ENSG[0-9]+\\.[0-9]+", rownames(txi$counts))),
             first_gene_id       = rownames(txi$counts)[1])
}
show_check("txi_structure", rbind(validate_txi(par_txi, lab_P), validate_txi(tse_txi, lab_T)))

# Gene IDs must be identical (same annotation, same order) before cbind()
ids_P <- rownames(par_txi$counts)
ids_T <- rownames(tse_txi$counts)
if (identical(ids_P, ids_T)) {
  say("Gene IDs: identical and in the same order in both datasets (", length(ids_P), " genes).")
} else if (setequal(ids_P, ids_T)) {
  for (m in c("abundance", "counts", "length")) tse_txi[[m]] <- tse_txi[[m]][ids_P, , drop = FALSE]
  say("Gene IDs: same set but different order -> Tseng matrices reordered to the Paryani order.")
  flags <- c(flags, "Gene order differed between the datasets (Tseng reordered).")
} else {
  say("Genes only in Paryani: ", length(setdiff(ids_P, ids_T)),
      " (e.g. ", paste(head(setdiff(ids_P, ids_T), 5), collapse = ", "), ")")
  say("Genes only in Tseng:   ", length(setdiff(ids_T, ids_P)),
      " (e.g. ", paste(head(setdiff(ids_T, ids_P), 5), collapse = ", "), ")")
  stop("Gene ID sets differ between the datasets - stop here and send me this output.")
}

# countsFromAbundance: the stored label ...
if (!identical(as.character(par_txi$countsFromAbundance), "no"))
  stop("Paryani txi$countsFromAbundance is not 'no'")
if (!is.null(tse_txi$countsFromAbundance) && !identical(as.character(tse_txi$countsFromAbundance), "no"))
  stop("Tseng txi$countsFromAbundance is not 'no'")

# ... and an empirical check that works even if the label is not stored.
# Within one sample, tximport countsFromAbundance = "no" implies
#   counts = constant * abundance * length                (exact up to Salmon rounding)
# "lengthScaledTPM" implies counts = constant * abundance * rowMeans(length),
# "scaledTPM"       implies counts = constant * abundance.
# The hypothesis whose log-ratio has ~zero spread (MAD) across genes is the true one.
cfa_diag <- function(txi, label, min_count = 100) {
  avg_len <- rowMeans(txi$length)
  per_sample <- vapply(seq_len(ncol(txi$counts)), function(j) {
    cnt <- txi$counts[, j]
    ab  <- txi$abundance[, j]
    len <- txi$length[, j]
    k   <- cnt >= min_count & ab > 0
    c(no              = mad(log(cnt[k] / (ab[k] * len[k]))),
      lengthScaledTPM = mad(log(cnt[k] / (ab[k] * avg_len[k]))),
      scaledTPM       = mad(log(cnt[k] / ab[k])))
  }, numeric(3))
  med <- apply(per_sample, 1, median)
  data.frame(dataset                = label,
             MAD_if_no              = signif(med[["no"]], 3),
             MAD_if_lengthScaledTPM = signif(med[["lengthScaledTPM"]], 3),
             MAD_if_scaledTPM       = signif(med[["scaledTPM"]], 3),
             best_match             = names(med)[which.min(med)])
}
cfa_tab <- rbind(cfa_diag(par_txi, lab_P), cfa_diag(tse_txi, lab_T))
show_check("count_type_empirical_check (median over samples)", cfa_tab)
if (any(cfa_tab$best_match != "no" | cfa_tab$MAD_if_no > 0.01))
  stop("The empirical check does not support countsFromAbundance = 'no' for both datasets - send me this output.")
if (is.null(tse_txi$countsFromAbundance)) {
  say("Tseng txi has no countsFromAbundance element; the data are consistent with 'no' -> set to 'no' in the merged object.")
  flags <- c(flags, "Tseng txi did not store countsFromAbundance; 'no' was inferred from the data.")
}


# ---- 3. Harmonise the sample metadata ---------------------------------------
section("3. Harmonise sample metadata")

need_P <- c("title.x", "Sample.Name", "tissue", "group", "batch", "sex", "grade", "age_binned", "cag_repeat")
need_T <- c("run_accession", "geo_accession", "sample_name", "group")
if (length(setdiff(need_P, names(par_meta))) > 0)
  stop("Paryani metadata lacks column(s): ", paste(setdiff(need_P, names(par_meta)), collapse = ", "))
if (length(setdiff(need_T, names(tse_meta))) > 0)
  stop("Tseng metadata lacks column(s): ", paste(setdiff(need_T, names(tse_meta)), collapse = ", "))

# Find the metadata column that holds the txi column names and put the rows in txi order
align_meta <- function(meta, txi, preferred, label) {
  ids_txi <- colnames(txi$counts)
  hits <- names(meta)[vapply(meta, function(x) all(ids_txi %in% as.character(x)), logical(1))]
  if (length(hits) == 0) stop(label, ": no metadata column contains all txi column names")
  idcol <- if (preferred %in% hits) preferred else hits[1]
  ids <- as.character(meta[[idcol]])
  if (anyDuplicated(ids) > 0) stop(label, ": duplicated values in metadata column ", idcol)
  say(label, ": txi columns matched to metadata column '", idcol, "' (",
      nrow(meta), " metadata rows, ", length(ids_txi), " txi columns)")
  out <- meta[match(ids_txi, ids), , drop = FALSE]
  stopifnot(identical(as.character(out[[idcol]]), ids_txi))
  out
}
par_meta <- align_meta(par_meta, par_txi, "Sample.Name", lab_P)
tse_meta <- align_meta(tse_meta, tse_txi, "run_accession", lab_T)

# Paryani: clean and validate the key columns
P_title  <- trimws(as.character(par_meta$title.x))
P_region <- trimws(as.character(par_meta$tissue))
P_cond   <- toupper(trimws(as.character(par_meta$group)))
P_batch  <- gsub("[^0-9]", "", as.character(par_meta$batch))
P_sex    <- tolower(trimws(as.character(par_meta$sex)))
P_grade  <- trimws(as.character(par_meta$grade))
assert_values(P_region, c("Accumbens", "Caudate", "Cingulate"), "Paryani tissue")
assert_values(P_cond,   c("CTRL", "HD"),                        "Paryani group")
assert_values(P_batch,  c("1", "2", "3", "4"),                  "Paryani batch")
assert_values(P_sex,    c("female", "male"),                    "Paryani sex")

P_donor_raw <- sub("_.*$", "", P_title)          # "H4161_Caud" -> "H4161"
P_prefix    <- substr(P_donor_raw, 1, 1)
show_check("Paryani_donor_prefix_vs_group", table(prefix = P_prefix, group = P_cond))
n_bad_prefix <- sum(!((P_prefix == "H" & P_cond == "HD") | (P_prefix == "C" & P_cond == "CTRL")))
if (n_bad_prefix > 0)
  flags <- c(flags, paste0(n_bad_prefix, " Paryani sample(s): donor prefix (H/C) does not match group."))

# Tseng: clean and validate
T_name <- trimws(as.character(tse_meta$sample_name))
T_cond <- toupper(trimws(as.character(tse_meta$group)))
assert_values(T_cond, c("CTRL", "HD"), "Tseng group")
T_prefix <- toupper(sub("-.*$", "", T_name))    # "HD-57" -> "HD"
show_check("Tseng_sample_name_prefix_vs_group", table(prefix = T_prefix, group = T_cond))
if (any(T_prefix != T_cond))
  flags <- c(flags, "Tseng: sample_name prefix does not match group for some samples.")

# Combined colData (rows in the same order as the txi columns: Paryani, then Tseng)
cd <- rbind(
  data.frame(sample_id    = as.character(par_meta$Sample.Name),
             txi_colname  = colnames(par_txi$counts),
             sample_label = P_title,
             study        = lab_P,
             geo_series   = "GSE242195",
             donor        = paste0("P_", P_donor_raw),
             region       = P_region,
             condition    = P_cond,
             batch_study  = paste0(lab_P, "_b", P_batch),
             sex_reported = P_sex,
             grade        = P_grade,
             age_binned   = as.character(par_meta$age_binned),
             cag_repeat   = as.character(par_meta$cag_repeat),
             stringsAsFactors = FALSE),
  data.frame(sample_id    = as.character(tse_meta$geo_accession),
             txi_colname  = colnames(tse_txi$counts),
             sample_label = T_name,
             study        = lab_T,
             geo_series   = "GSE79666",
             donor        = paste0("T_", T_name),
             region       = "BA4",
             condition    = T_cond,
             batch_study  = lab_T,
             sex_reported = NA_character_,
             grade        = NA_character_,
             age_binned   = NA_character_,
             cag_repeat   = NA_character_,
             stringsAsFactors = FALSE))
if (anyNA(cd$sample_id) || anyDuplicated(cd$sample_id) > 0) stop("sample_id (GSM) is missing or duplicated")
if (!all(grepl("^GSM[0-9]+$", cd$sample_id))) flags <- c(flags, "Some sample_id values are not GSM accessions.")

cd$tissue         <- ifelse(cd$region %in% c("Caudate", "Accumbens"), "Striatum", "Cortex")
cd$group          <- paste(cd$tissue, cd$condition, sep = "_")
cd$juvenile_onset <- ifelse(cd$study == lab_P, cd$grade == "J", NA)

cd$study       <- factor(cd$study,       levels = c(lab_P, lab_T))
cd$region      <- factor(cd$region,      levels = c("Caudate", "Accumbens", "Cingulate", "BA4"))
cd$tissue      <- factor(cd$tissue,      levels = c("Striatum", "Cortex"))
cd$condition   <- factor(cd$condition,   levels = c("CTRL", "HD"))        # CTRL = reference
cd$group       <- factor(cd$group,       levels = c("Striatum_CTRL", "Striatum_HD", "Cortex_CTRL", "Cortex_HD"))
cd$batch_study <- factor(cd$batch_study, levels = c(paste0(lab_P, "_b", 1:4), lab_T))
rownames(cd)   <- cd$sample_id
stopifnot(!anyNA(cd$study), !anyNA(cd$region), !anyNA(cd$tissue),
          !anyNA(cd$condition), !anyNA(cd$group), !anyNA(cd$batch_study))
say("colData built: ", nrow(cd), " samples")


# ---- 4. Merge the tximport objects ------------------------------------------
section("4. Merge the tximport objects")

stopifnot(identical(cd$txi_colname, c(colnames(par_txi$counts), colnames(tse_txi$counts))))
merged_txi <- list()
for (m in c("abundance", "counts", "length")) {
  merged_txi[[m]] <- cbind(par_txi[[m]], tse_txi[[m]])
  colnames(merged_txi[[m]]) <- cd$sample_id
}
merged_txi$countsFromAbundance <- "no"

isP <- cd$study == lab_P
isT <- cd$study == lab_T
stopifnot(identical(colnames(merged_txi$counts), rownames(cd)),
          identical(dimnames(merged_txi$abundance), dimnames(merged_txi$counts)),
          identical(dimnames(merged_txi$length),    dimnames(merged_txi$counts)),
          all(merged_txi$counts[, isP] == par_txi$counts),
          all(merged_txi$counts[, isT] == tse_txi$counts),
          all(merged_txi$length > 0))
say("Merged txi: ", nrow(merged_txi$counts), " genes x ", ncol(merged_txi$counts),
    " samples; countsFromAbundance = 'no'")

cd$lib_size_counts   <- round(colSums(merged_txi$counts))   # estimated counts assigned to genes
cd$genes_counts_ge10 <- colSums(merged_txi$counts >= 10)
show_check("library_size_millions_by_study",
           sapply(split(cd$lib_size_counts / 1e6, cd$study), function(x) round(summary(x), 1)))
show_check("genes_with_counts_ge10_by_study",
           sapply(split(cd$genes_counts_ge10, cd$study), function(x) round(summary(x))))


# ---- 5. Sample-structure checks ---------------------------------------------
section("5. Sample-structure checks")

grp_tab <- table(group = cd$group, study = cd$study)
show_check("group_x_study", grp_tab)
if (!all(rowSums(grp_tab)[names(expected_groups)] == expected_groups))
  flags <- c(flags, "Group sizes differ from the expected 10 / 31 / 17 / 32.")
show_check("region_x_condition", table(region = cd$region, condition = cd$condition))
show_check("batch_x_group", table(batch = cd$batch_study, group = cd$group))
show_check("batch_x_region", table(batch = cd$batch_study, region = cd$region))
show_check("Paryani_sex_reported_x_condition", table(sex = cd$sex_reported[isP], condition = cd$condition[isP]))
show_check("Paryani_grade_x_group", table(grade = cd$grade[isP], group = droplevels(cd$group[isP])))
show_check("Paryani_grade_x_batch", table(grade = cd$grade[isP], batch = droplevels(cd$batch_study[isP])))

# Donor-level summary
attr_cols <- c("condition", "sex_reported", "grade", "age_binned", "cag_repeat")
donor_tab <- do.call(rbind, lapply(split(cd, cd$donor), function(d) {
  n_distinct <- vapply(attr_cols, function(a) length(unique(as.character(d[[a]]))), integer(1))
  data.frame(donor               = d$donor[1],
             study               = as.character(d$study[1]),
             condition           = as.character(d$condition[1]),
             n_samples           = nrow(d),
             regions             = paste(as.character(sort(d$region)), collapse = "+"),
             n_regions           = length(unique(d$region)),
             batches             = paste(unique(as.character(sort(d$batch_study))), collapse = "+"),
             n_batches           = length(unique(d$batch_study)),
             sex_reported        = as.character(d$sex_reported[1]),
             grade               = as.character(d$grade[1]),
             age_binned          = as.character(d$age_binned[1]),
             cag_repeat          = as.character(d$cag_repeat[1]),
             inconsistent_fields = paste(attr_cols[n_distinct > 1], collapse = ","),
             stringsAsFactors = FALSE)
}))
rownames(donor_tab) <- NULL
dP <- donor_tab$study == lab_P

show_check("donors_per_study_x_condition", table(study = donor_tab$study, condition = donor_tab$condition))
show_check("samples_per_donor", table(study_condition = paste(donor_tab$study, donor_tab$condition),
                                      n_samples = donor_tab$n_samples))
show_check("Paryani_region_combinations_per_donor",
           table(regions = donor_tab$regions[dP], condition = donor_tab$condition[dP]))

dup_region <- donor_tab[donor_tab$n_samples > donor_tab$n_regions, ]
show_check("donors_with_more_than_one_sample_from_the_same_region",
           if (nrow(dup_region) > 0) cd[cd$donor %in% dup_region$donor, c("sample_label", "donor", "region", "batch_study", "condition")] else "none")
multi_batch <- donor_tab[donor_tab$n_batches > 1, ]
show_check("donors_spanning_more_than_one_batch",
           if (nrow(multi_batch) > 0) multi_batch[, c("donor", "condition", "regions", "batches")] else "none")
inconsistent <- donor_tab[donor_tab$inconsistent_fields != "", ]
show_check("donors_with_inconsistent_donor_level_metadata",
           if (nrow(inconsistent) > 0) cd[cd$donor %in% inconsistent$donor, c("sample_label", "donor", attr_cols)] else "none")

n_donor_P <- table(factor(donor_tab$condition[dP], levels = c("CTRL", "HD")))
if (!all(n_donor_P == expected_donors_P))
  flags <- c(flags, sprintf("Paryani donors derived from title.x: %d CTRL / %d HD (paper: %d / %d).",
                            n_donor_P[["CTRL"]], n_donor_P[["HD"]], expected_donors_P[["CTRL"]], expected_donors_P[["HD"]]))
if (nrow(dup_region) > 0)
  flags <- c(flags, paste0(nrow(dup_region), " donor(s) with >1 sample from the same region (technical replicates?)."))
if (nrow(inconsistent) > 0)
  flags <- c(flags, paste0(nrow(inconsistent), " donor(s) with inconsistent donor-level metadata."))

show_check("Paryani_age_binned_x_condition (DONOR level)",
           table(age_binned = donor_tab$age_binned[dP], condition = donor_tab$condition[dP]))
show_check("Paryani_age_binned_x_group (SAMPLE level)",
           table(age_binned = cd$age_binned[isP], group = droplevels(cd$group[isP])))


# ---- 6. Sex inference from XIST and Y-linked genes --------------------------
section("6. Sex inference from XIST and Y-linked genes")

tx2g <- as.data.frame(tse_bundle$tx2gene, stringsAsFactors = FALSE)
if (!all(c("gene_id", "gene_symbol") %in% names(tx2g))) stop("Tseng tx2gene lacks gene_id / gene_symbol columns")
ann <- unique(data.frame(gene_id     = as.character(tx2g$gene_id),
                         gene_symbol = as.character(tx2g$gene_symbol),
                         stringsAsFactors = FALSE))
n_multi_symbol <- sum(duplicated(ann$gene_id))
ann <- ann[!duplicated(ann$gene_id), ]
gene_annot <- data.frame(gene_id     = rownames(merged_txi$counts),
                         gene_symbol = ann$gene_symbol[match(rownames(merged_txi$counts), ann$gene_id)],
                         stringsAsFactors = FALSE)
say("Gene annotation (Tseng tx2gene): ", sum(!is.na(gene_annot$gene_symbol)), " of ", nrow(gene_annot),
    " genes have a symbol; gene IDs with >1 symbol (first kept): ", n_multi_symbol)

# Marker genes; the unversioned Ensembl IDs are only used to cross-check the symbol mapping
sex_marker_ids <- c(XIST = "ENSG00000229807", RPS4Y1 = "ENSG00000129824", DDX3Y = "ENSG00000067048",
                    KDM5D = "ENSG00000012817", UTY = "ENSG00000183878")
y_genes <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY")
sex_ids <- vapply(names(sex_marker_ids), function(s) {
  hit <- gene_annot$gene_id[gene_annot$gene_symbol %in% s & !grepl("_PAR_Y$", gene_annot$gene_id)]
  if (length(hit) != 1) stop("Gene symbol ", s, " maps to ", length(hit), " gene IDs: ", paste(hit, collapse = ", "))
  hit
}, character(1))
id_ok <- unname(sub("\\..*$", "", sex_ids) == sex_marker_ids)
show_check("sex_marker_gene_ids", data.frame(symbol = names(sex_ids), gene_id = unname(sex_ids),
                                             matches_expected_Ensembl_ID = id_ok))
if (!all(id_ok)) flags <- c(flags, "A sex-marker symbol mapped to an unexpected Ensembl ID.")

# log2(CPM + 1) of the marker genes (library size = all estimated counts of the sample)
lib_size   <- colSums(merged_txi$counts)
sex_logcpm <- log2(sweep(merged_txi$counts[sex_ids, , drop = FALSE], 2, lib_size, "/") * 1e6 + 1)
rownames(sex_logcpm) <- names(sex_ids)
y_mean    <- colMeans(sex_logcpm[y_genes, , drop = FALSE])
sex_score <- y_mean - sex_logcpm["XIST", ]   # > 0: male-like, < 0: female-like

# Decision rule: |score| > 1 log2 unit -> call; otherwise "ambiguous" (inspect manually)
score_margin <- 1
sex_tab <- data.frame(cd[, c("sample_id", "sample_label", "study", "donor", "region", "condition", "sex_reported")],
                      round(t(sex_logcpm), 3), check.names = FALSE, stringsAsFactors = FALSE)
sex_tab$Y_mean       <- round(y_mean, 3)
sex_tab$sex_score    <- round(sex_score, 3)
sex_tab$sex_inferred <- ifelse(sex_score >  score_margin, "male",
                               ifelse(sex_score < -score_margin, "female", "ambiguous"))
# Informational: both markers high (e.g. mixture, XXY) or both low (e.g. low quality, loss of Y)
sex_tab$marker_flag  <- ifelse(sex_tab$XIST > 2 & sex_tab$Y_mean > 2, "XIST and Y both high",
                               ifelse(sex_tab$XIST < 1 & sex_tab$Y_mean < 1, "XIST and Y both low", ""))
sex_tab$concordant_with_reported <- ifelse(is.na(sex_tab$sex_reported), NA,
                                           sex_tab$sex_reported == sex_tab$sex_inferred)

show_check("marker_separation_by_inferred_sex (log2 CPM+1 ranges)",
           do.call(rbind, lapply(split(sex_tab, sex_tab$sex_inferred), function(d) data.frame(
             sex_inferred = d$sex_inferred[1], n = nrow(d),
             XIST_range   = paste(range(d$XIST), collapse = " to "),
             Ymean_range  = paste(range(d$Y_mean), collapse = " to "),
             score_range  = paste(range(d$sex_score), collapse = " to ")))))
show_check("Paryani_reported_vs_inferred_sex",
           table(reported = sex_tab$sex_reported[isP], inferred = sex_tab$sex_inferred[isP]))
show_check("Tseng_marker_expression_and_inferred_sex",
           sex_tab[isT, c("sample_label", "condition", "XIST", y_genes, "Y_mean", "sex_score", "sex_inferred")])
sexT_tab <- table(condition = sex_tab$condition[isT],
                  inferred  = factor(sex_tab$sex_inferred[isT], levels = c("female", "male", "ambiguous")))
show_check("Tseng_inferred_sex_x_condition", sexT_tab)
show_check("Tseng_expected_sex_from_paper_Table1", expected_sex_T)
if (!all(sexT_tab[c("CTRL", "HD"), c("female", "male")] == expected_sex_T[c("CTRL", "HD"), c("female", "male")]))
  flags <- c(flags, "Tseng inferred sex distribution differs from the paper's Table 1 (CTRL 6M:1F, HD 4M:3F).")

to_review <- sex_tab[sex_tab$sex_inferred == "ambiguous" | sex_tab$concordant_with_reported %in% FALSE |
                       sex_tab$marker_flag != "", ]
show_check("samples_to_review_for_sex",
           if (nrow(to_review) > 0) to_review[, c("sample_label", "study", "donor", "sex_reported", "XIST",
                                                  "Y_mean", "sex_score", "sex_inferred", "marker_flag")] else "none")

n_sex_per_donor <- tapply(sex_tab$sex_inferred, sex_tab$donor, function(x) length(unique(x)))
bad_donors <- names(n_sex_per_donor)[n_sex_per_donor > 1]
show_check("donors_whose_samples_differ_in_inferred_sex",
           if (length(bad_donors) > 0) sex_tab[sex_tab$donor %in% bad_donors,
                                               c("sample_label", "donor", "region", "sex_reported", "sex_inferred")] else "none")

if (sum(sex_tab$concordant_with_reported %in% FALSE) > 0)
  flags <- c(flags, paste0(sum(sex_tab$concordant_with_reported %in% FALSE),
                           " Paryani sample(s): reported sex discordant with expression."))
if (sum(sex_tab$sex_inferred == "ambiguous") > 0)
  flags <- c(flags, paste0(sum(sex_tab$sex_inferred == "ambiguous"), " sample(s) with ambiguous inferred sex."))
if (length(bad_donors) > 0)
  flags <- c(flags, paste0(length(bad_donors), " donor(s) whose samples differ in inferred sex (possible sample swap)."))
if (any(sex_tab$marker_flag != ""))
  flags <- c(flags, "Some samples show an atypical XIST / Y-marker combination.")

# Plot: XIST vs mean Y-marker expression
sex_tab$reported_for_plot <- ifelse(is.na(sex_tab$sex_reported), paste0("not reported (", lab_T, ")"),
                                    sex_tab$sex_reported)
lab_dat <- sex_tab[isT | sex_tab$sex_inferred == "ambiguous" | sex_tab$concordant_with_reported %in% FALSE, ]
p_sex <- ggplot(sex_tab, aes(x = XIST, y = Y_mean)) +
  geom_abline(intercept = c(-score_margin, score_margin), slope = 1, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = reported_for_plot, shape = study), size = 2.6, alpha = 0.85) +
  geom_text(data = lab_dat, aes(label = sample_label), size = 2.3, vjust = -1, colour = "grey20") +
  labs(x = "XIST  [log2(CPM + 1)]",
       y = "Mean of RPS4Y1, DDX3Y, KDM5D, UTY  [log2(CPM + 1)]",
       colour = "Reported sex", shape = "Study",
       title = "Sex inference from sex-chromosome marker genes",
       subtitle = paste0("Dashed lines: Y mean - XIST = +/-", score_margin,
                         " (samples between the lines = ambiguous)")) +
  theme_bw(base_size = 11)
tryCatch(print(p_sex), error = function(e)
  say("Plot not shown in the Plots pane (it is still saved to file): ", conditionMessage(e)))


# ---- 7. Final sex covariate and design check --------------------------------
section("7. Final sex covariate and design check")

stopifnot(identical(sex_tab$sample_id, cd$sample_id))
cd$sex_inferred <- sex_tab$sex_inferred
cd$sex_score    <- sex_tab$sex_score
cd$sex_source   <- ifelse(isP, "reported", "inferred_from_expression")
cd$sex <- factor(ifelse(isP, cd$sex_reported,
                        ifelse(cd$sex_inferred %in% c("female", "male"), cd$sex_inferred, NA)),
                 levels = c("female", "male"))
show_check("final_sex_x_group_x_study",
           ftable(table(study = cd$study, group = cd$group, sex = cd$sex, useNA = "ifany")))

if (anyNA(cd$sex)) {
  flags <- c(flags, paste0("sex is NA for ", sum(is.na(cd$sex)), " sample(s) - must be resolved before Stage 2."))
  say("Design check skipped: sex is NA for ", paste(cd$sample_label[is.na(cd$sex)], collapse = ", "))
} else {
  mm <- model.matrix(~ batch_study + sex + group, data = cd)
  rk <- qr(mm)$rank
  say("Design ~ batch_study + sex + group: ", ncol(mm), " coefficients, rank ", rk,
      if (rk == ncol(mm)) " -> full rank (all coefficients estimable)" else " -> NOT full rank")
  say("Coefficients: ", paste(colnames(mm), collapse = ", "))
  if (rk < ncol(mm)) flags <- c(flags, "The planned design matrix is not full rank.")
}

col_order <- c("sample_id", "sample_label", "txi_colname", "study", "geo_series", "donor",
               "region", "tissue", "condition", "group", "batch_study",
               "sex", "sex_source", "sex_reported", "sex_inferred", "sex_score",
               "grade", "juvenile_onset", "age_binned", "cag_repeat",
               "lib_size_counts", "genes_counts_ge10")
stopifnot(setequal(col_order, names(cd)))
cd <- cd[, col_order]
emit(c("", "-- final colData structure --", capture.output(str(cd, vec.len = 2))))


# ---- 8. (Optional) Salmon run metadata --------------------------------------
section("8. (Optional) Salmon run metadata from meta_info.json")

salmon_tab <- tryCatch({
  meta_json <- scan_files[tolower(basename(scan_files)) == "meta_info.json"]
  if (length(meta_json) == 0) {
    say("No meta_info.json files found under base_dir - skipped.")
    NULL
  } else if (!requireNamespace("jsonlite", quietly = TRUE)) {
    say("Package jsonlite is not installed - skipped.")
    NULL
  } else {
    get1 <- function(x, k) if (is.null(x[[k]])) NA_character_ else paste(x[[k]], collapse = ",")
    st <- do.call(rbind, lapply(meta_json, function(f) {
      mi <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
      if (is.null(mi)) return(NULL)
      sdir <- basename(dirname(dirname(f)))            # .../<sample>/aux_info/meta_info.json
      sid  <- c(sdir, sub("[_.].*$", "", sdir))
      data.frame(dataset           = if (any(sid %in% colnames(par_txi$counts))) lab_P
                 else if (any(sid %in% colnames(tse_txi$counts))) lab_T
                 else if (grepl("GSE242195|paryani", f, ignore.case = TRUE)) lab_P
                 else if (grepl("GSE79666|tseng", f, ignore.case = TRUE)) lab_T
                 else "other/unknown",
                 sample_dir        = sdir,
                 salmon_version    = get1(mi, "salmon_version"),
                 library_types     = get1(mi, "library_types"),
                 gc_bias_correct   = get1(mi, "gc_bias_correct"),
                 seq_bias_correct  = get1(mi, "seq_bias_correct"),
                 num_decoy_targets = get1(mi, "num_decoy_targets"),
                 index_seq_hash    = substr(get1(mi, "index_seq_hash"), 1, 12),
                 num_processed     = get1(mi, "num_processed"),
                 percent_mapped    = get1(mi, "percent_mapped"),
                 path              = f,
                 stringsAsFactors  = FALSE)
    }))
    if (!is.null(st) && nrow(st) > 0) {
      show_check("Salmon_run_summary_by_dataset", do.call(rbind, lapply(split(st, st$dataset), function(d) {
        pm <- suppressWarnings(as.numeric(d$percent_mapped))
        data.frame(dataset           = d$dataset[1],
                   n_runs            = nrow(d),
                   salmon_version    = paste(unique(d$salmon_version), collapse = " | "),
                   index_seq_hash    = paste(unique(d$index_seq_hash), collapse = " | "),
                   library_types     = paste(unique(d$library_types), collapse = " | "),
                   gc_bias           = paste(unique(d$gc_bias_correct), collapse = " | "),
                   decoys            = paste(unique(d$num_decoy_targets), collapse = " | "),
                   pct_mapped_median = if (all(is.na(pm))) NA else round(median(pm, na.rm = TRUE), 1),
                   pct_mapped_range  = if (all(is.na(pm))) NA else paste(round(range(pm, na.rm = TRUE), 1), collapse = " - "))
      })))
    }
    st
  }
}, error = function(e) {
  say("Salmon metadata section failed (not critical): ", conditionMessage(e))
  NULL
})


# ---- 9. Save outputs (never overwrite) --------------------------------------
section("9. Save outputs")

out <- c(merged_txi_rds  = safe_path(obj_dir, "merged_txi_Paryani76_Tseng14.rds"),
         colData_rds     = safe_path(obj_dir, "merged_colData_Paryani76_Tseng14.rds"),
         colData_csv     = safe_path(tab_dir, "merged_colData_Paryani76_Tseng14.csv"),
         gene_annot_rds  = safe_path(obj_dir, "gene_annotation_GENCODEv44.rds"),
         gene_annot_csv  = safe_path(tab_dir, "gene_annotation_GENCODEv44.csv"),
         sex_table_csv   = safe_path(tab_dir, "sex_markers_log2CPM_inferred_sex.csv"),
         donor_table_csv = safe_path(tab_dir, "donor_summary_Paryani_Tseng.csv"),
         checks_rds      = safe_path(obj_dir, "stage1_checks.rds"),
         sex_plot_pdf    = safe_path(fig_dir, "stage1_sex_inference_XIST_vs_Ymarkers.pdf"),
         sex_plot_png    = safe_path(fig_dir, "stage1_sex_inference_XIST_vs_Ymarkers.png"),
         log_txt         = safe_path(log_dir, "stage1_log.txt"))
if (!is.null(salmon_tab)) out["salmon_meta_csv"] <- safe_path(tab_dir, "salmon_run_metadata.csv")

saveRDS(merged_txi, out[["merged_txi_rds"]])
saveRDS(cd, out[["colData_rds"]])
write.csv(cd, out[["colData_csv"]], row.names = FALSE)
saveRDS(gene_annot, out[["gene_annot_rds"]])
write.csv(gene_annot, out[["gene_annot_csv"]], row.names = FALSE)
write.csv(sex_tab[, setdiff(names(sex_tab), "reported_for_plot")], out[["sex_table_csv"]], row.names = FALSE)
write.csv(donor_tab, out[["donor_table_csv"]], row.names = FALSE)
if (!is.null(salmon_tab)) write.csv(salmon_tab, out[["salmon_meta_csv"]], row.names = FALSE)
tryCatch({
  ggsave(out[["sex_plot_pdf"]], p_sex, width = 8, height = 6)
  ggsave(out[["sex_plot_png"]], p_sex, width = 8, height = 6, dpi = 300)
}, error = function(e) {
  say("Saving the sex plot failed: ", conditionMessage(e))
  flags <<- c(flags, "Sex plot could not be saved.")
})

section("FLAGS TO REVIEW")
if (length(flags) == 0) say("None.") else say(paste0(" - ", flags))

saveRDS(c(checks, list(flags = flags, input_paths = input_paths, outputs = out)), out[["checks_rds"]])
show_check("saved_files", data.frame(output = names(out), path = unname(out)))

writeLines(c(paste0("Stage 1 log - ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_lines,
             "", "-- sessionInfo --", capture.output(sessionInfo())),
           out[["log_txt"]])
cat("\nStage 1 finished. Log written to:", out[["log_txt"]], "\n")