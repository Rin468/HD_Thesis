# =============================================================================
# GSE242195: create a gene-level count bundle from Salmon quantifications
# =============================================================================
#
# Run this script once to rebuild the 76-sample tximport object and gene-count
# matrix. In later R sessions, load the saved bundle with readRDS() instead of
# rerunning tximport.

options(timeout = 900)

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tximport)
})

# 1. File locations ----------------------------------------------------------

project_dir <- "C:/!Virag.ORI/Bioinfo/szakdoga/Paryani"
quant_dir <- file.path(project_dir, "quant")

tx2gene_file <- file.path(project_dir, "Reference", "tx2gene_v44.csv")
gsm_srx_file <- file.path(
  project_dir,
  "Metadata",
  "gsm_title_srx_group_map.csv"
)
sra_file <- file.path(project_dir, "Metadata", "SraRunTable.csv")

bundle_file <- file.path(project_dir, "GSE242195_analysis_bundle.rds")
metadata_file <- file.path(
  project_dir,
  "sample_metadata_76_GSE242195.csv"
)

stopifnot(
  dir.exists(project_dir),
  dir.exists(quant_dir),
  file.exists(tx2gene_file),
  file.exists(gsm_srx_file),
  file.exists(sra_file)
)

# 2. Construct the 76-sample metadata table ---------------------------------

gse <- getGEO(GEO = "GSE242195", GSEMatrix = TRUE)
geo_metadata <- pData(phenoData(gse[[1]]))

# Match the spelling used in the local mapping table.
geo_metadata$title <- gsub(
  "H5608-Cing",
  "H5608_Cing",
  geo_metadata$title,
  fixed = TRUE
)

gsm_srx_map <- read.csv(gsm_srx_file, stringsAsFactors = FALSE)
sra_metadata <- read.csv(sra_file, stringsAsFactors = FALSE)

stopifnot(
  all(c("gsm", "title", "srx") %in% names(gsm_srx_map)),
  all(c("Experiment", "Run", "grade") %in% names(sra_metadata))
)

# Rename overlapping fields before joining so that title.x/title.y are not
# created and every column has an unambiguous meaning.
gsm_srx_map_clean <- gsm_srx_map %>%
  transmute(
    geo_accession = gsm,
    mapping_title = title,
    srx = srx
  )

metadata_76 <- geo_metadata %>%
  transmute(
    sample_title = title,
    geo_accession = geo_accession
  ) %>%
  filter(!grepl("^HMC3", sample_title)) %>%
  left_join(gsm_srx_map_clean, by = "geo_accession") %>%
  left_join(sra_metadata, by = c("srx" = "Experiment")) %>%
  dplyr::rename(srr_id = Run)

# Derive disease status from the neuropathological grade. The group column in
# the original local mapping file is not used because it labelled all samples
# as controls.
metadata_76 <- metadata_76 %>%
  mutate(
    grade = trimws(as.character(grade)),
    group = factor(
      if_else(grade == "Control", "CTRL", "HD"),
      levels = c("CTRL", "HD")
    )
  )

stopifnot(
  nrow(metadata_76) == 76L,
  !anyNA(metadata_76$srr_id),
  !anyNA(metadata_76$geo_accession),
  !anyNA(metadata_76$grade),
  !anyNA(metadata_76$group),
  !anyDuplicated(metadata_76$srr_id),
  !anyDuplicated(metadata_76$geo_accession),
  sum(metadata_76$group == "CTRL") == 20L,
  sum(metadata_76$group == "HD") == 56L
)

# 3. Locate and order the 76 Salmon quant.sf files --------------------------

all_files <- list.files(
  path = quant_dir,
  pattern = "^quant\\.sf$",
  recursive = TRUE,
  full.names = TRUE
)

names(all_files) <- sub(
  "_quant$",
  "",
  basename(dirname(all_files))
)

stopifnot(
  length(all_files) > 0L,
  !anyDuplicated(names(all_files)),
  all(metadata_76$srr_id %in% names(all_files))
)

# Subset and order files according to the metadata, then use GSM accessions as
# the final sample names.
files_76 <- all_files[metadata_76$srr_id]

stopifnot(
  identical(names(files_76), metadata_76$srr_id),
  all(file.exists(files_76))
)

names(files_76) <- metadata_76$geo_accession

# 4. Load the GENCODE v44 transcript-to-gene mapping ------------------------

tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)

stopifnot(
  all(c("transcript_id", "gene_id") %in% names(tx2gene))
)

tx2gene <- tx2gene %>%
  select(transcript_id, gene_id) %>%
  distinct()

stopifnot(
  nrow(tx2gene) > 0L,
  !anyNA(tx2gene$transcript_id),
  !anyNA(tx2gene$gene_id)
)

# 5. Aggregate transcript estimates to genes with tximport -----------------

txi_76 <- tximport(
  files = files_76,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreAfterBar = TRUE
)

count_matrix_76 <- txi_76$counts

stopifnot(
  identical(dim(count_matrix_76), c(62266L, 76L)),
  identical(colnames(count_matrix_76), metadata_76$geo_accession),
  identical(dim(txi_76$length), dim(count_matrix_76)),
  identical(dim(txi_76$abundance), dim(count_matrix_76)),
  txi_76$countsFromAbundance == "no",
  !anyNA(count_matrix_76),
  all(is.finite(count_matrix_76)),
  all(count_matrix_76 >= 0),
  !anyDuplicated(rownames(count_matrix_76)),
  !anyDuplicated(colnames(count_matrix_76))
)

# 6. Save one complete, reproducible analysis bundle ------------------------

analysis_bundle <- list(
  txi = txi_76,
  count_matrix = count_matrix_76,
  sample_metadata = metadata_76,
  provenance = list(
    geo_accession = "GSE242195",
    annotation = "GENCODE v44",
    quantifier = "Salmon",
    counts_from_abundance = txi_76$countsFromAbundance
  )
)

saveRDS(analysis_bundle, bundle_file, compress = "xz")
write.csv(metadata_76, metadata_file, row.names = FALSE)

# 7. Reload the saved file and verify it ------------------------------------

check_bundle <- readRDS(bundle_file)

stopifnot(
  identical(check_bundle$count_matrix, check_bundle$txi$counts),
  identical(dim(check_bundle$count_matrix), c(62266L, 76L)),
  identical(
    colnames(check_bundle$count_matrix),
    check_bundle$sample_metadata$geo_accession
  ),
  sum(check_bundle$sample_metadata$group == "CTRL") == 20L,
  sum(check_bundle$sample_metadata$group == "HD") == 56L,
  !anyNA(check_bundle$sample_metadata$group)
)

cat("\nVALIDATION PASSED\n")
cat("Genes:", nrow(check_bundle$count_matrix), "\n")
cat("Samples:", ncol(check_bundle$count_matrix), "\n")
cat("Groups:\n")
print(table(check_bundle$sample_metadata$group, useNA = "always"))

# In a later R session, load the finished analysis with:
# analysis_bundle <- readRDS(bundle_file)
# txi_76 <- analysis_bundle$txi
# count_matrix_76 <- analysis_bundle$count_matrix
# metadata_76 <- analysis_bundle$sample_metadata

#for creating the multiqc with only the 76 samples
writeLines(metadata_76$srr_id, "C:/!Virag.ORI/Bioinfo/szakdoga/Paryani/Metadata/srr_list_76.txt")
