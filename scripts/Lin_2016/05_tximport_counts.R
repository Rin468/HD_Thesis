# ============================================================
# Tseng et al. 2016 - GSE79666
# Huntington's disease motor cortex RNA-seq analysis
# ============================================================


# ============================================================
# 1. Set the Tseng analysis directory
# ============================================================

setwd("/Users/ronalab/Virag/RStudio/Tseng")

getwd()

# ============================================================
# 2. Package setup
# ============================================================

# Install BiocManager only if it is missing.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install tximport only if it is missing.
if (!requireNamespace("tximport", quietly = TRUE)) {
  BiocManager::install(
    "tximport",
    update = FALSE,
    ask = FALSE
  )
}

# Install DESeq2 only if it is missing.
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install(
    "DESeq2",
    update = FALSE,
    ask = FALSE
  )
}

library(tximport)
library(DESeq2)

cat("Required packages loaded successfully!\n")

# ============================================================
# 3. Create basic Tseng sample metadata
# ============================================================

tseng_metadata <- data.frame(
  run_accession = paste0("SRR", 3306823:3306836),
  
  geo_accession = paste0("GSM", 2100621:2100634),
  
  sample_name = c(
    "HD-57",
    "HD-58",
    "HD-61",
    "HD-68",
    "HD-283",
    "HD-288",
    "HD-290",
    "CTRL-3",
    "CTRL-4",
    "CTRL-8",
    "CTRL-9",
    "CTRL-11",
    "CTRL-12",
    "CTRL-13"
  ),
  
  group = c(
    rep("HD", 7),
    rep("CTRL", 7)
  ),
  
  stringsAsFactors = FALSE
)

# Set CTRL as the reference group.
tseng_metadata$group <- relevel(
  factor(tseng_metadata$group),
  ref = "CTRL"
)

# Use SRR run accessions as metadata row names.
rownames(tseng_metadata) <- tseng_metadata$run_accession

tseng_metadata
table(tseng_metadata$group)

# ============================================================
# 4. Locate the 14 Salmon quantification files
# ============================================================

tseng_files <- file.path(
  "Outputs",
  "salmon_quant",
  paste0(tseng_metadata$run_accession, "_quant"),
  "quant.sf"
)

# Give every path its corresponding SRR sample name.
names(tseng_files) <- tseng_metadata$run_accession

# Display all paths.
tseng_files

# Check whether every file exists.
file.exists(tseng_files)

# Count the files R can find.
sum(file.exists(tseng_files))

# Stop immediately if even one file is missing.
stopifnot(all(file.exists(tseng_files)))

cat("All 14 Salmon quantification files were found!\n")

# ============================================================
# 5. Create transcript-to-gene mapping (Step 5 teaches R which transcripts 
# belong to which genes, so Step 6 can turn Salmon’s transcript
# counts from all 14 samples into gene counts.)
# ============================================================

# Read one representative Salmon quantification file.
example_quant <- read.delim(
  tseng_files[1],
  header = TRUE,
  check.names = FALSE
)

# Inspect the beginning of the file.
head(example_quant)

# Separate every pipe-delimited transcript identifier.
identifier_parts <- strsplit(
  example_quant$Name,
  split = "\\|"
)

# Confirm that every identifier has enough components.
identifier_lengths <- lengths(identifier_parts)

table(identifier_lengths)

stopifnot(all(identifier_lengths >= 6))

# Extract the transcript ID, gene ID, and gene symbol
# from the Salmon identifiers.
tseng_tx2gene <- data.frame(
  transcript_id = vapply(
    identifier_parts,
    function(x) x[1],
    character(1)
  ),
  
  gene_id = vapply(
    identifier_parts,
    function(x) x[2],
    character(1)
  ),
  
  gene_symbol = vapply(
    identifier_parts,
    function(x) x[6],
    character(1)
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# 6. Import all Salmon expression estimates at gene level
# ============================================================

# Read all 14 quant.sf files.
#
# tximport uses tseng_tx2gene to determine which transcripts
# belong to the same gene and summarizes them at gene level.
tseng_tximport <- tximport(
  files = tseng_files,
  type = "salmon",
  
  tx2gene = tseng_tx2gene[
    ,
    c("transcript_id", "gene_id")
  ],
  
  # Use the transcript ID before the first "|" in Salmon's Name field.
  ignoreAfterBar = TRUE
)

# ============================================================
# 7. Create a readable gene-level count matrix
# ============================================================

# Keep the gene ID and gene symbol from the transcript mapping.
tseng_gene_annotation <- unique(
  tseng_tx2gene[, c("gene_id", "gene_symbol")]
)

# Keep only one annotation row for each gene.
tseng_gene_annotation <- tseng_gene_annotation[
  !duplicated(tseng_gene_annotation$gene_id),
]

tseng_count_matrix <- as.data.frame(
  tseng_tximport$counts
)

# Turn the row names into a regular gene-ID column.
tseng_count_matrix$gene_id <- rownames(
  tseng_count_matrix
)

tseng_count_matrix$gene_symbol <- tseng_gene_annotation$gene_symbol[
  match(
    tseng_count_matrix$gene_id,
    tseng_gene_annotation$gene_id
  )
]

tseng_count_matrix <- tseng_count_matrix[
  ,
  c(
    "gene_id",
    "gene_symbol",
    tseng_metadata$run_accession
  )
]

#save results
dir.create(
  "Results",
  showWarnings = FALSE
)

write.csv(
  tseng_count_matrix,
  "Results/Tseng_gene_count_matrix.csv",
  row.names = FALSE
)

write.csv(
  tseng_metadata,
  "Results/Tseng_sample_metadata.csv",
  row.names = FALSE
)

saveRDS(
  tseng_tximport,
  "Results/Tseng_tximport.rds"
)

saveRDS(
  tseng_metadata,
  "Results/Tseng_metadata.rds"
)

saveRDS(
  tseng_tx2gene,
  "Results/Tseng_tx2gene.rds"
)

# ============================================================
# 8. Save the complete Tseng analysis bundle
# ============================================================

# Store all important objects together in one list.
tseng_analysis_bundle <- list(
  txi = tseng_tximport,
  count_matrix = tseng_count_matrix,
  metadata = tseng_metadata,
  tx2gene = tseng_tx2gene,
  
  study_information = list(
    study = "Tseng et al. 2016",
    accession = "GSE79666",
    tissue = "BA4 motor cortex",
    comparison = "HD versus CTRL",
    reference = "GENCODE v44",
    salmon_version = "2.4.1",
    library_type = "paired-end",
    input_reads = "fastp-trimmed reads"
  )
)

saveRDS(
  tseng_analysis_bundle,
  "Results/GSE79666_analysis_bundle.rds"
)

cat("Tseng analysis bundle saved successfully!\n")