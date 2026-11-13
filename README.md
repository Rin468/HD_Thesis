# Comparative analysis of transcriptomic profiles in Huntington's disease and healthy brain tissues based on bulk RNA-seq data

Analysis code, metadata and result tables for my student research paper and BSc thesis.
Gene expression of Huntington's disease (HD) and control post-mortem brain samples is compared
between the striatum and the cortex, using two public datasets processed with one workflow.

## Data

| Dataset | Study | Brain regions | Samples analysed |
|---|---|---|---|
| GSE242195 | Paryani et al., 2024 | caudate nucleus, nucleus accumbens, cingulate cortex | 76 bulk RNA-seq samples |
| GSE79666 | Lin et al., 2016 | motor cortex (BA4) | 14 (7 HD, 7 control) |

Raw reads are not stored here; they are public in GEO/SRA (run lists in `metadata/`).
In file and object names the Lin et al. (2016) dataset is also called "Tseng".

## Methods (short)

Trimming with fastp, quantification with Salmon 2.4.1 against GENCODE release 44, gene-level
import with tximport, differential expression with DESeq2 (one model per study), donor-aware
limma-voom as a check for repeated donors, inverse-variance pooling of the two cortical regions,
GO enrichment with clusterProfiler.

## Folder structure

- `scripts/Paryani_2024/` – download, QC, trimming and Salmon on the server (steps 01–06) and gene counts in R (07)
- `scripts/Lin_2016/` – the same for Lin et al. (2016) (steps 01–05)
- `scripts/Merged_analysis/` – the main analysis in R; `additional/` holds earlier versions and checks
- `reference/` – GENCODE 44 transcript-to-gene table and Salmon index settings
- `metadata/` – SRA run lists and sample tables for both datasets
- `results/` – result tables of the final runs and package versions (`results/session_info/`)

Each folder has its own README with details.

## How to run

1. Server scripts in `scripts/Paryani_2024` and `scripts/Lin_2016`, in numerical order
2. The tximport scripts (`scripts/Paryani_2024/07_tximport_counts.R`, `scripts/Lin_2016/05_tximport_counts.R`)
3. `scripts/Merged_analysis`, in numerical order

Paths are those of the original server and computers and must be changed to rerun the scripts.

## Requirements

- Server: SRA Toolkit 3.4.1, FastQC, MultiQC, fastp, Salmon 2.4.1, SLURM
- R 4.5.3, Bioconductor 3.22; package versions in `results/session_info/`

