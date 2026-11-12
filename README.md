# Huntington's disease bulk RNA-seq: comparison of two studies

R analysis scripts

The work compares gene expression between Huntington's disease (HD) and control (CTRL)
post-mortem brain samples, using two public datasets analysed separately and then compared.

Raw data are not stored here. Download them from GEO: https://www.ncbi.nlm.nih.gov/geo/

# Paryani et al. (2024), GSE242195: from raw reads to gene counts

Steps 1–6 are SLURM scripts run on the lab's computing cluster; step 7 is an R script.

| Order | Script | What it does |
|---|---|---|
| 1 | `01_download_fastq2.sh` | Downloads each run in `SRR_Acc_List.txt` from SRA (prefetch and fasterq-dump, SRA Toolkit 3.4.1) and compresses the FASTQ files |
| 2 | `02_classify_safe.sh` | Checks that each downloaded file contains real base-quality values and lists the complete runs |
| 3 | `03_build_qc_set.sh` | Collects the complete runs for quality control; if a run was downloaded from both NCBI and ENA, keeps the NCBI files |
| 4 | `04_run_qc.sh` | Read quality check with FastQC, summarised with MultiQC |
| 5 | `05_trim_fastp.sh` | Trimming with fastp: automatic adapter detection for paired-end reads, quality threshold Q20, minimum read length 25 |
| 6 | `06_run_salmon.sh` | Quantification with Salmon against the GENCODE release 44 index: automatic library type, no bias correction |
| 7 | `07_tximport_counts.R` | Imports the Salmon results of the 76 bulk tissue samples and sums transcripts to genes with tximport, using `reference/tx2gene_v44.csv` |

## Samples

`SRR_Acc_List.txt` lists all 325 sequencing runs of GSE242195, and all of them were processed by steps 1–6:
- 76 bulk tissue runs (35 cingulate cortex, 25 caudate nucleus, 16 nucleus accumbens), analysed in the thesis
- 189 runs labelled "TRANSCRIPTOMIC SINGLE CELL" in SRA, not used
- 60 runs from the HMC3 microglia cell line, not used

Step 7 selects the 76 bulk samples (listed in `metadata/Paryani_2024/srr_list_76.txt`) using
`metadata/Paryani_2024/SraRunTable.csv` and `gsm_title_srx_group_map.csv`. Disease status is taken
from the `grade` column of `SraRunTable.csv`; the `group` column of `gsm_title_srx_group_map.csv` is
incorrect (all samples labelled CTRL) and is not used.


## Source files (not stored here)
Download from GENCODE release 44: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
- `gencode.v44.transcripts.fa.gz`
- `gencode.v44.annotation.gtf.gz`
- 
## Notes

- 28 runs that did not finish in the first Salmon job were rerun with identical settings.
- The files of one run (SRR25877074) were compressed manually after download.
- Paths are those of the original server and must be changed to rerun the scripts elsewhere.
  The July scripts use `~/Paryani_2024`, the later ones `~/arcticles/Paryani_2024`
  (the project folder was moved).
