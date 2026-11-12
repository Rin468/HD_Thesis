# Paryani et al. (2024), GSE242195: read processing on the server

SLURM scripts run on the lab's computing cluster, in this order:

`01_download_fastq.sh` | Downloads each run in `SRR_Acc_List.txt` from SRA (prefetch and fasterq-dump, SRA Toolkit 3.4.1) and compresses the FASTQ files |
`02_classify_downloads.sh` | Checks that each downloaded file contains real base-quality values and lists the complete runs |
`03_build_qc_set.sh` | Collects the complete runs for quality control; if a run was downloaded from both NCBI and ENA, keeps the NCBI files |
`04_run_qc.sh` | Read quality check with FastQC, summarised with MultiQC |
`05_trim_fastp.sh` | Trimming with fastp: automatic adapter detection for paired-end reads, quality threshold Q20, minimum read length 25 |
`06_run_salmon.sh` | Quantification with Salmon against a GENCODE release 44 index: automatic library type, no bias correction |

`SRR_Acc_List.txt` lists the 325 sequencing runs that were processed. Only the 76 bulk tissue
samples used in the thesis are imported into R and analysed.

- Paths are those of the original server and must be changed to rerun the scripts elsewhere.
