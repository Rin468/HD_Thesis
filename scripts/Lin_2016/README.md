# Lin et al. (2016), GSE79666: from raw reads to gene counts

The server folder for this dataset is named `Tseng_2016`, and the data objects use the label "Tseng".
Steps 1–4 are SLURM scripts run on the lab's computing cluster; step 5 is an R script.

| Order | Script | What it does |
|---|---|---|
| 1 | `01_sra_download.sh` | Downloads each run in `metadata/Lin_2016/SRR_Acc_List.txt` from SRA (prefetch and fasterq-dump, SRA Toolkit 3.4.1) and compresses the FASTQ files |
| 2 | `02_run_qc.sh` | Read quality check with FastQC, summarised with MultiQC |
| 3 | `03_trim_fastp.sh` | Trimming with fastp: automatic adapter detection for paired-end reads, quality threshold Q20, minimum read length 25 |
| 4 | `04_run_salmon.sh` | Quantification with Salmon against the GENCODE release 44 index: automatic library type, no bias correction |
| 5 | `05_tximport_counts.R` | R script (run on the Mac): defines the 14 samples (7 HD, 7 control) and sums Salmon transcript estimates to genes with tximport |

`SRR_Acc_List.txt` lists the 14 sequencing runs (one per sample, paired-end).

The fastp and Salmon settings and the Salmon index are identical to those used for
Paryani et al. (2024); see `scripts/Paryani_2024` and `reference/`.

Paths are those of the original server and computer and must be changed to rerun the scripts elsewhere.
