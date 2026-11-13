   # Metadata for Lin et al. (2016), GSE79666

   The server folder for this dataset is named `Tseng_2016`.

   - `SRR_Acc_List.txt` – the 14 SRA runs (one per sample, paired-end)
   - `SraRunTable.txt` – ENA file report
   - `gsm_title_srx_map.csv` – GEO sample IDs, sample names and SRA experiment IDs
   - `Tseng_sample_metadata.csv` – sample table written by `scripts/Lin_2016/05_tximport_counts.R` (7 HD, 7 control)

- `run_experiment_sample_ENA.tsv` – run (SRR), experiment (SRX) and GEO sample (GSM) for each run, retrieved from the ENA API; confirms the sample order used in `05_tximport_counts.R`
