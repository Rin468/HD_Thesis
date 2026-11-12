# Main analysis (R)

Input: the gene-level tximport objects of the two studies, made by
`scripts/Paryani_2024/07_tximport_counts.R` and `scripts/Lin_2016/05_tximport_counts.R`.
In the data objects the Lin et al. (2016) dataset is labelled "Tseng".

Run the scripts in numerical order. Each script writes to a new time-stamped folder
(`final_runs/<name>_<date>_<time>/`) and does not overwrite earlier results.
The analysis of record is the main run `bystudy_20261015_151051`; later scripts read it.

| Script | What it does | Thesis section |
|---|---|---|
| `01_build_merged_dataset.R` | Combines the two tximport objects and sample tables; infers sex from XIST and Y-chromosome genes | 3.1, 3.2 |
| `01b_revise_donor_and_sex.R` | Reconstructs donor IDs of the Paryani batch-4 samples; uses expression-inferred sex as the sex covariate | 3.2 |
| `04_main_DESeq2_by_study.R` | Tissue check with marker scores and relabelling; main DESeq2 models (one per study); all contrasts, QC plots, striatum vs cingulate and cingulate vs motor cortex comparisons | 3.2, 3.3 (Figures 3, 4; Table 2) |
| `05_sensitivity_robustness_S1-S9.R` | Alternative analyses S1–S9 and limma-voom with and without donor blocking | 3.6 (Table 4) |
| `05c_S9_followup_cingulate.R` | Follow-up of S9 (sex-corrected samples) | 3.6 (Table 4) |
| `05d_S10_age_GSE79666.R` | Age and grade of the Lin samples from GEO; age-adjusted model (S10) | 3.6 (Table 4) |
| `05e_sex_check_GSE79666.R` | Lin model without the sex covariate | 3.6 (Table 4) |
| `06_null_calibration_Paryani.R` | Random splits of the striatal HD samples (20 per type) to test calibration for repeated donors | 3.4 |
| `07_thesis_requests_20261017.R` | Volcano plots and DEG counts per region; GO enrichment per region and direction | 3.6 (Figures 2, 7) |
| `08_regional_effect_crossstudy.R` | Striatum vs motor cortex across the two studies for the interaction genes | 3.5 |
| `09_high_confidence_DEGs.R` | Genes significant in both DESeq2 and donor-aware limma-voom | 3.4 |
| `10_composition_markers.R` | Cell-type marker scores per sample; neuronal score as covariate | 3.6 |
| `11_pooled_cortex.R` | Pooled cortex (inverse-variance meta-analysis of cingulate and motor cortex) and region-dependent HD effects | 3.5 (Table 3) |
| `11b_pooled_cortex_heterogeneity.R` | Heterogeneity between the two cortical HD effects | 3.5 |
| `11c_plot_striatum_vs_cortex_summary.R` | Summary figure of region-dependent HD effects | Figure 5 |
| `12_pooled_cortex_robustness.R` | Region-dependent genes repeated with donor-aware limma-voom and without batch 4 | 3.5 |
| `15_voom_exact_true_donors.R` | Exact donor-aware limma-voom results with the reconstructed donors; final high-confidence and robust gene lists | 3.4, 3.5 |
| `16_verify_annotation_errors.R` | Independent check of the sex and tissue corrections (whole-transcriptome similarity) | 3.2 |
| `18_GO_region_dependent.R` | GO enrichment of the region-dependent and robust genes | 3.6 |
| `18b_Figure_GO_robust_region_dependent_genes.R` | GO figure of the 704 robust genes | Figure 6 |
| `19_core_genes_all_three_regions.R` | Genes significant in all three regions (core genes) | 3.5 |
| `20_Figure_workflow_overview.R` | Workflow overview figure | Figure 1 |

`07_thesis_requests_20261017.R` was run from the console and recovered verbatim from the R history.
Script numbers 02, 03, 13, 14 and 17 belong to earlier or exploratory steps that are not part of the
final analysis. Paths are those of the original laptop and must be changed to rerun the scripts.
