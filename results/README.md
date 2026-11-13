# Result tables

Copied from the final analysis runs (folder names in brackets). Contrast codes used in file
names and columns: A = striatum HD vs control, B_P = cingulate cortex HD vs control,
B_T = motor cortex HD vs control (Lin et al.), C = control striatum vs control cingulate,
D = HD striatum vs HD cingulate, E_P = difference between the striatal and cingulate HD effects.

| Folder | Content | Thesis |
|---|---|---|
| `01_samples_and_annotation` | Final sample table, tissue corrections, donor structure, independent check of the corrections (`bystudy_20261015_151051`, `annotation_check_20261101_155123`) | 3.1, 3.2, Table 1 |
| `02_main_DESeq2` | All genes and contrasts of the main models (`MASTER_DESeq2_by_study_compact.csv`), DEG counts, model designs, dispersion comparison, agreement between regions (`bystudy_20261015_151051`) | 3.3, Table 2 |
| `03_repeated_donors` | Null calibration (20 splits per type) and high-confidence genes from DESeq2 and donor-aware limma-voom (`nullcal_20261031_225350`, `voomexact_20261101_153230`) | 3.4 |
| `04_region_dependent` | Pooled cortex, region-dependent and robust genes, heterogeneity, cross-study check, core genes (`pooled_cortex_20261101_123913`, `voomexact_20261101_153230`, `heterogeneity_20261110_170130`, `regional_crossstudy_20261031_221615`, `coregenes_20261111_084547`) | 3.5, Table 3 |
| `05_cell_composition` | Marker scores per sample and HD vs control tests (`composition_20261101_003043`) | 3.6 |
| `06_GO_enrichment` | GO results per region and direction (`GO_20261017_223957`) and for the region-dependent genes (`GO_regional_20261101_163253`) | 3.6, Figures 6, 7 |
| `07_robustness` | Alternative analyses compared with the main result (`sensitivity_20261017_163558`, `S9followup_20261017_171514`, `results/sensitivity_age_GSE79666`, `results/sensitivity_sex_GSE79666`, `composition_20261101_003043`) | 3.6, Table 4 |
| `session_info` | R and package versions recorded by the main runs | 3.7 |
