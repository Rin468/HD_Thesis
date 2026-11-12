# =============================================================================
#  01b_revise_donor_and_sex.R
#  Stage 1b - Revise the merged colData from Stage 1:
#   (1) reconstruct donor IDs for the Paryani batch-4 ("OF") samples
#   (2) use expression-inferred sex as the sex covariate for ALL samples
#  Input : latest merged_colData_Paryani76_Tseng14*.rds (Stage 1 output)
#  Output: revised colData (RDS + CSV), OF donor map, revised donor summary
#  The merged txi object is NOT changed (same samples, same order).
#  Save as: Merged_analysis/scripts/01b_revise_donor_and_sex.R  -> run with Source
# =============================================================================

base_dir <- "C:/!Virag.ORI/Bioinfo/szakdoga"
obj_dir  <- file.path(base_dir, "Merged_analysis", "results", "objects")
tab_dir  <- file.path(base_dir, "Merged_analysis", "results", "tables")

safe_path <- function(dir, filename) {
  path <- file.path(dir, filename)
  if (!file.exists(path)) return(path)
  stem <- tools::file_path_sans_ext(filename); ext <- tools::file_ext(filename)
  i <- 2
  while (file.exists(file.path(dir, sprintf("%s_v%d.%s", stem, i, ext)))) i <- i + 1
  file.path(dir, sprintf("%s_v%d.%s", stem, i, ext))
}

# ---- 0. Load the most recent Stage 1 colData --------------------------------
cd_files <- list.files(obj_dir, pattern = "^merged_colData_Paryani76_Tseng14(_v[0-9]+)?\\.rds$",
                       full.names = TRUE)
if (length(cd_files) == 0) stop("Stage 1 colData not found in ", obj_dir)
cd_file <- cd_files[which.max(file.info(cd_files)$mtime)]
cat("Input colData:", cd_file, "\n")
cd <- readRDS(cd_file)
stopifnot(nrow(cd) == 90)

# ---- 1. Donor reconstruction for batch-4 "OF" samples -----------------------
# Rule: consecutive OF IDs with identical sex, grade, age bin and CAG genotype
# and no repeated region are treated as one donor. Explicit, documented map.
of_map <- data.frame(
  sample_label = sprintf("OF%03d", 19:52),
  donor_new = c("OF019", "OF020", "OF021",
                "OF022-023", "OF022-023",
                "OF024-025", "OF024-025",
                rep("OF026-028", 3),
                "OF029-030", "OF029-030",
                rep("OF031-033", 3),
                "OF034",
                rep("OF035-037", 3),
                "OF038-039", "OF038-039",
                rep("OF040-042", 3),
                rep("OF043-045", 3),
                rep("OF046-048", 3),
                "OF049-050", "OF049-050",
                "OF051-052", "OF051-052"),
  stringsAsFactors = FALSE)
of_map$evidence <- "consecutive IDs; identical sex, grade, age bin and CAG"
of_map$evidence[of_map$donor_new %in% c("OF019", "OF020", "OF021", "OF034")] <- "single sample"
of_map$evidence[of_map$donor_new %in% c("OF022-023", "OF029-030")] <-
  "identical sex, grade, age bin; CAG missing for >=1 sample (less certain)"
of_map$evidence[of_map$donor_new == "OF043-045"] <-
  "identical sex, grade, age bin; CAG 42-15 vs 41-15 (1-repeat difference)"
stopifnot(nrow(of_map) == 34,
          setequal(of_map$sample_label, cd$sample_label[grepl("^P_OF", cd$donor)]))

cd$donor_from_title <- cd$donor                       # keep the original ID
i <- match(cd$sample_label, of_map$sample_label)
cd$donor[!is.na(i)] <- paste0("P_", of_map$donor_new[i[!is.na(i)]])

# Every reconstructed donor must have one sex, grade and age bin, and no repeated region
chk <- do.call(rbind, lapply(split(cd, cd$donor), function(d) data.frame(
  donor = d$donor[1], n = nrow(d),
  n_sex = length(unique(d$sex_inferred)), n_grade = length(unique(d$grade)),
  n_age = length(unique(d$age_binned)),
  dup_region = anyDuplicated(as.character(d$region)) > 0)))
bad <- chk[chk$n_sex > 1 | chk$n_grade > 1 | chk$n_age > 1 | chk$dup_region, ]
if (nrow(bad) > 0) { print(bad); stop("Inconsistent reconstructed donor(s) - send me this output.") }
cat("Donor consistency check passed.\n")

# ---- 2. Sex covariate = expression-inferred sex for all samples -------------
stopifnot(all(cd$sex_inferred %in% c("female", "male")))
cd$sex_discordant_with_GEO <- !is.na(cd$sex_reported) & cd$sex_reported != cd$sex_inferred
cd$sex        <- factor(cd$sex_inferred, levels = c("female", "male"))
cd$sex_source <- "inferred_from_expression"
cat("\n-- Samples whose GEO sex label disagrees with expression --\n")
print(cd[cd$sex_discordant_with_GEO, c("sample_label", "donor", "sex_reported", "sex_inferred")])

# ---- 3. Revised checks ------------------------------------------------------
don_sum <- do.call(rbind, lapply(split(cd, cd$donor), function(d) data.frame(
  donor = d$donor[1], study = as.character(d$study[1]), condition = as.character(d$condition[1]),
  n_samples = nrow(d), regions = paste(as.character(sort(d$region)), collapse = "+"),
  batches = paste(unique(as.character(sort(d$batch_study))), collapse = "+"),
  sex = as.character(d$sex[1]), grade = d$grade[1], age_binned = d$age_binned[1],
  stringsAsFactors = FALSE)))
rownames(don_sum) <- NULL

cat("\n-- Donors per study x condition --\n")
print(table(study = don_sum$study, condition = don_sum$condition))
cat("\n-- Samples per donor --\n")
print(table(study_condition = paste(don_sum$study, don_sum$condition), n_samples = don_sum$n_samples))
cat("\nJuvenile-onset (grade J) donors:", sum(don_sum$grade %in% "J"), "\n")
cat("\n-- Sex (expression-based) x group x study, samples --\n")
print(ftable(table(study = cd$study, group = cd$group, sex = cd$sex)))

mm <- model.matrix(~ batch_study + sex + group, data = cd)
cat("\nDesign ~ batch_study + sex + group:", ncol(mm), "coefficients, rank", qr(mm)$rank, "\n")

# ---- 4. Save (never overwrite) ----------------------------------------------
out <- c(colData_rds = safe_path(obj_dir, "merged_colData_Paryani76_Tseng14_donorSexRevised.rds"),
         colData_csv = safe_path(tab_dir, "merged_colData_Paryani76_Tseng14_donorSexRevised.csv"),
         of_map_csv  = safe_path(tab_dir, "Paryani_batch4_OF_donor_reconstruction.csv"),
         donors_csv  = safe_path(tab_dir, "donor_summary_Paryani_Tseng_revised.csv"))
saveRDS(cd, out[["colData_rds"]])
write.csv(cd, out[["colData_csv"]], row.names = FALSE)
write.csv(of_map, out[["of_map_csv"]], row.names = FALSE)
write.csv(don_sum, out[["donors_csv"]], row.names = FALSE)
cat("\nSaved:\n"); print(data.frame(output = names(out), path = unname(out)))

