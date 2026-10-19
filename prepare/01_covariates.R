# Participant-level covariates from the release tables: demographics, HDP exposure and
# pregnancy complications, BIBSNet total brain volume, and CSV copies of the maternal
# questionnaires read by the sensitivity analyses.

source("R/common.R")

num <- function(v) suppressWarnings(as.numeric(v))

D <- "sed_basic_demographics_"
demo <- read_release("demographics/sed_basic_demographics.tsv")
demo <- demo[order(demo$session_id), ]
first <- function(col) {
  v <- num(demo[[paste0(D, col)]])
  if (col == "sex") v[!v %in% 0:1] <- NA
  tapply(v, demo[[ID]], function(x) x[!is.na(x)][1])
}
person <- data.frame(names(first("sex")), first("mother_age_delivery"), first("gestational_age_delivery"),
                     first("rc_mother_education"), first("rc_mother_income"),
                     first("rc_mother_ethnoracial_aou_race_ethnicity"), first("sex"))
names(person) <- c(ID, "mother_age", GA, "education", "income", "race_eth", "child_sex")
# Recruitment site as recorded at the most recent visit.
visits <- read_release("demographics/par_visit_data.tsv")
visits <- visits[order(visits$session_id, decreasing = TRUE), ]
site <- visits[!duplicated(visits[[ID]]), c(ID, "par_visit_data_site")]
names(site)[2] <- "site"
write_input(merge(person, site, by = ID, all.x = TRUE), "covariates/hdp_missingness_person_level.csv")

pregnancy <- read_release("pregnancy_exposure/pex_bm_healthv2_preg.tsv")
pregnancy <- pregnancy[!duplicated(pregnancy[[ID]]), ]
items <- c("gestational_diabetes", "gestational_hypertension", "preeclampsia_eclampsia",
           "rare_complication_4", "rare_complication_5", "rare_complication_6", "chorioamnionitis",
           "placenta_previa", "oligohydramnios", "polyhydramnios", "fetal_growth_restriction")
exposure <- data.frame(pregnancy[[ID]], sapply(paste0("pex_bm_healthv2_preg__compl_001___", 1:11),
                                               function(x) num(pregnancy[[x]])))
names(exposure) <- c(ID, items)
other <- setdiff(items, c("gestational_hypertension", "preeclampsia_eclampsia"))
exposure$other_complication <- as.integer(rowSums(exposure[other] == 1, na.rm = TRUE) > 0)

chronic <- read_release("pregnancy_exposure/pex_bm_health_preg__chroncond.tsv")
chronic <- chronic[chronic$session_id == "ses-V01", ]
chronic <- chronic[!duplicated(chronic[[ID]]), c(ID, "pex_bm_health_preg__chroncond_001___1")]
names(chronic)[2] <- "chronic_hypertension"
chronic$chronic_hypertension <- num(chronic$chronic_hypertension)
exposure <- merge(exposure, chronic, by = ID, all.x = TRUE)

pregnancy_hdp <- exposure$gestational_hypertension == 1 | exposure$preeclampsia_eclampsia == 1
groups <- c("0 0" = "normotensive", "1 0" = "gestational", "1 1" = "superimposed", "0 1" = "chronic_only")
group <- groups[paste(as.integer(pregnancy_hdp), exposure$chronic_hypertension)]
stopifnot(table(factor(group, groups)) == c(1564, 306, 189, 116))
exposure[[HDP]] <- ifelse(group == "chronic_only", NA, as.integer(group != "normotensive"))
write_input(exposure[, c(ID, HDP, "chronic_hypertension", "gestational_hypertension",
                         "preeclampsia_eclampsia", "fetal_growth_restriction", "other_complication")],
            "covariates/v02_pregnancy_complications.tsv")

# Total brain volume: cerebral, cerebellar, subcortical and brainstem tissue, excluding
# ventricles and CSF. Missing or negative components, or an empty hemisphere, give NA.
bibsnet <- read_release("imaging/struct_bibsnet/img_bibsnet_space-T2w_desc-aseg_volumes.tsv")
prefix <- "img_bibsnet_space-T2w_desc-aseg_volumes_"
bibsnet <- bibsnet[bibsnet$session_id == "ses-V02", ]
tissue <- c("Cerebral-White-Matter", "Cerebral-Cortex", "Cerebellum-Cortex", "Cerebellum-White-Matter",
            "Thalamus-Proper*", "Caudate", "Putamen", "Pallidum", "Hippocampus", "Amygdala",
            "Accumbens-area", "VentralDC")
hemi <- lapply(c("Left", "Right"), function(s) sapply(bibsnet[paste0(prefix, s, "-", tissue)], num))
parts <- cbind(hemi[[1]], hemi[[2]], num(bibsnet[[paste0(prefix, "Brain-Stem")]]))
tbv <- rowSums(parts)
tbv[(rowSums(hemi[[1]]) <= 0 | rowSums(hemi[[2]]) <= 0 | tbv <= 0 | rowSums(parts < 0) > 0) %in% TRUE] <- NA
write_input(data.frame(participant_id = bibsnet[[ID]], total_brain_volume_mm3 = tbv),
            "covariates/v02_bibsnet_brain_volume.tsv")

stopifnot(file.copy(file.path(RELEASE, "imaging/qc_mriqc/img_mriqc_T2w.tsv"),
                    file.path(DATA, "covariates/img_mriqc_T2w.tsv"), overwrite = TRUE))

for (table in c("pex_bm_healthv2_inf", "pex_bm_healthv2_preg", "pex_bm_health_preg__healthhx",
                "pex_bm_health_preg__chroncond", "pex_bm_health_preg__meds", "pex_ch_tlfb")) {
  x <- read_release(paste0("pregnancy_exposure/", table, ".tsv"), colClasses = "character",
                    na.strings = character(0))
  write_input(x, paste0("phenotype/", table, ".csv"))
}
