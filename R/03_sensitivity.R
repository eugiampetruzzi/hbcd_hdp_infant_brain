# Sensitivity analyses of the global cortical thickness difference (Table 4). Each check is
# re-harmonized with the covariates of its own analysis model.

source("R/common.R")

m <- "Structural MRI"
Y <- "global_area_weighted_thickness"
QC <- MODALITIES[[m]]$qc
ys <- names(MODALITIES[[m]]$outcomes)

frames <- load_modalities()
final <- read_table("analytic_sample.csv")

# Birth weight for GA: sex-specific residual of caregiver-reported weight given s(GA).
infant <- read_input("phenotype/pex_bm_healthv2_inf.csv")
infant <- infant[!duplicated(infant[[ID]]), ]
infant$bw_g <- (as.numeric(infant$pex_bm_healthv2_inf_001__02) *
  16 +
  as.numeric(infant$pex_bm_healthv2_inf_001__01)) *
  28.3495
infant$bw_g[!(infant$bw_g > 400 & infant$bw_g < 6000)] <- NA
infant$nicu <- suppressWarnings(as.numeric(infant$pex_bm_healthv2_inf_005))
infant$nicu[!infant$nicu %in% c(0, 1)] <- NA
births <- do.call(rbind, lapply(frames, function(x) x[, c(ID, GA, "sex_m")]))
births <- merge(
  births[!duplicated(births[[ID]]), ],
  infant[, c(ID, "bw_g")],
  by = ID
)
births <- births[complete.cases(births) & births[[ID]] %in% final[[ID]], ]
bw_fit <- gam(
  as.formula(paste0("bw_g ~ sex_m + s(", GA, ", by = factor(sex_m), k = 5)")),
  data = births
)
births$bw_ga_z <- residuals(bw_fit) / sd(residuals(bw_fit))

complications <- read_input("covariates/v02_pregnancy_complications.tsv")
complications$fgr <- complications$fetal_growth_restriction
volume <- read_input("covariates/v02_bibsnet_brain_volume.tsv")
volume$log_tbv <- log(volume$total_brain_volume_mm3)
person <- read_input("covariates/hdp_missingness_person_level.csv")
person <- person[
  !duplicated(person[[ID]]),
  c(ID, "mother_age", "education", "income", "race_eth")
]

num <- function(v) suppressWarnings(as.numeric(v))
binary <- function(v) ifelse(num(v) %in% c(0, 1), num(v), NA)
first_row <- function(x) x[!duplicated(x[[ID]]), ]

# Any reported prenatal use of each substance on the Timeline Followback (Visits 1 and 2).
# TLFB weeks 1-2 precede the last menstrual period, so only weeks 3-9 count as prenatal.
tlfb <- read_input("phenotype/pex_ch_tlfb.csv")
substance <- data.frame(unique(tlfb[[ID]]))
names(substance) <- ID
for (s in c("nic", "alc", "thc", "opd")) {
  weeks <- grep(
    paste0("^pex_ch_tlfb_", s, "_wk_0[3-9]"),
    names(tlfb),
    value = TRUE
  )
  used <- tapply(
    rowSums(sapply(tlfb[, weeks], num) > 0, na.rm = TRUE) > 0,
    tlfb[[ID]],
    any
  )
  substance[[paste0("tlfb_", s)]] <- as.numeric(used[substance[[ID]]])
}
substance$any_substance <- as.numeric(rowSums(substance[, -1]) > 0)

history <- first_row(read_input("phenotype/pex_bm_health_preg__healthhx.csv"))
# Live births (preghx_002) is skipped for first pregnancies (preghx_001 = 1): nulliparous.
pregnancies <- num(history$pex_bm_health_preg__healthhx__preghx_001)
live_births <- num(history$pex_bm_health_preg__healthhx__preghx_002)
live_births[pregnancies %in% 1 & is.na(live_births)] <- 0
maternal <- data.frame(
  history[[ID]],
  num(history$pex_bm_health_preg__healthhx_011),
  ifelse(live_births %in% 0:20, as.numeric(live_births == 0), NA),
  binary(history$pex_bm_health_preg__healthhx__exp__pnv_012)
)
names(maternal) <- c(ID, "bmi", "nulliparous", "aspirin_v1")
maternal$bmi[!(maternal$bmi > 13 & maternal$bmi < 70)] <- NA
chronic <- first_row(read_input("phenotype/pex_bm_health_preg__chroncond.csv"))
chronic$pregest_diabetes <- as.numeric(
  binary(chronic$pex_bm_health_preg__chroncond_001___3) %in%
    1 |
    binary(chronic$pex_bm_health_preg__chroncond_001___4) %in% 1
)
pregnancy <- first_row(read_input("phenotype/pex_bm_healthv2_preg.csv"))
labor <- num(pregnancy$pex_bm_healthv2_preg__outcome_001)
delivery <- data.frame(
  pregnancy[[ID]],
  binary(pregnancy$pex_bm_healthv2_preg__compl_001___1),
  binary(pregnancy$pex_bm_healthv2_preg__exp__pnv_012),
  ifelse(labor %in% 0:2, as.numeric(labor == 1), NA),
  binary(pregnancy$pex_bm_healthv2_preg__outcome_002),
  binary(pregnancy$pex_bm_healthv2_preg__outcome_002__01___13)
)
names(delivery) <- c(
  ID,
  "gest_diabetes",
  "aspirin_v2",
  "induced",
  "cesarean",
  "multiple_gestation"
)

# Antihypertensive use: RxNorm codes recorded at Visit 1 or 2 that Visit 1 names as an antihypertensive.
meds <- first_row(read_input("phenotype/pex_bm_health_preg__meds.csv"))
named <- grep("_name$", names(meds), value = TRUE)
pairs <- do.call(
  rbind,
  lapply(named, function(n) {
    data.frame(
      code = as.character(meds[[sub("_name$", "", n)]]),
      name = as.character(meds[[n]])
    )
  })
)
antihypertensive_names <- "labetalol|nifedipine|methyldopa|hydralazine|amlodipine|metoprolol|procardia|trandate|aldomet"
codes <- unique(c(
  pairs$code[grepl(antihypertensive_names, pairs$name, ignore.case = TRUE)],
  "6185",
  "7417",
  "6876",
  "5470",
  "17767",
  "6918"
))
any_code <- function(x, cols) {
  apply(x[, cols, drop = FALSE], 1, function(r) {
    any(sub("\\.0$", "", r) %in% codes)
  })
}
v1_codes <- grep(
  "^pex_bm_health_preg__meds_0[0-9]{2}__01$",
  names(meds),
  value = TRUE
)
v2_codes <- grep(
  "^pex_bm_healthv2_preg__meds_0(4[7-9]|[5-9][0-9])__01$",
  names(pregnancy),
  value = TRUE
)
treated <- merge(
  data.frame(participant_id = meds[[ID]], v1 = any_code(meds, v1_codes)),
  data.frame(
    participant_id = pregnancy[[ID]],
    v2 = any_code(pregnancy, v2_codes)
  ),
  by = ID,
  all = TRUE
)
treated$antihypertensive <- as.numeric(
  treated$v1 %in% TRUE | treated$v2 %in% TRUE
)

# T2-weighted gray/white contrast from MRIQC.
mriqc <- read_input("covariates/img_mriqc_T2w.tsv")
mriqc <- first_row(mriqc[mriqc$session_id == "ses-V02", ])
gm <- mriqc$img_mriqc_T2w_summary_gm_median
wm <- mriqc$img_mriqc_T2w_summary_wm_median
contrast <- data.frame(
  mriqc[[ID]],
  as.numeric(scale((wm - gm) / ((wm + gm) / 2))),
  as.numeric(scale(mriqc$img_mriqc_T2w_cjv)),
  as.numeric(scale(mriqc$img_mriqc_T2w_cnr))
)
names(contrast) <- c(ID, "gw_contrast", "cjv", "cnr")

d0 <- modality_sample(frames, m, extra = "chronic_hypertension")
d0 <- d0[d0[[ID]] %in% final[[ID]][final$in_structural], ]
for (x in list(
  births[, c(ID, "bw_ga_z")],
  infant[, c(ID, "nicu")],
  complications[, c(ID, "other_complication", "preeclampsia_eclampsia", "fgr")],
  volume[, c(ID, "log_tbv")],
  person,
  substance,
  maternal,
  chronic[, c(ID, "pregest_diabetes")],
  delivery,
  treated[, c(ID, "antihypertensive")],
  contrast
)) {
  d0 <- merge(d0, x, by = ID, all.x = TRUE)
}
d0$diabetes <- as.numeric(d0$pregest_diabetes %in% 1 | d0$gest_diabetes %in% 1)
d0$aspirin <- ifelse(
  is.na(d0$aspirin_v1) & is.na(d0$aspirin_v2),
  NA,
  as.numeric(d0$aspirin_v1 %in% 1 | d0$aspirin_v2 %in% 1)
)
d0$hdp_x_male <- d0[[HDP]] * d0$sex_m
for (v in c("race_eth", "education", "income")) {
  dummies <- model.matrix(~ factor(d0[[v]]))[, -1, drop = FALSE]
  colnames(dummies) <- paste0(v, "_", seq_len(ncol(dummies)))
  d0[rownames(dummies), colnames(dummies)] <- dummies
}
dummies_of <- function(v) {
  grep(paste0("^", v, "_[0-9]+$"), names(d0), value = TRUE)
}

# Re-harmonize `d` preserving `extra`, then fit the primary model with `extra` added.
check <- function(
  group,
  label,
  d = d0,
  extra = character(0),
  model_extra = character(0),
  term = HDP,
  harmonized = TRUE,
  qc = QC
) {
  needed <- all.vars(reformulate(c("1", extra, model_extra)))
  d <- d[complete.cases(d[, needed, drop = FALSE]), ]
  d <- drop_small_batches(d)
  if (harmonized) {
    d <- harmonize(d, ys, harmonization_rhs(extra))
  }
  rhs <- paste(
    c(term, "sex_m", qc, extra, model_extra, smooth(PMA)),
    collapse = " + "
  )
  cbind(group = group, analysis = label, fit_term(d, Y, rhs, term))
}

race <- c("mother_age", dummies_of("race_eth"))
socio <- c(race, dummies_of("education"), dummies_of("income"))
perinatal <- c("bw_ga_z", "fgr", "nicu")
substances <- c("tlfb_nic", "tlfb_alc", "tlfb_thc", "tlfb_opd")
maternal_health <- c("bmi", "nulliparous", "diabetes")
best_qc <- d0[[QC]] >= median(d0[[QC]])

rows <- list(
  check("Primary model", "HDP vs normotensive, matched postmenstrual age"),
  check(
    "Imaging-related factors",
    "Scanner fixed effects",
    model_extra = "factor(batch)"
  ),
  check(
    "Imaging-related factors",
    "Site fixed effects",
    model_extra = "factor(site)"
  ),
  check(
    "Imaging-related factors",
    "No harmonization, scanner fixed effects",
    model_extra = "factor(batch)",
    harmonized = FALSE
  ),
  check(
    "Imaging-related factors",
    "Image quality as a smooth term",
    qc = smooth(QC, 5)
  ),
  check(
    "Imaging-related factors",
    "Highest-quality half of scans",
    d = d0[best_qc, ]
  ),
  check(
    "Imaging-related factors",
    "Gray/white contrast (T2-weighted contrast, CJV, CNR)",
    extra = c("gw_contrast", "cjv", "cnr")
  ),
  check("Imaging-related factors", "Total brain volume", extra = "log_tbv"),
  check(
    "Perinatal factors",
    "Birth weight for gestational age",
    extra = "bw_ga_z"
  ),
  check("Perinatal factors", "Fetal growth restriction", extra = "fgr"),
  check("Perinatal factors", "NICU admission", extra = "nicu"),
  check(
    "Perinatal factors",
    "All perinatal factors and brain volume",
    extra = c(perinatal, "log_tbv")
  ),
  check(
    "Perinatal factors",
    "Excluding NICU admission or growth restriction",
    d = d0[d0$nicu %in% 0 & d0$fgr %in% 0, ]
  ),
  check(
    "Perinatal factors",
    "Excluding births < 32 weeks",
    d = d0[d0[[GA]] >= 32, ]
  ),
  check(
    "Perinatal factors",
    "Excluding births < 37 weeks",
    d = d0[d0[[GA]] >= 37, ]
  ),
  check(
    "Perinatal factors",
    "Excluding multiple gestations",
    d = d0[!(d0$multiple_gestation %in% 1), ]
  ),
  check(
    "Perinatal factors",
    "Induction of labor and cesarean delivery",
    extra = c("induced", "cesarean")
  ),
  check(
    "Co-occurring pregnancy complications",
    "Adjusted for co-occurring complications",
    extra = "other_complication"
  ),
  check(
    "Co-occurring pregnancy complications",
    "Excluding co-occurring complications",
    d = d0[d0$other_complication %in% 0, ]
  ),
  check(
    "Maternal health and substance use",
    "Any prenatal substance use",
    extra = "any_substance"
  ),
  check(
    "Maternal health and substance use",
    "Nicotine, alcohol, cannabis, and opioid use",
    extra = substances
  ),
  check(
    "Maternal health and substance use",
    "Excluding any prenatal substance use",
    d = d0[d0$any_substance %in% 0, ]
  ),
  check(
    "Maternal health and substance use",
    "Pre-pregnancy body mass index",
    extra = "bmi"
  ),
  check(
    "Maternal health and substance use",
    "Nulliparity",
    extra = "nulliparous"
  ),
  check(
    "Maternal health and substance use",
    "Diabetes (pre-existing or gestational)",
    extra = "diabetes"
  ),
  check(
    "Maternal health and substance use",
    "Low-dose aspirin",
    extra = "aspirin"
  ),
  check(
    "Sociodemographic factors",
    "Maternal age and race and ethnicity",
    extra = race
  ),
  check(
    "Sociodemographic factors",
    "+ Maternal education",
    extra = c(race, dummies_of("education"))
  ),
  check(
    "Sociodemographic factors",
    "+ Household income",
    extra = c(race, dummies_of("income"))
  ),
  check(
    "Combined adjustment",
    "All sociodemographic, perinatal, maternal health, and substance-use factors and brain volume",
    extra = c(
      socio,
      "other_complication",
      perinatal,
      "induced",
      "cesarean",
      "log_tbv",
      substances,
      maternal_health
    )
  ),
  check(
    "Effect modification",
    "HDP × infant sex",
    extra = "hdp_x_male",
    model_extra = HDP,
    term = "hdp_x_male"
  )
)

# Within HDP: subtype contrasts, harmonized on the full sample with HDP preserved.
within_hdp <- function(label, exposure) {
  d <- d0[complete.cases(d0[, exposure, drop = FALSE]), ]
  d$exposure <- d[[HDP]] * d[[exposure]]
  h <- harmonize(
    drop_small_batches(d),
    ys,
    paste("exposure", harmonization_rhs(), sep = " + ")
  )
  h <- h[h[[HDP]] == 1, ]
  out <- fit_term(
    h,
    Y,
    paste("exposure + sex_m", QC, smooth(PMA), sep = " + "),
    "exposure"
  )
  out$n_hdp <- sum(h$exposure)
  cbind(group = "Clinical presentation (within HDP)", analysis = label, out)
}
rows <- c(
  rows,
  list(
    within_hdp(
      "Pre-eclampsia or eclampsia vs gestational hypertension",
      "preeclampsia_eclampsia"
    ),
    within_hdp(
      "With vs without history of hypertension",
      "chronic_hypertension"
    ),
    within_hdp(
      "With vs without antihypertensive medication",
      "antihypertensive"
    )
  )
)
sensitivity <- do.call(rbind, rows)
write_table(sensitivity, "sensitivity.csv")

loso <- do.call(
  rbind,
  lapply(sort(unique(d0$batch)), function(b) {
    cbind(
      scanner_left_out = b,
      check("Leave one scanner out", b, d = d0[d0$batch != b, ])[, -(1:2)]
    )
  })
)
write_table(loso, "leave_one_scanner_out.csv")

# HDP x scanner interaction, likelihood-ratio test on ML fits of the harmonized data.
h <- harmonize(d0, ys, harmonization_rhs())
h$z <- zscore(h[[Y]])
base <- paste("z ~", primary_rhs(m), "+ factor(batch)")
lr <- anova(
  gam(as.formula(base), data = h, method = "ML"),
  gam(
    as.formula(paste(base, "+", HDP, ": factor(batch)")),
    data = h,
    method = "ML"
  ),
  test = "Chisq"
)
write_table(
  data.frame(df = lr$Df[2], deviance = lr$Deviance[2], p = lr$`Pr(>Chi)`[2]),
  "scanner_interaction.csv"
)

print(
  sensitivity[, c("analysis", "n", "n_hdp", "beta", "ci_lo", "ci_hi", "p")],
  digits = 2
)
cat(sprintf(
  "Leave-one-scanner-out: %.2f to %.2f; HDP x scanner p = %.2f\n",
  min(loso$beta),
  max(loso$beta),
  lr$`Pr(>Chi)`[2]
))
