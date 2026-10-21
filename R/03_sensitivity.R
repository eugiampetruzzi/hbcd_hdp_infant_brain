# Sensitivity analyses of the global cortical thickness difference (Table 4). Each check is
# re-harmonized with the covariates of its own analysis model.

source("R/common.R")
source("R/sensitivity_covariates.R")

m <- "Structural MRI"
Y <- "global_area_weighted_thickness"
QC <- MODALITIES[[m]]$qc
ys <- names(MODALITIES[[m]]$outcomes)

frames <- load_modalities()
final <- read_table("analytic_sample.csv")

# Birth weight for GA: sex-specific residual of caregiver-reported weight given s(GA).
infant <- read_input("phenotype/pex_bm_healthv2_inf.csv")
infant <- infant[!duplicated(infant[[ID]]), ]
infant <- birth_weight_covariate(infant)
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
volume <- brain_volume_covariate(read_input(
  "covariates/v02_bibsnet_brain_volume.tsv"
))
person <- read_input("covariates/hdp_missingness_person_level.csv")
person <- person[
  !duplicated(person[[ID]]),
  c(ID, "mother_age", "education", "income", "race_eth")
]

num <- function(v) suppressWarnings(as.numeric(v))
binary <- function(v) ifelse(num(v) %in% c(0, 1), num(v), NA)
# 1 if any item is 1, NA if every item is unknown, otherwise 0.
any_yes <- function(...) {
  x <- cbind(...)
  ifelse(
    rowSums(x == 1, na.rm = TRUE) > 0,
    1,
    ifelse(rowSums(!is.na(x)) == 0, NA, 0)
  )
}
first_row <- function(x) x[!duplicated(x[[ID]]), ]

# Any reported prenatal use of each substance on the Timeline Followback (Visits 1 and 2).
# TLFB weeks 1-2 precede the last menstrual period, so only weeks 3-9 count as prenatal.
tlfb <- read_input("phenotype/pex_ch_tlfb.csv")
substance <- substance_exposure(tlfb)

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
chronic$pregest_diabetes <- any_yes(
  binary(chronic$pex_bm_health_preg__chroncond_001___3),
  binary(chronic$pex_bm_health_preg__chroncond_001___4)
)
pregnancy <- first_row(read_input("phenotype/pex_bm_healthv2_preg.csv"))
labor <- num(pregnancy$pex_bm_healthv2_preg__outcome_001)
delivery <- data.frame(
  pregnancy[[ID]],
  binary(pregnancy$pex_bm_healthv2_preg__compl_001___1),
  binary(pregnancy$pex_bm_healthv2_preg__exp__pnv_012),
  ifelse(labor %in% 0:2, as.numeric(labor == 1), NA),
  ifelse(labor %in% 0:2, as.numeric(labor == 2), NA),
  binary(pregnancy$pex_bm_healthv2_preg__outcome_002),
  binary(pregnancy$pex_bm_healthv2_preg__outcome_002__01___13)
)
# Labor onset: spontaneous (reference), induced, or no labor (planned cesarean).
names(delivery) <- c(
  ID,
  "gest_diabetes",
  "aspirin_v2",
  "induced",
  "no_labor",
  "cesarean",
  "multiple_gestation"
)

# Antihypertensive use: RxNorm codes recorded at Visit 1 or 2 that Visit 1 names as an antihypertensive.
meds <- first_row(read_input("phenotype/pex_bm_health_preg__meds.csv"))
treated <- medication_exposure(meds, pregnancy)

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
  treated[, c(ID, "antihypertensive")]
)) {
  d0 <- merge(d0, x, by = ID, all.x = TRUE)
}
d0$diabetes <- any_yes(d0$pregest_diabetes, d0$gest_diabetes)
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
labor_onset <- c("induced", "no_labor")
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
    "Labor onset (spontaneous, induced, no labor)",
    extra = labor_onset
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
      labor_onset,
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

# HDP x labor onset: likelihood-ratio test on ML fits, and the HDP difference within each onset.
d <- d0[complete.cases(d0[, labor_onset]), ]
d$hdp_induced <- d[[HDP]] * d$induced
d$hdp_no_labor <- d[[HDP]] * d$no_labor
interaction <- c("hdp_induced", "hdp_no_labor")
h <- harmonize(
  drop_small_batches(d),
  ys,
  harmonization_rhs(c(labor_onset, interaction))
)
h$z <- zscore(h[[Y]])
f_main <- paste("z ~", primary_rhs(m, extra = labor_onset))
f_int <- paste(f_main, "+", paste(interaction, collapse = " + "))
lr_labor <- anova(
  gam(as.formula(f_main), data = h, method = "ML"),
  gam(as.formula(f_int), data = h, method = "ML"),
  test = "Chisq"
)
g <- gam(as.formula(f_int), data = h, method = "REML")
b <- coef(g)
V <- vcov(g, sandwich = TRUE)
onset <- list(
  Spontaneous = character(0),
  Induced = "hdp_induced",
  "No labor" = "hdp_no_labor"
)
labor_rows <- do.call(
  rbind,
  lapply(names(onset), function(o) {
    w <- setNames(as.numeric(names(b) %in% c(HDP, onset[[o]])), names(b))
    est <- sum(w * b)
    se <- sqrt(drop(t(w) %*% V %*% w))
    in_onset <- if (o == "Spontaneous") {
      h$induced == 0 & h$no_labor == 0
    } else {
      h[[sub("hdp_", "", onset[[o]])]] == 1
    }
    data.frame(
      onset = o,
      n = sum(in_onset),
      n_hdp = sum(in_onset & h[[HDP]] == 1),
      beta = est,
      se = se,
      ci_lo = est - 1.96 * se,
      ci_hi = est + 1.96 * se,
      p = 2 * pnorm(-abs(est / se))
    )
  })
)
labor_rows$interaction_df <- lr_labor$Df[2]
labor_rows$interaction_p <- lr_labor$`Pr(>Chi)`[2]
write_table(labor_rows, "labor_onset.csv")

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
print(labor_rows, digits = 2)
