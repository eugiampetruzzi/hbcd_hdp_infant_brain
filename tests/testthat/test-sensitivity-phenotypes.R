test_that("documented pounds and ounces yield the correct participant birth weights", {
  # https://docs.hbcdstudy.org/latest/instruments/pregexp/pex/
  # __01 is ounces; __02 is pounds. Unknown weight is not a zero-weight infant.
  e <- local_sensitivity_helpers()
  d <- release_frame(
    "infant_health",
    data.frame(
      participant_id = c("seven-half", "six-three-quarter", "unknown"),
      session_id = "ses-V02",
      pex_bm_healthv2_inf_001__01 = c("8", "12", NA),
      pex_bm_healthv2_inf_001__02 = c("7", "6", NA)
    )
  )
  got <- e$birth_weight_covariate(d[3:1, ])
  expect_equal(
    got$bw_g[match(d$participant_id, got$participant_id)],
    c(3401.94, 3061.746, NA)
  )
})

test_that("explicit unknown or declined birth weight cannot become a measured weight", {
  e <- local_sensitivity_helpers()
  d <- release_frame(
    "infant_health",
    data.frame(
      participant_id = c("known", "unknown", "declined"),
      session_id = "ses-V02",
      pex_bm_healthv2_inf_001__01 = 8,
      pex_bm_healthv2_inf_001__02 = 7,
      pex_bm_healthv2_inf_001__dk = c(NA, "999", "777")
    )
  )
  got <- e$birth_weight_covariate(d)
  expect_equal(got$bw_g, c(3401.94, NA, NA))
  expect_equal(got$participant_id[!is.na(got$bw_g)], "known")
})

test_that("birth-weight dropdown codes outside the supplied ranges are not measurements", {
  e <- local_sensitivity_helpers()
  d <- release_frame(
    "infant_health",
    data.frame(
      participant_id = c(
        "known",
        "bad-ounce",
        "bad-pound",
        "negative-ounce",
        "fractional-code"
      ),
      session_id = "ses-V02",
      pex_bm_healthv2_inf_001__01 = c(8, 16, 0, -1, 8.5),
      pex_bm_healthv2_inf_001__02 = c(7, 7, 21, 7, 7)
    )
  )
  expect_equal(e$birth_weight_covariate(d)$bw_g, c(3401.94, NA, NA, NA, NA))
})

test_that("prenatal exposure uses aggregate weeks 3-9 and retains positive evidence", {
  e <- local_sensitivity_helpers()
  ids <- c("pre", "early", "late", "partial", "zero", "unknown")
  d <- release_frame(
    "substance_use",
    data.frame(participant_id = ids, session_id = "ses-V01")
  )
  aggregate <- unlist(lapply(c("nic", "alc", "thc", "opd"), function(s) {
    paste0("pex_ch_tlfb_", s, "_wk_", sprintf("%02d", 3:9))
  }))
  d[c(1:3, 5), aggregate] <- 0
  # Actual competing subtype/classes and preconception fields must not leak in.
  d$pex_ch_tlfb_nic_cig_wk_03 <- 9
  d$pex_ch_tlfb_stim_wk_03 <- 9
  d$pex_ch_tlfb_nic_wk_02 <- 9
  d$pex_ch_tlfb_nic_wk_03[c(2, 4)] <- 1
  d$pex_ch_tlfb_opd_wk_09[3] <- 1
  got <- e$substance_exposure(d)
  got <- got[match(ids, got$participant_id), ]
  expect_equal(unname(got$tlfb_nic), c(0, 1, 0, 1, 0, NA))
  expect_equal(unname(got$tlfb_opd), c(0, 0, 1, NA, 0, NA))
  expect_equal(got$any_substance, c(0, 1, 1, 1, 0, NA))
  expect_equal(got$participant_id[got$any_substance %in% 0], c("pre", "zero"))
})

test_that("prenatal exposure is invariant to visit order and uses either visit", {
  e <- local_sensitivity_helpers()
  d <- release_frame(
    "substance_use",
    data.frame(
      participant_id = c("b", "a", "a"),
      session_id = c("ses-V01", "ses-V02", "ses-V01")
    )
  )
  d$pex_ch_tlfb_thc_wk_03 <- c(NA, 1, 0)
  first <- e$substance_exposure(d)
  second <- e$substance_exposure(d[3:1, ])
  first <- first[order(first$participant_id), ]
  second <- second[order(second$participant_id), ]
  rownames(first) <- rownames(second) <- NULL
  expect_equal(first, second)
  expect_equal(
    first$any_substance[match(c("a", "b"), first$participant_id)],
    c(1, NA)
  )
})

test_that("medication treatment follows paired names and codes, not dose fields", {
  e <- local_sensitivity_helpers()
  ids <- c("discovery", "visit2", "dose_only", "decimal")
  meds <- release_frame(
    "medications",
    data.frame(participant_id = ids, session_id = "ses-V01")
  )
  meds$pex_bm_health_preg__meds_003__01 <- c(
    "900001.0",
    "999999",
    "999999",
    "900001"
  )
  meds$pex_bm_health_preg__meds_003__01_name <- c(
    "labetalol",
    "unrelated",
    "unrelated",
    NA
  )
  meds$pex_bm_health_preg__meds_004__01 <- "800001"
  meds$pex_bm_health_preg__meds_004__01_name <- "unrelated"
  pregnancy <- release_frame(
    "pregnancy",
    data.frame(participant_id = ids[4:1], session_id = "ses-V02")
  )
  pregnancy$pex_bm_healthv2_preg__meds_047__01 <- c(
    "999999",
    "999999",
    "900001",
    "999999"
  )
  pregnancy$pex_bm_healthv2_preg__meds_047__02 <- "900001"
  got <- e$medication_exposure(meds, pregnancy)
  expect_equal(
    got$antihypertensive[match(ids, got$participant_id)],
    c(1, 1, 0, 1)
  )
})
