test_that("prenatal exposure uses aggregate weeks 3-9 and retains positive evidence", {
  e <- local_sensitivity_helpers()
  ids <- c("pre", "early", "late", "partial", "zero", "unknown")
  d <- release_frame("substance_use", data.frame(participant_id = ids, session_id = "ses-V01"))
  aggregate <- unlist(lapply(c("nic", "alc", "thc", "opd"), function(s)
    paste0("pex_ch_tlfb_", s, "_wk_", sprintf("%02d", 3:9))))
  d[c(1:3, 5), aggregate] <- 0
  # Actual competing subtype/classes and preconception fields must not leak in.
  d$pex_ch_tlfb_nic_cig_wk_03 <- 9
  d$pex_ch_tlfb_stim_wk_03 <- 9
  d$pex_ch_tlfb_nic_wk_02 <- 9
  d$pex_ch_tlfb_nic_wk_03[c(2, 4)] <- 1
  d$pex_ch_tlfb_opd_wk_09[3] <- 1
  got <- e$substance_exposure(d)
  got <- got[match(ids, got$participant_id), ]
  expect_equal(got$tlfb_nic, c(0, 1, 0, 1, 0, NA))
  expect_equal(got$tlfb_opd, c(0, 0, 1, NA, 0, NA))
  expect_equal(got$any_substance, c(0, 1, 1, 1, 0, NA))
  expect_equal(got$participant_id[got$any_substance %in% 0], c("pre", "zero"))
})

test_that("prenatal exposure is invariant to visit order and uses either visit", {
  e <- local_sensitivity_helpers()
  d <- release_frame("substance_use", data.frame(participant_id = c("b", "a", "a"),
                                                session_id = c("ses-V01", "ses-V02", "ses-V01")))
  d$pex_ch_tlfb_thc_wk_03 <- c(NA, 1, 0)
  first <- e$substance_exposure(d)
  second <- e$substance_exposure(d[3:1, ])
  expect_equal(first[order(first$participant_id), ], second[order(second$participant_id), ],
               ignore_attr = TRUE)
  expect_equal(first$any_substance[match(c("a", "b"), first$participant_id)], c(1, NA))
})

test_that("medication treatment follows paired names and codes, not dose fields", {
  e <- local_sensitivity_helpers()
  ids <- c("discovery", "visit2", "dose_only", "decimal")
  meds <- release_frame("medications", data.frame(participant_id = ids, session_id = "ses-V01"))
  meds$pex_bm_health_preg__meds_003__01 <- c("900001.0", "999999", "999999", "900001")
  meds$pex_bm_health_preg__meds_003__01_name <- c("labetalol", "unrelated", "unrelated", NA)
  meds$pex_bm_health_preg__meds_004__01 <- "800001"
  meds$pex_bm_health_preg__meds_004__01_name <- "unrelated"
  pregnancy <- release_frame("pregnancy", data.frame(participant_id = ids[4:1], session_id = "ses-V02"))
  pregnancy$pex_bm_healthv2_preg__meds_047__01 <- c("999999", "999999", "900001", "999999")
  pregnancy$pex_bm_healthv2_preg__meds_047__02 <- "900001"
  got <- e$medication_exposure(meds, pregnancy)
  expect_equal(got$antihypertensive[match(ids, got$participant_id)], c(1, 1, 0, 1))
})
