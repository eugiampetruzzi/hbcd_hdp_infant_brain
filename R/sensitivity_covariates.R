# Pure transformations shared by the sensitivity analysis and synthetic tests.
birth_weight_covariate <- function(infant) {
  # HBCD Health V2 Infancy: __01 ounces, __02 pounds.
  # https://docs.hbcdstudy.org/latest/instruments/pregexp/pex/
  ounces <- suppressWarnings(as.numeric(infant$pex_bm_healthv2_inf_001__01))
  pounds <- suppressWarnings(as.numeric(infant$pex_bm_healthv2_inf_001__02))
  unknown <- suppressWarnings(as.numeric(infant$pex_bm_healthv2_inf_001__dk))
  infant$bw_g <- (pounds * 16 + ounces) * 28.3495
  valid <- ounces %in% 0:15 & pounds %in% 0:20
  infant$bw_g[!valid | unknown %in% c(999, 777)] <- NA
  infant$bw_g[!(infant$bw_g > 400 & infant$bw_g < 6000)] <- NA
  infant
}

brain_volume_covariate <- function(volume) {
  volume$log_tbv <- log(volume$total_brain_volume_mm3)
  volume
}

substance_exposure <- function(tlfb) {
  substance <- data.frame(participant_id = unique(tlfb$participant_id))
  for (s in c("nic", "alc", "thc", "opd")) {
    weeks <- paste0("pex_ch_tlfb_", s, "_wk_", sprintf("%02d", 3:9))
    x <- sapply(tlfb[, weeks, drop = FALSE], function(v) {
      suppressWarnings(as.numeric(v))
    })
    used <- tapply(rowSums(x > 0, na.rm = TRUE) > 0, tlfb$participant_id, any)
    answered <- tapply(rowSums(!is.na(x)) > 0, tlfb$participant_id, any)
    substance[[paste0("tlfb_", s)]] <- as.numeric(ifelse(
      answered[substance$participant_id],
      as.numeric(used[substance$participant_id]),
      NA
    ))
  }
  x <- substance[, -1, drop = FALSE]
  # Positive evidence is conclusive; preserve unknown status for incomplete negatives.
  substance$any_substance <- ifelse(
    rowSums(x == 1, na.rm = TRUE) > 0,
    1,
    as.numeric(rowSums(x) > 0)
  )
  substance
}

medication_codes <- function(meds) {
  named <- grep("_name$", names(meds), value = TRUE)
  drug_names <- "labetalol|nifedipine|methyldopa|hydralazine|amlodipine|metoprolol|procardia|trandate|aldomet"
  found <- unlist(lapply(named, function(n) {
    meds[[sub("_name$", "", n)]][grepl(
      drug_names,
      meds[[n]],
      ignore.case = TRUE
    )]
  }))
  unique(sub(
    "\\.0$",
    "",
    as.character(c(found, "6185", "7417", "6876", "5470", "17767", "6918"))
  ))
}

medication_exposure <- function(meds, pregnancy) {
  codes <- medication_codes(meds)
  any_code <- function(x, pattern) {
    cols <- grep(pattern, names(x), value = TRUE)
    apply(x[, cols, drop = FALSE], 1, function(r) {
      any(sub("\\.0$", "", r) %in% codes)
    })
  }
  treated <- merge(
    data.frame(
      participant_id = meds$participant_id,
      v1 = any_code(meds, "^pex_bm_health_preg__meds_0[0-9]{2}__01$")
    ),
    data.frame(
      participant_id = pregnancy$participant_id,
      v2 = any_code(
        pregnancy,
        "^pex_bm_healthv2_preg__meds_0(4[7-9]|[5-9][0-9])__01$"
      )
    ),
    by = "participant_id",
    all = TRUE
  )
  # Ascertainment/skip-code policy remains unchanged pending the data dictionary.
  treated$antihypertensive <- as.numeric(
    treated$v1 %in% TRUE | treated$v2 %in% TRUE
  )
  treated
}
