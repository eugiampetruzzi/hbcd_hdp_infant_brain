test_that("real scanner schemas preserve software batches and cannot inflate cohort thresholds", {
  a <- local_analysis()
  frames <- synthetic_frames(22)
  write_modality_inputs(a, frames)
  withr::local_envvar(c(PYTHONPATH = file.path(a$repo, "python")))
  ids <- frames[[1]]$participant_id
  for (known in c(9L, 10L)) {
    for (i in seq_along(ids)) {
      variant <- if (i <= known) 1L else if (i <= known + 10L) 2L else 3L
      scans <- release_frame("scanner_metadata", data.frame(
        filename = c("dwi/a_dwi.nii.gz", "anat/a_T2w.nii.gz", "func/a_bold.nii.gz")), variant)
      scans$site <- "same-site"
      if (variant < 3L) {
        scans$scanner_serial_number <- "S1"
        scans$scanner_software <- if (variant == 1L) "v1" else "v2"
        scans$scanner_manufacturer <- "Siemens"
        scans$scanner_model <- "synthetic"
      }
      write_release_fixture(a, paste0("imaging/qc_mriqc/hbcd/rawdata/", ids[i],
                                      "/ses-V02/", ids[i], "_ses-V02_scans.tsv"), scans)
    }
    output <- system2("uv", c("run", "--project", shQuote(a$repo), "python",
                              shQuote(file.path(a$repo, "prepare/02_scanners.py"))),
                      stdout = TRUE, stderr = TRUE)
    expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
    loaded <- a$env$load_modalities()
    for (modality in names(loaded)) {
      d <- loaded[[modality]]
      expect_true(all(is.na(d$batch[match(ids[(known + 11):22], d$participant_id)])))
      sample <- a$env$modality_sample(loaded, modality)
      expected <- if (known == 9L) ids[10:19] else ids[1:20]
      expect_equal(sample$participant_id, expected)
      expect_equal(unique(sample$batch), if (known == 9L) "S1 | v2" else c("S1 | v1", "S1 | v2"))
    }
  }
})
