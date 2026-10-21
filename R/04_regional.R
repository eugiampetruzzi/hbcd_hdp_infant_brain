# Bundle- and parcel-level HDP differences, and whether the parcel thickness map tracks
# regional effects of image quality or scanner.

source("R/common.R")
source("R/regional_summary.R")

COVERAGE <- 0.80
REGIONAL <- list(
  list(
    metric = "Bundle FA",
    modality = "Diffusion MRI",
    file = "outcomes/diffusion/v02_bundle_fa_masked_median.tsv",
    prefix = "fa_"
  ),
  list(
    metric = "Bundle MD",
    modality = "Diffusion MRI",
    file = "outcomes/diffusion/v02_bundle_md_masked_median.tsv",
    prefix = "md_"
  ),
  list(
    metric = "Bundle MKT",
    modality = "Diffusion MRI",
    file = "outcomes/diffusion/v02_bundle_mkt_masked_median.tsv",
    prefix = "mkt_"
  ),
  list(
    metric = "Parcel cortical thickness",
    modality = "Structural MRI",
    file = "outcomes/structural/v02_myerslabonte_thickness.tsv",
    prefix = ""
  ),
  list(
    metric = "Parcel sulcal depth",
    modality = "Structural MRI",
    file = "outcomes/structural/v02_myerslabonte_sulc.tsv",
    prefix = ""
  )
)

frames <- load_modalities()
final <- read_table("analytic_sample.csv")

regional_data <- function(r) {
  m <- r$modality
  base <- modality_sample(frames, m)
  base <- base[
    base[[ID]] %in%
      final[[ID]][final[[paste0("in_", tolower(sub(" MRI", "", m)))]]],
  ]
  x <- read_input(r$file)
  x <- x[!duplicated(x[[ID]]), ]
  feats <- setdiff(
    names(x)[sapply(x, is.numeric)],
    c(ID, "session_id", "source_run")
  )
  if (nzchar(r$prefix)) {
    feats <- feats[startsWith(feats, r$prefix)]
  }
  d <- merge(
    base[, c(ID, HDP, PMA, "sex_m", MODALITIES[[m]]$qc, "batch")],
    x[, c(ID, feats)],
    by = ID
  )
  feats <- feats[colMeans(!is.na(d[, feats])) >= COVERAGE]
  list(d = d[rowSums(!is.na(d[, feats])) > 0, ], feats = feats)
}

# Missing values are filled with the batch median for harmonization only, then restored.
harmonize_regional <- function(d, feats) {
  missing <- is.na(d[, feats])
  for (f in feats) {
    med <- ave(d[[f]], d$batch, FUN = function(v) median(v, na.rm = TRUE))
    d[[f]][missing[, f]] <- med[missing[, f]]
  }
  raw <- d
  d <- harmonize(d, feats, harmonization_rhs(), eb = TRUE)
  d[, feats][missing] <- NA
  raw[, feats][missing] <- NA
  list(harmonized = d, raw = raw)
}

results <- list()
for (r in REGIONAL) {
  rd <- regional_data(r)
  hz <- harmonize_regional(rd$d, rd$feats)
  m <- r$modality
  for (f in rd$feats) {
    h <- hz$harmonized[!is.na(hz$harmonized[[f]]), ]
    est <- fit_term(h, f, primary_rhs(m))
    row <- data.frame(
      metric = r$metric,
      feature = sub(paste0("^", r$prefix), "", f),
      est,
      qc_beta = NA,
      scanner_eta2 = NA
    )
    if (r$metric == "Parcel cortical thickness") {
      row$qc_beta <- fit_term(
        h,
        f,
        primary_rhs(m),
        term = MODALITIES[[m]]$qc
      )$beta
      raw <- hz$raw[!is.na(hz$raw[[f]]), ]
      raw$y <- raw[[f]]
      e <- residuals(gam(
        as.formula(paste("y ~", harmonization_rhs())),
        data = raw,
        method = "REML"
      ))
      row$scanner_eta2 <- summary(lm(e ~ factor(raw$batch)))$r.squared
    }
    results[[length(results) + 1]] <- row
  }
}
results <- do.call(rbind, results)
results$q <- ave(results$p, results$metric, FUN = function(p) p.adjust(p, "BH"))
parcels <- results$metric == "Parcel cortical thickness"
# Myers-Labonte parcel names are <hemisphere>_<network>_<index>; parcel 1 has no prefix
# and lies in the left hemisphere.
results$network <- ifelse(
  parcels,
  sub("^([LR]_)?(.*)_[0-9]+$", "\\2", results$feature),
  NA
)
results$hemisphere <- ifelse(
  parcels,
  ifelse(grepl("^R_", results$feature), "R", "L"),
  NA
)
write_table(results, "regional.csv")

thick <- results[parcels, ]
# Descriptive only: parcels are spatially dependent and share infants, so no p-values.
artifact <- data.frame(
  comparison = c("image quality", "scanner"),
  rho = c(
    cor(thick$beta, thick$qc_beta, method = "spearman"),
    cor(thick$beta, thick$scanner_eta2, method = "spearman")
  )
)
write_table(artifact, "parcel_artifact_check.csv")

networks <- parcel_network_summary(results)
write_table(networks, "parcel_networks.csv")

print(aggregate(
  cbind(significant = q < .05, total = 1) ~ metric,
  data = results,
  FUN = sum
))
print(artifact, digits = 2)
