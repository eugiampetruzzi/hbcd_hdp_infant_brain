# Gestational age at birth at matched PMA: smooth GA associations, HDP x GA difference
# smooths (full range and 35-40 wk), 37-38 vs >=39 wk among term births, and the HDP
# thickness difference within GA bands, and a bootstrap of the portion accounted for by GA.

source("R/common.R")

frames <- load_modalities()
final <- read_table("analytic_sample.csv")
in_final <- function(d, m) {
  d[
    d[[ID]] %in%
      final[[ID]][final[[paste0("in_", tolower(sub(" MRI", "", m)))]]],
  ]
}

# HDP as an ordered factor so s(GA, by = hdp_o) is the HDP-minus-control difference curve.
timing_models <- function(range = NULL) {
  tests <- list()
  curves <- list()
  k_ga <- if (is.null(range)) K else 5
  k_diff <- if (is.null(range)) 6 else 4
  for (m in names(MODALITIES)) {
    ys <- names(MODALITIES[[m]]$outcomes)
    d <- in_final(modality_sample(frames, m), m)
    if (!is.null(range)) {
      d <- drop_small_batches(d[d[[GA]] >= range[1] & d[[GA]] <= range[2], ])
    }
    d$hdp_o <- factor(d[[HDP]], levels = c(0, 1), ordered = TRUE)
    terms <- paste(smooth(PMA), smooth(GA, k_ga), sep = " + ")
    diff_term <- paste0("s(", GA, ", by = hdp_o, k = ", k_diff, ")")
    h <- harmonize(d, ys, paste("hdp_o + sex_m", terms, diff_term, sep = " + "))
    for (y in ys) {
      h$z <- zscore(h[[y]])
      rhs <- paste(
        "hdp_o + sex_m",
        paste(MODALITIES[[m]]$qc, collapse = " + "),
        terms,
        sep = " + "
      )
      main <- gam(as.formula(paste("z ~", rhs)), data = h, method = "REML")
      inter <- gam(
        as.formula(paste("z ~", rhs, "+", diff_term)),
        data = h,
        method = "REML"
      )
      s_main <- summary(main)$s.table
      s_inter <- summary(inter)$s.table
      ga_row <- paste0("s(", GA, ")")
      tests[[length(tests) + 1]] <- data.frame(
        range = if (is.null(range)) "full" else paste(range, collapse = "-"),
        modality = m,
        outcome = MODALITIES[[m]]$outcomes[[y]],
        n = nrow(h),
        n_hdp = sum(h[[HDP]]),
        ga_edf = s_main[ga_row, "edf"],
        ga_p = s_main[ga_row, "p-value"],
        interaction_edf = s_inter[grep("hdp_o", rownames(s_inter)), "edf"],
        interaction_p = s_inter[grep("hdp_o", rownames(s_inter)), "p-value"]
      )
      if (!is.null(range)) {
        grid <- expand.grid(
          ga = range[1]:range[2],
          hdp_o = levels(h$hdp_o),
          sex_m = c(0, 1)
        )
        names(grid)[1] <- GA
        grid$hdp_o <- factor(
          grid$hdp_o,
          levels = levels(h$hdp_o),
          ordered = TRUE
        )
        grid[[PMA]] <- median(h[[PMA]])
        for (q in MODALITIES[[m]]$qc) {
          grid[[q]] <- mean(h[[q]])
        }
        pr <- predict(inter, grid, se.fit = TRUE)
        grid$fit <- pr$fit
        grid$se <- pr$se.fit
        avg <- aggregate(
          cbind(fit, se) ~ hdp_o + grid[[GA]],
          data = grid,
          FUN = mean
        )
        names(avg)[2] <- "ga"
        curves[[length(curves) + 1]] <- cbind(
          outcome = MODALITIES[[m]]$outcomes[[y]],
          avg
        )
      }
    }
  }
  tests <- do.call(rbind, tests)
  tests$ga_q <- p.adjust(tests$ga_p, "BH")
  tests$interaction_q <- p.adjust(tests$interaction_p, "BH")
  list(tests = tests, curves = if (length(curves)) do.call(rbind, curves))
}

full <- timing_models()
supported <- timing_models(c(35, 40))
write_table(rbind(full$tests, supported$tests), "delivery_timing_tests.csv")
write_table(supported$curves, "delivery_timing_curves.csv")

# Term births only: GA is linear in harmonization because few distinct weeks remain.
term <- list()
for (m in c("Diffusion MRI", "Structural MRI")) {
  ys <- names(MODALITIES[[m]]$outcomes)
  d <- in_final(modality_sample(frames, m), m)
  d <- drop_small_batches(d[d[[GA]] >= 37, ])
  d$early_term <- as.numeric(d[[GA]] < 39)
  h <- harmonize(
    d,
    ys,
    paste(HDP, "sex_m", "early_term", smooth(PMA), sep = " + ")
  )
  for (y in ys) {
    term[[length(term) + 1]] <- cbind(
      outcome = MODALITIES[[m]]$outcomes[[y]],
      n_early_term = sum(h$early_term),
      fit_term(h, y, primary_rhs(m, extra = "early_term"), term = "early_term")
    )
  }
}
term <- do.call(rbind, term)
term$q <- p.adjust(term$p, "BH")
write_table(term, "term_contrast.csv")

# HDP thickness difference within GA bands, with the primary batch and smoothing rules.
m <- "Structural MRI"
bands <- list("<37" = c(0, 36), "37-38" = c(37, 38), ">=39" = c(39, 99))
band <- do.call(
  rbind,
  lapply(names(bands), function(b) {
    d <- in_final(modality_sample(frames, m), m)
    d <- drop_small_batches(d[
      d[[GA]] >= bands[[b]][1] & d[[GA]] <= bands[[b]][2],
    ])
    if (sum(d[[HDP]]) < 20) {
      return(NULL)
    }
    h <- harmonize(d, names(MODALITIES[[m]]$outcomes), harmonization_rhs())
    cbind(
      band = b,
      fit_term(h, "global_area_weighted_thickness", primary_rhs(m))
    )
  })
)
write_table(band, "thickness_by_ga_band.csv")

# Portion of the HDP thickness difference accounted for by GA: 10,000 bootstrap resamples, refitting
# the primary and GA-adjusted models on data harmonized once in the full sample.
d <- in_final(modality_sample(frames, m), m)
ys <- names(MODALITIES[[m]]$outcomes)
y <- "global_area_weighted_thickness"
h_total <- harmonize(d, ys, harmonization_rhs())
h_ga <- harmonize(d, ys, harmonization_rhs(ga = TRUE))
hdp_beta <- function(h, ga) {
  h$z <- zscore(h[[y]])
  coef(gam(
    as.formula(paste("z ~", primary_rhs(m, ga = ga))),
    data = h,
    method = "REML"
  ))[[HDP]]
}
set.seed(20261018)
total <- hdp_beta(h_total, FALSE)
remaining <- hdp_beta(h_ga, TRUE)
boot <- replicate(10000, {
  i <- sample(nrow(h_total), replace = TRUE)
  c(hdp_beta(h_total[i, ], FALSE), hdp_beta(h_ga[i, ], TRUE))
})
accounted <- boot[1, ] - boot[2, ]
portion <- data.frame(
  total = total,
  remaining = remaining,
  accounted = total - remaining,
  accounted_lo = quantile(accounted, .025),
  accounted_hi = quantile(accounted, .975),
  proportion = (total - remaining) / total,
  proportion_lo = quantile(accounted / boot[1, ], .025),
  proportion_hi = quantile(accounted / boot[1, ], .975),
  row.names = NULL
)
write_table(portion, "ga_portion_bootstrap.csv")

print(
  rbind(full$tests, supported$tests)[, c(
    "range",
    "outcome",
    "ga_q",
    "interaction_q"
  )],
  digits = 2
)
print(term[, c("outcome", "beta", "ci_lo", "ci_hi", "q")], digits = 2)
print(
  band[, c("band", "n", "n_hdp", "beta", "ci_lo", "ci_hi", "p")],
  digits = 2
)
print(portion, digits = 2)
