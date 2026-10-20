# Gestational age at birth at matched PMA: smooth GA associations, HDP x GA difference
# smooths (full range and 35-40 wk), 37-38 vs >=39 wk among term births, and the HDP
# thickness difference within GA bands, and a bootstrap (re-harmonized per resample) of the
# portion accounted for by GA.

source("R/common.R")

frames <- load_modalities()
final <- read_table("analytic_sample.csv")
in_final <- function(d, m) d[d[[ID]] %in% final[[ID]][final[[paste0("in_", tolower(sub(" MRI", "", m)))]]], ]

# HDP as an ordered factor so s(GA, by = hdp_o) is the HDP-minus-control difference curve.
timing_models <- function(range = NULL) {
  tests <- list(); curves <- list()
  k_ga <- if (is.null(range)) K else 5
  k_diff <- if (is.null(range)) 6 else 4
  for (m in names(MODALITIES)) {
    ys <- names(MODALITIES[[m]]$outcomes)
    d <- in_final(modality_sample(frames, m), m)
    if (!is.null(range)) d <- drop_small_batches(d[d[[GA]] >= range[1] & d[[GA]] <= range[2], ])
    d$hdp_o <- factor(d[[HDP]], levels = c(0, 1), ordered = TRUE)
    terms <- paste(smooth(PMA), smooth(GA, k_ga), sep = " + ")
    diff_term <- paste0("s(", GA, ", by = hdp_o, k = ", k_diff, ")")
    h <- harmonize(d, ys, paste("hdp_o + sex_m", terms, diff_term, sep = " + "))
    for (y in ys) {
      h$z <- zscore(h[[y]])
      rhs <- paste("hdp_o + sex_m", paste(MODALITIES[[m]]$qc, collapse = " + "), terms, sep = " + ")
      main <- gam(as.formula(paste("z ~", rhs)), data = h, method = "REML")
      inter <- gam(as.formula(paste("z ~", rhs, "+", diff_term)), data = h, method = "REML")
      s_main <- summary(main)$s.table; s_inter <- summary(inter)$s.table
      ga_row <- paste0("s(", GA, ")")
      tests[[length(tests) + 1]] <- data.frame(
        range = if (is.null(range)) "full" else paste(range, collapse = "-"),
        modality = m, outcome = MODALITIES[[m]]$outcomes[[y]], n = nrow(h), n_hdp = sum(h[[HDP]]),
        ga_edf = s_main[ga_row, "edf"], ga_p = s_main[ga_row, "p-value"],
        interaction_edf = s_inter[grep("hdp_o", rownames(s_inter)), "edf"],
        interaction_p = s_inter[grep("hdp_o", rownames(s_inter)), "p-value"])
      if (!is.null(range)) {
        grid <- expand.grid(ga = range[1]:range[2], hdp_o = levels(h$hdp_o), sex_m = c(0, 1))
        names(grid)[1] <- GA
        grid$hdp_o <- factor(grid$hdp_o, levels = levels(h$hdp_o), ordered = TRUE)
        grid[[PMA]] <- median(h[[PMA]])
        for (q in MODALITIES[[m]]$qc) grid[[q]] <- mean(h[[q]])
        # Average the prediction matrix over sex, so the SE accounts for covariance.
        X <- predict(inter, grid, type = "lpmatrix")
        cell <- interaction(grid$hdp_o, grid[[GA]], drop = TRUE, lex.order = TRUE)
        L <- rowsum(X, cell) / as.vector(table(cell))
        avg <- unique(data.frame(hdp_o = grid$hdp_o, ga = grid[[GA]], cell = cell))
        avg <- avg[match(rownames(L), avg$cell), c("hdp_o", "ga")]
        avg$fit <- drop(L %*% coef(inter))
        avg$se <- sqrt(rowSums((L %*% vcov(inter)) * L))
        curves[[length(curves) + 1]] <- cbind(outcome = MODALITIES[[m]]$outcomes[[y]], avg)
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
  h <- harmonize(d, ys, paste(HDP, "sex_m", "early_term", smooth(PMA), sep = " + "))
  for (y in ys) {
    term[[length(term) + 1]] <- cbind(
      outcome = MODALITIES[[m]]$outcomes[[y]], n_early_term = sum(h$early_term),
      fit_term(h, y, primary_rhs(m, extra = "early_term"), term = "early_term"))
  }
}
term <- do.call(rbind, term)
term$q <- p.adjust(term$p, "BH")
write_table(term, "term_contrast.csv")

# HDP thickness difference within GA bands, with the primary batch and smoothing rules.
m <- "Structural MRI"
bands <- list("<37" = c(0, 36), "37-38" = c(37, 38), ">=39" = c(39, 99))
band <- do.call(rbind, lapply(names(bands), function(b) {
  d <- in_final(modality_sample(frames, m), m)
  d <- drop_small_batches(d[d[[GA]] >= bands[[b]][1] & d[[GA]] <= bands[[b]][2], ])
  if (sum(d[[HDP]]) < 20) return(NULL)
  h <- harmonize(d, names(MODALITIES[[m]]$outcomes), harmonization_rhs())
  cbind(band = b, fit_term(h, "global_area_weighted_thickness", primary_rhs(m)))
}))
write_table(band, "thickness_by_ga_band.csv")

# Portion of the HDP thickness difference accounted for by GA. Both coefficients are on one
# scale (mm, re-expressed in SDs of the primary harmonized outcome). The bootstrap resamples
# infants within scanner batch and re-harmonizes in every resample, so intervals include
# harmonization uncertainty.
B <- 10000
d <- in_final(modality_sample(frames, m), m)
ys <- names(MODALITIES[[m]]$outcomes)
y <- "global_area_weighted_thickness"
sd_pop <- function(v) sqrt(mean((v - mean(v))^2))
hdp_mm <- function(h, ga) {
  coef(gam(as.formula(paste(y, "~", primary_rhs(m, ga = ga))), data = h, method = "REML"))[[HDP]]
}
decompose <- function(d) {
  h_total <- harmonize(d, ys, harmonization_rhs()); h_ga <- harmonize(d, ys, harmonization_rhs(ga = TRUE))
  s <- sd_pop(h_total[[y]])
  c(total = hdp_mm(h_total, FALSE) / s, remaining = hdp_mm(h_ga, TRUE) / s)
}
set.seed(20261018)
est <- decompose(d)
by_batch <- split(seq_len(nrow(d)), d$batch)
boot <- replicate(B, decompose(d[unlist(lapply(by_batch, function(i) i[sample.int(length(i), replace = TRUE)])), ]))
accounted <- boot["total", ] - boot["remaining", ]
ci <- function(v) quantile(v, c(.025, .975), names = FALSE)
portion <- data.frame(resamples = B, total = est[["total"]], total_lo = ci(boot["total", ])[1],
                      total_hi = ci(boot["total", ])[2], remaining = est[["remaining"]],
                      remaining_lo = ci(boot["remaining", ])[1], remaining_hi = ci(boot["remaining", ])[2],
                      accounted = est[["total"]] - est[["remaining"]],
                      accounted_lo = ci(accounted)[1], accounted_hi = ci(accounted)[2],
                      proportion = 1 - est[["remaining"]] / est[["total"]],
                      proportion_lo = ci(accounted / boot["total", ])[1],
                      proportion_hi = ci(accounted / boot["total", ])[2])
write_table(portion, "ga_portion_bootstrap.csv")

print(rbind(full$tests, supported$tests)[, c("range", "outcome", "ga_q", "interaction_q")], digits = 2)
print(term[, c("outcome", "beta", "ci_lo", "ci_hi", "q")], digits = 2)
print(band[, c("band", "n", "n_hdp", "beta", "ci_lo", "ci_hi", "p")], digits = 2)
print(portion, digits = 2)
