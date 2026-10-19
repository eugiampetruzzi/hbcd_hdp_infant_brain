# HDP differences in the seven global measures at matched PMA (primary), and after adding
# a smooth of gestational age at birth. Each model set uses its own harmonized data.

source("R/common.R")

partial_r2 <- function(d, y, m, ga) {
  d$z <- zscore(d[[y]])
  full <- gam(as.formula(paste("z ~", primary_rhs(m, ga = ga))), data = d, method = "REML")
  reduced <- gam(as.formula(paste("z ~", sub(paste0(HDP, " \\+ "), "", primary_rhs(m, ga = ga)))),
                 data = d, method = "REML")
  1 - sum(residuals(full)^2) / sum(residuals(reduced)^2)
}

# Share of covariate-adjusted variance explained by batch, before or after harmonization.
batch_variance <- function(d, ys, ga) {
  do.call(rbind, lapply(ys, function(y) {
    d$y <- d[[y]]
    e <- residuals(gam(as.formula(paste("y ~", harmonization_rhs(ga = ga))), data = d, method = "REML"))
    eta2 <- summary(lm(e ~ factor(d$batch)))$r.squared
    data.frame(outcome = y, eta2 = eta2, cohens_f = sqrt(eta2 / (1 - eta2)))
  }))
}

frames <- load_modalities()
results <- list(); diagnostics <- list(); batches <- list(); harmonized <- list()

for (m in names(MODALITIES)) {
  ys <- names(MODALITIES[[m]]$outcomes)
  raw <- frames[[m]]
  complete <- raw[complete.cases(raw[, c(ID, HDP, GA, PMA, "sex_m", MODALITIES[[m]]$qc, ys)]), ]
  d <- modality_sample(frames, m)
  batches[[m]] <- data.frame(modality = m, n_complete = nrow(complete),
                             multi_scanner = sum(is.na(complete$batch)),
                             small_batch = sum(!is.na(complete$batch)) - nrow(d),
                             n = nrow(d), n_batches = length(unique(d$batch)))

  for (ga in c(FALSE, TRUE)) {
    h <- harmonize(d, ys, harmonization_rhs(ga = ga))
    if (!ga) harmonized[[m]] <- h[, c(ID, ys)]
    set <- if (ga) "ga_adjusted" else "primary"
    diagnostics[[length(diagnostics) + 1]] <- rbind(
      cbind(model_set = set, modality = m, stage = "before", batch_variance(d, ys, ga)),
      cbind(model_set = set, modality = m, stage = "after", batch_variance(h, ys, ga)))
    for (y in ys) {
      results[[length(results) + 1]] <- cbind(
        model_set = set, modality = m, outcome = MODALITIES[[m]]$outcomes[[y]], outcome_var = y,
        fit_term(h, y, primary_rhs(m, ga = ga)), partial_r2 = partial_r2(h, y, m, ga))
    }
  }
}

results <- do.call(rbind, results)
results$n_control <- results$n - results$n_hdp
results$q <- ave(results$p, results$model_set, FUN = function(p) p.adjust(p, "BH"))
# Minimum detectable effect: 80% power at two-sided alpha = .05.
results$mde <- 2.8 * results$se
write_table(results, "global_models.csv")
write_table(do.call(rbind, diagnostics), "harmonization_diagnostics.csv")
write_table(do.call(rbind, batches), "analytic_samples.csv")

# Participant-level file (protected data): sample membership and harmonized measures.
sample <- Reduce(function(a, b) merge(a, b, by = ID, all = TRUE), harmonized)
for (m in names(harmonized)) {
  sample[[paste0("in_", tolower(sub(" MRI", "", m)))]] <- sample[[ID]] %in% harmonized[[m]][[ID]]
}
write_table(sample, "analytic_sample.csv")

print(results[, c("model_set", "outcome", "n", "n_hdp", "beta", "ci_lo", "ci_hi", "p", "q")], digits = 3)
