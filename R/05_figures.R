# Figures 2, 3 and 4b. Figure 4a (parcel map) is drawn by python/figure4_map.py.

source("R/common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

CONTROL <- "#8f8d86"; HDP_RED <- "#C8102E"; INK <- "#000000"; BLUE <- "#2E6DB4"
theme_paper <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = INK),
        axis.line = element_line(colour = INK, linewidth = .4),
        axis.ticks = element_line(colour = INK, linewidth = .4),
        legend.title = element_blank(), strip.background = element_blank(),
        plot.tag = element_text(face = "bold", size = 11))
save_figure <- function(p, name, width, height) {
  for (ext in c("png", "pdf")) {
    ggsave(file.path(OUT, "figures", paste0(name, ".", ext)), p, width = width, height = height,
           units = "mm", dpi = 300, bg = "white", device = if (ext == "pdf") cairo_pdf)
  }
}
group_factor <- function(hdp) factor(ifelse(hdp == 1, "HDP", "Normotensive"), levels = c("Normotensive", "HDP"))

# Figure 2: age distributions in the final sample.
final <- read_table("analytic_sample.csv")
ages <- do.call(rbind, lapply(load_modalities(), function(x) x[, c(ID, HDP, GA, "candidate_age_weeks", PMA)]))
ages <- ages[!duplicated(ages[[ID]]) & ages[[ID]] %in% final[[ID]] & complete.cases(ages), ]
ages$group <- group_factor(ages[[HDP]])

age_histogram <- function(weeks, lo, hi, xlab, breaks, labels) {
  binned <- factor(pmin(pmax(floor(weeks), lo), hi), levels = lo:hi)
  counts <- as.data.frame(table(group = ages$group, week = binned))
  counts$week <- as.numeric(as.character(counts$week))
  counts$percent <- 100 * counts$Freq / ave(counts$Freq, counts$group, FUN = sum)
  ggplot(counts, aes(week, percent, fill = group)) +
    geom_col(position = position_dodge(width = .85), width = .8) +
    scale_fill_manual(values = c(CONTROL, HDP_RED)) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    scale_y_continuous(expand = expansion(mult = c(0, .06))) +
    labs(x = xlab, y = "Infants (%)") + theme_paper
}
fig2 <- (age_histogram(ages[[GA]], 32, 42, "Gestational age at birth (weeks)", seq(32, 42, 2),
                       c("≤32", 34, 36, 38, 40, "≥42")) |
         age_histogram(ages$candidate_age_weeks, 0, 14, "Postnatal age at scan (weeks)", seq(0, 14, 2),
                       c(0, 2, 4, 6, 8, 10, 12, "≥14")) |
         age_histogram(ages[[PMA]], 40, 50, "Postmenstrual age at scan (weeks)", seq(40, 50, 2),
                       c("≤40", 42, 44, 46, 48, "≥50"))) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "top")
save_figure(fig2, "figure2_delivery_timing", 183, 68)

# Figure 3: a-e GA curves at matched PMA (35-40 wk); f thickness difference split by GA.
curves <- read_table("delivery_timing_curves.csv")
curves$group <- group_factor(curves$hdp_o)
panels <- c("Mean diffusivity" = "Mean diffusivity", "Fractional anisotropy" = "Fractional anisotropy",
            "Mean kurtosis tensor" = "Mean kurtosis", "Cortical thickness" = "Cortical thickness",
            "Absolute sulcal depth" = "Sulcal depth")
curve_panel <- function(outcome) {
  ggplot(curves[curves$outcome == outcome, ], aes(ga, fit, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = .18, colour = NA) +
    geom_line(linewidth = .8) +
    scale_colour_manual(values = c(CONTROL, HDP_RED)) + scale_fill_manual(values = c(CONTROL, HDP_RED)) +
    scale_x_continuous(breaks = 35:40) + coord_cartesian(ylim = c(-.7, 1.05)) +
    labs(x = "Gestational age at birth (weeks)", y = paste(panels[[outcome]], "(SD)")) + theme_paper
}
global <- read_table("global_models.csv")
total <- global[global$model_set == "primary" & global$outcome == "Cortical thickness", ]
remaining <- global[global$model_set == "ga_adjusted" & global$outcome == "Cortical thickness", ]
split <- data.frame(part = factor(c("Remaining", "Accounted for"), levels = c("Accounted for", "Remaining")),
                    ymin = c(0, remaining$beta), ymax = c(remaining$beta, total$beta))
split_panel <- ggplot(split) +
  geom_rect(aes(xmin = .6, xmax = 1.4, ymin = ymin, ymax = ymax, fill = part)) +
  annotate("errorbar", x = 1, ymin = total$ci_lo, ymax = total$ci_hi, width = .15, linewidth = .4) +
  annotate("text", x = 1.47, y = (remaining$beta + total$beta) / 2, label = "Accounted for by\ngestational age",
           hjust = 0, size = 2.5, lineheight = .9) +
  annotate("text", x = 1.47, y = remaining$beta / 2, label = "Remaining", hjust = 0, size = 2.5) +
  geom_hline(yintercept = 0, linewidth = .4) +
  scale_fill_manual(values = c("#BDBDBD", "#008300"), guide = "none") +
  scale_x_continuous(limits = c(.5, 2.3), breaks = NULL) +
  scale_y_continuous(limits = c(0, .5), breaks = seq(0, .5, .1), expand = c(0, 0)) +
  labs(x = NULL, y = "Thicker cortex in HDP (SD)") + theme_paper + theme(axis.line.x = element_blank())
fig3 <- (Reduce(`+`, lapply(names(panels), curve_panel)) + split_panel) +
  plot_layout(ncol = 3, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "top")
save_figure(fig3, "figure3_delivery_timing_curves", 183, 120)

# Figure 4b: share of parcels per network with thicker cortex (networks with >= 5 parcels).
labels <- c(InfMotor = "Inferior motor", SupMotor = "Superior motor", PostPreMotor = "Posterior premotor",
            AntPreMotor = "Anterior premotor", DorsalAttn = "Dorsal attention", VentralAttn = "Ventral attention",
            AntDefault = "Anterior default", PostDefault = "Posterior default",
            LatOrbitoFrontal = "Lateral orbitofrontal", MedOrbitoFrontal = "Medial orbitofrontal",
            InfFrontoParietal = "Inferior frontoparietal", AntFrontoParietal = "Anterior frontoparietal",
            PostFrontoParietal = "Posterior frontoparietal", PreCentFrontoParietal = "Precentral frontoparietal",
            SupPrecuneus = "Superior precuneus", AntPrecuneus = "Anterior precuneus",
            InfPrecuneus = "Inferior precuneus", Auditory = "Auditory", Insula = "Insula", Salience = "Salience",
            Cingulate = "Cingulate", LateralVisual = "Lateral visual", MedialVisual = "Medial visual", None = "Unassigned")
networks <- read_table("parcel_networks.csv")
networks <- networks[networks$parcels >= 5, ]
networks$percent <- 100 * networks$significant / networks$parcels
networks$label <- factor(labels[networks$network], levels = labels[networks$network][order(networks$percent)])
fig4b <- ggplot(networks, aes(percent, label)) +
  geom_col(fill = BLUE, width = .7) +
  geom_text(aes(label = paste0(significant, "/", parcels)), hjust = -.15, size = 2.6) +
  scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Parcels with thicker cortex in HDP infants (%)", y = NULL) + theme_paper
save_figure(fig4b, "figure4b_networks", 95, 85)
