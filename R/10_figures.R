# ---- Figures ----------------------------------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages(library(ggplot2))
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

cal  <- readRDS("output/calibration.rds")
diag <- utils::read.csv(file.path(DIR_TABLES, "estimator_diagnostics.csv"))

label_of <- function(d) {
  ifelse(is.na(d$ps_spec) & is.na(d$om_spec), d$estimator,
  ifelse(is.na(d$om_spec), paste0(d$estimator, " [", d$ps_spec, "]"),
  ifelse(is.na(d$ps_spec), paste0(d$estimator, " [", d$om_spec, "]"),
         paste0(d$estimator, " [", d$ps_spec, "/", d$om_spec, "]"))))
}

th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, colour = "grey35"))

# ---- 1. Relative bias against the guardrail band ----------------------------
# Relative bias, not share-of-gap-closed. The naive readout is 350-590% off, so
# closing 90% of that gap still leaves a 50% overstatement - a framing that
# reads as success while the answer remains unusable.
d1 <- cal[cal$outcome == "conversion" & cal$alpha > 0 & cal$estimand == "ATE" &
          (is.na(cal$om_spec) | cal$om_spec == "rich"), ]
d1$label <- label_of(d1)
d1$arm <- factor(d1$arm, levels = c("observed", "hidden"),
                 labels = c("observed confounding\n(recoverable in principle)",
                            "hidden confounder\n(f0 withheld from estimators)"))
band <- 100 * GUARDRAILS$max_relative_bias
p1 <- ggplot(d1, aes(factor(alpha), 100 * relative_bias, group = label, colour = label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -band, ymax = band,
           fill = "seagreen", alpha = 0.15) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~arm) +
  scale_y_continuous(breaks = scales::pretty_breaks(8)) +
  labs(title = "Relative bias against the randomized benchmark",
       subtitle = paste0("Conversion outcome. Green band is the pre-specified ",
                         band, "% tolerance. 0 = recovered the experiment."),
       x = "targeting strength (alpha)", y = "relative bias (%)", colour = NULL,
       caption = paste(
         "The naive line is the unadjusted readout, 350-590% too high.",
         "The oracle propensity score is known only because this study wrote the selection rule;",
         "it is an upper bound, not something a real analyst can compute.", sep = "\n")) + th +
  theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
message("  figure 1/3"); ggsave(file.path(DIR_FIGURES, "relative_bias.png"), p1,
                                width = 11, height = 7.5, dpi = 150)

# ---- 2. Estimates against the benchmark, with intervals ----------------------
d2 <- cal[cal$outcome == "conversion" & cal$arm == "observed" &
          cal$estimand == "ATE" & (is.na(cal$om_spec) | cal$om_spec == "rich"), ]
d2$label <- label_of(d2)
p2 <- ggplot(d2, aes(x = 100 * estimate, y = reorder(label, estimate))) +
  geom_vline(aes(xintercept = 100 * benchmark), colour = "firebrick", linewidth = 0.7) +
  geom_errorbar(aes(xmin = 100 * lower, xmax = 100 * upper),
                orientation = "y", colour = "grey40", linewidth = 0.4, width = 0) +
  geom_point(size = 2) +
  facet_wrap(~paste0("alpha = ", alpha), scales = "free_x", nrow = 1) +
  labs(title = "Estimates against the randomized benchmark, observed-confounding arm",
       subtitle = "Red line is the benchmark for that subsample. Bars are 95% intervals.",
       x = "effect on conversion (percentage points)", y = NULL) + th
message("  figure 2/3"); ggsave(file.path(DIR_FIGURES, "calibration.png"), p2, width = 12, height = 5.5, dpi = 150)

# ---- 3. Balance before and after weighting ----------------------------------
d3 <- diag[diag$alpha > 0, ]
d3$arm <- factor(d3$arm, levels = c("observed", "hidden"))
p3 <- ggplot(d3, aes(x = factor(alpha))) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "firebrick") +
  geom_point(aes(y = max_abs_smd_unadjusted, shape = "before weighting"),
             size = 2.5, colour = "grey30") +
  geom_point(aes(y = max_abs_smd_weighted, colour = ps_spec, shape = "after weighting"),
             size = 2.5) +
  facet_wrap(~arm) +
  scale_shape_manual(values = c("before weighting" = 4, "after weighting" = 16)) +
  labs(title = "Worst covariate imbalance before and after weighting",
       subtitle = "Dashed line is the pre-specified 0.10 threshold.",
       x = "targeting strength (alpha)", y = "max |standardized mean difference|",
       colour = "propensity spec", shape = NULL) + th
message("  figure 3/3"); ggsave(file.path(DIR_FIGURES, "balance.png"), p3, width = 10, height = 5, dpi = 150)

message("Wrote 3 figures to ", DIR_FIGURES)
