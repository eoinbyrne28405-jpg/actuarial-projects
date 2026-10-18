# =============================================================================
# 02_frequency_model.R
# Poisson GLM for claim frequency with log(exposure) offset.
#
#   ClaimNb ~ Poisson(lambda * Exposure),  log(lambda) = x'beta
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
d     <- readRDS("data/model_data.rds")
train <- d[split == "train"]
test  <- d[split == "test"]

poisson_deviance <- function(y, mu) {
  2 * mean(ifelse(y == 0, 0, y * log(y / mu)) - (y - mu))
}

# --- 1. One-way analysis: observed frequency by factor -------------------------
oneway <- function(var) {
  d[, .(freq = sum(ClaimNb) / sum(Exposure), exposure = sum(Exposure)), by = var][
    , .(level = as.character(get(var)), freq, exposure, variable = var)]
}
ow <- rbindlist(lapply(c("DrivAgeBand", "VehAgeBand", "VehPowerF", "Area", "VehGas"), oneway))
# Order levels naturally: bands by lower bound, numbers numerically, else A-Z
sort_key <- function(x) {
  num <- suppressWarnings(as.numeric(sub("^[(\\[]?(-?[0-9.]+).*$", "\\1", x)))
  ifelse(is.na(num), rank(x), num)
}
ow[, key := sort_key(level), by = variable]
ow[, plot_level := paste(variable, level, sep = "::")]
ow[, plot_level := factor(plot_level, levels = plot_level[order(variable, key)])]
p_ow <- ggplot(ow, aes(plot_level, freq)) +
  scale_x_discrete(labels = function(x) sub("^.*::", "", x)) +
  geom_col(fill = "#6baed6", width = 0.7) +
  facet_wrap(~ variable, scales = "free_x", ncol = 2) +
  labs(title = "Observed claim frequency by rating factor (one-way)",
       x = NULL, y = "Claims per policy-year") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/figures/oneway_frequency.png", p_ow, width = 9, height = 8, dpi = 150)

# Bonus-malus deserves its own plot: it's the strongest single predictor
bm <- d[, .(freq = sum(ClaimNb) / sum(Exposure), exposure = sum(Exposure)),
        by = .(bm = pmin(10 * floor(BonusMalus / 10), 150))]
p_bm <- ggplot(bm[exposure > 200], aes(bm, freq)) +
  geom_point(aes(size = exposure), colour = "#08519c") +
  geom_smooth(method = "loess", se = FALSE, colour = "grey40", formula = y ~ x) +
  scale_size_area(max_size = 6, guide = "none") +
  labs(title = "Claim frequency rises steeply with the bonus-malus level",
       x = "Bonus-malus level (50 = best, 100 = neutral)", y = "Claims per policy-year") +
  theme_minimal(base_size = 11)
ggsave("output/figures/bonus_malus_frequency.png", p_bm, width = 8, height = 5, dpi = 150)

# --- 2. Models -----------------------------------------------------------------
f_null <- ClaimNb ~ 1
f_full <- ClaimNb ~ DrivAgeBand + VehAgeBand + VehPowerF + VehBrand + VehGas +
                    BonusMalusC + LogDensity + Region + Area

m_null <- glm(f_null, family = poisson(), data = train, offset = log(Exposure))
m_full <- glm(f_full, family = poisson(), data = train, offset = log(Exposure))

# Area and LogDensity both measure urbanisation; test whether Area adds anything
m_noarea <- update(m_full, . ~ . - Area)
lr <- anova(m_noarea, m_full, test = "Chisq")
print(lr)
area_p <- lr$`Pr(>Chi)`[2]
m_final <- if (area_p > 0.05) m_noarea else m_full
cat("Area p-value:", signif(area_p, 3), "->", if (area_p > 0.05) "dropped" else "kept", "\n")

# Over-dispersion check: Pearson chi-square / residual df ~ 1 for a Poisson
disp <- sum(residuals(m_final, type = "pearson")^2) / m_final$df.residual
cat("Dispersion estimate:", round(disp, 3), "\n")

# --- 3. Out-of-sample performance ----------------------------------------------
test[, mu_null  := predict(m_null,  newdata = test, type = "response")]
test[, mu_final := predict(m_final, newdata = test, type = "response")]
perf <- data.table(
  model = c("Intercept only", "Final GLM"),
  n_parameters = c(length(coef(m_null)), length(coef(m_final))),
  AIC_train = round(c(AIC(m_null), AIC(m_final))),
  poisson_deviance_test = round(c(poisson_deviance(test$ClaimNb, test$mu_null),
                                  poisson_deviance(test$ClaimNb, test$mu_final)), 5)
)
perf[, deviance_reduction_pct := round(100 * (1 - poisson_deviance_test / poisson_deviance_test[1]), 2)]
print(perf)
fwrite(perf, "output/tables/frequency_model_performance.csv")

# --- 4. Relativities -----------------------------------------------------------
# The data are over-dispersed (dispersion > 1), so Poisson standard errors are
# too narrow. Quasi-Poisson keeps the same estimates and scales the SEs by
# sqrt(dispersion), which gives more honest confidence intervals.
co <- summary(m_final, dispersion = disp)$coefficients
rel <- data.table(term = rownames(co), estimate = co[, 1], se = co[, 2], p_value = co[, 4])
rel[, `:=`(relativity = exp(estimate), lower95 = exp(estimate - 1.96 * se),
           upper95 = exp(estimate + 1.96 * se))]
fwrite(rel[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
       "output/tables/frequency_relativities.csv")

plot_rel <- function(prefix, base, title, order = NULL) {
  r <- rel[startsWith(term, prefix)]
  r[, level := sub(prefix, "", term, fixed = TRUE)]
  r <- rbind(r, data.table(level = base, relativity = 1, lower95 = 1, upper95 = 1), fill = TRUE)
  if (!is.null(order)) r[, level := factor(level, levels = order)]
  ggplot(r, aes(level, relativity)) +
    geom_hline(yintercept = 1, colour = "grey60") +
    geom_pointrange(aes(ymin = lower95, ymax = upper95), colour = "#08519c") +
    labs(title = title, x = NULL, y = "Frequency relativity (base = 1)") +
    theme_minimal(base_size = 11)
}
p_age <- plot_rel("DrivAgeBand", "(40,50]", "Driver age relativities (95% CI)",
                  levels(cut(1, c(17, 20, 25, 30, 40, 50, 60, 70, 80, Inf))))
ggsave("output/figures/relativities_driver_age.png", p_age, width = 8, height = 4.5, dpi = 150)
p_veh <- plot_rel("VehAgeBand", "(5,10]", "Vehicle age relativities (95% CI)",
                  levels(cut(1, c(-1, 0, 2, 5, 10, 15, Inf))))
ggsave("output/figures/relativities_vehicle_age.png", p_veh, width = 8, height = 4.5, dpi = 150)

# --- 5. Lift chart on the test set ---------------------------------------------
test[, decile := cut(frank(mu_final / Exposure, ties.method = "random"),
                     10, labels = 1:10)]
lift <- test[, .(observed = sum(ClaimNb) / sum(Exposure),
                 predicted = sum(mu_final) / sum(Exposure)), by = decile][order(decile)]
fwrite(lift, "output/tables/frequency_lift.csv")
p_lift <- ggplot(melt(lift, id.vars = "decile"), aes(decile, value, colour = variable, group = variable)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c(observed = "#cb181d", predicted = "#08519c")) +
  labs(title = "Test-set lift chart: policies ranked by predicted frequency",
       x = "Predicted frequency decile (1 = lowest risk)", y = "Claims per policy-year",
       colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/frequency_lift_chart.png", p_lift, width = 8, height = 5, dpi = 150)
cat("Top/bottom decile observed frequency ratio:",
    round(lift$observed[10] / lift$observed[1], 2), "\n")

dir.create("output/models", showWarnings = FALSE)
saveRDS(m_final, "output/models/frequency_model.rds")
