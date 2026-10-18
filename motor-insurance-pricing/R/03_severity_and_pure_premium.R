# =============================================================================
# 03_severity_and_pure_premium.R
# Gamma GLM for average claim size, then pure premium = frequency x severity.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
d       <- readRDS("data/model_data.rds")
m_freq  <- readRDS("output/models/frequency_model.rds")
load("data/freMTPL2sev.rda")

# --- 1. Large claims ------------------------------------------------------------
# Claim sizes are extremely heavy-tailed (max 4m vs median ~1,170). A handful of
# large claims would dominate a gamma GLM, so individual claims are capped and
# the excess is spread across all policies as a flat loading.
amounts <- freMTPL2sev$ClaimAmount
cap     <- as.numeric(round(quantile(amounts, 0.995), -2))
excess_share <- sum(pmax(amounts - cap, 0)) / sum(amounts)
cat(sprintf("Cap at 99.5th percentile = %s; excess over cap = %.1f%% of total claim cost\n",
            format(cap, big.mark = ","), 100 * excess_share))

sev <- as.data.table(freMTPL2sev)[, IDpol := as.integer(as.character(IDpol))]
sev[, capped := pmin(ClaimAmount, cap)]
pol_sev <- sev[, .(capped_total = sum(capped), n = .N), by = IDpol]
s <- merge(d, pol_sev, by = "IDpol")
s[, avg_sev := capped_total / n]

p_sev <- ggplot(sev, aes(ClaimAmount)) +
  geom_histogram(bins = 80, fill = "#6baed6", colour = "white") +
  geom_vline(xintercept = cap, colour = "#cb181d") +
  scale_x_log10(labels = scales::label_comma()) +
  annotate("text", x = cap, y = Inf, vjust = 1.5, hjust = -0.05, colour = "#cb181d",
           label = paste0("Cap: ", format(cap, big.mark = ",")), size = 3.3) +
  labs(title = "Individual claim amounts (log scale)", x = "Claim amount (EUR)", y = "Count") +
  theme_minimal(base_size = 11)
ggsave("output/figures/claim_size_distribution.png", p_sev, width = 8, height = 5, dpi = 150)

# --- 2. Gamma GLM (log link), weighted by number of claims ----------------------
# Average of n claims has variance phi*mu^2/n, hence weights = n.
f_sev <- avg_sev ~ DrivAgeBand + VehAgeBand + VehGas + BonusMalusC + LogDensity + VehBrand
m_sev_full <- glm(f_sev, family = Gamma(link = "log"), data = s[split == "train"], weights = n)
# Drop terms that don't earn their place (F-test, appropriate for gamma)
print(drop1(m_sev_full, test = "F"))
keep <- drop1(m_sev_full, test = "F")
sig  <- rownames(keep)[!is.na(keep$`Pr(>F)`) & keep$`Pr(>F)` < 0.05]
m_sev <- glm(reformulate(sig, "avg_sev"), family = Gamma(link = "log"),
             data = s[split == "train"], weights = n)
cat("Severity terms kept:", paste(sig, collapse = ", "), "\n")

co <- summary(m_sev)$coefficients
fwrite(data.table(term = rownames(co), relativity = signif(exp(co[, 1]), 4),
                  p_value = signif(co[, 4], 3)),
       "output/tables/severity_relativities.csv")

# --- 3. Pure premium ------------------------------------------------------------
# Expected cost per policy-year =
#   frequency(x) * capped severity(x) * (1 + large-loss loading)
loading <- 1 / (1 - excess_share)
d[, freq_hat := predict(m_freq, newdata = d, type = "response") / Exposure]
d[, sev_hat  := predict(m_sev,  newdata = d, type = "response")]
d[, pure_premium := freq_hat * sev_hat * loading]

# Balance check: total predicted cost vs actual cost, train and test.
# Large claims are rare and lumpy, so a 20% test sample can easily contain
# more or fewer of them than average. Also check the capped (attritional) cost,
# which the GLMs model directly and should be much more stable.
sev_pol <- sev[, .(capped_cost = sum(capped)), by = IDpol]
d <- merge(d, sev_pol, by = "IDpol", all.x = TRUE)
d[is.na(capped_cost), capped_cost := 0]
bal <- d[, .(actual_total = sum(ClaimTotal),
             predicted_total = sum(pure_premium * Exposure),
             actual_capped = sum(capped_cost),
             predicted_capped = sum(pure_premium * Exposure) / loading), by = split]
bal[, `:=`(ratio_total = round(predicted_total / actual_total, 3),
           ratio_capped = round(predicted_capped / actual_capped, 3))]
print(bal)
fwrite(bal, "output/tables/balance_check.csv")
test <- d[split == "test"]

pp_summary <- d[, .(mean = mean(pure_premium), p10 = quantile(pure_premium, 0.1),
                    median = median(pure_premium), p90 = quantile(pure_premium, 0.9))]
cat("Pure premium distribution (EUR per policy-year):\n")
print(round(pp_summary, 1))

# Example customers
examples <- data.table(
  profile = c("Young driver, new powerful car, Paris area, BM 100",
              "Middle-aged driver, 8-year-old car, rural, BM 50",
              "Retired driver, 3-year-old car, mid-size town, BM 50"),
  DrivAge = c(19, 45, 72), DrivAgeBand = factor(c("(17,20]", "(40,50]", "(70,80]"), levels = levels(d$DrivAgeBand)),
  VehAgeBand = factor(c("(-1,0]", "(5,10]", "(2,5]"), levels = levels(d$VehAgeBand)),
  VehPowerF = factor(c("10", "6", "6"), levels = levels(d$VehPowerF)),
  VehBrand = factor("B1", levels = levels(d$VehBrand)),
  VehGas = factor(c("Regular", "Diesel", "Regular"), levels = levels(d$VehGas)),
  BonusMalusC = c(100, 50, 50), LogDensity = log(c(20000, 40, 800)),
  Region = factor(c("Ile-de-France", "Bretagne", "Rhone-Alpes"), levels = levels(d$Region)),
  Exposure = 1)
examples[, frequency := predict(m_freq, newdata = examples, type = "response")]
examples[, severity  := predict(m_sev,  newdata = examples, type = "response") * loading]
examples[, pure_premium := frequency * severity]
out <- examples[, .(profile, frequency = round(frequency, 4), severity = round(severity),
                    pure_premium = round(pure_premium))]
print(out)
fwrite(out, "output/tables/example_premiums.csv")

# --- 4. Plots -------------------------------------------------------------------
test[, decile := cut(frank(pure_premium, ties.method = "random"), 10, labels = 1:10)]
lift <- test[, .(observed = sum(ClaimTotal) / sum(Exposure),
                 predicted = sum(pure_premium * Exposure) / sum(Exposure)), by = decile][order(decile)]
fwrite(lift, "output/tables/pure_premium_lift.csv")
p_lift <- ggplot(melt(lift, id.vars = "decile"), aes(decile, value, fill = variable)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c(observed = "#fd8d3c", predicted = "#08519c")) +
  labs(title = "Test-set pure premium lift: actual vs predicted cost per policy-year",
       subtitle = "Observed cost is noisy because a few large claims fall in any decile",
       x = "Predicted pure premium decile (1 = lowest)", y = "EUR per policy-year", fill = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/pure_premium_lift.png", p_lift, width = 8, height = 5, dpi = 150)

p_pp <- ggplot(d, aes(pure_premium)) +
  geom_histogram(bins = 80, fill = "#9ecae1", colour = "white") +
  scale_x_log10(labels = scales::label_comma()) +
  labs(title = "Distribution of modelled pure premiums across the portfolio",
       x = "Pure premium (EUR per policy-year, log scale)", y = "Policies") +
  theme_minimal(base_size = 11)
ggsave("output/figures/pure_premium_distribution.png", p_pp, width = 8, height = 5, dpi = 150)
