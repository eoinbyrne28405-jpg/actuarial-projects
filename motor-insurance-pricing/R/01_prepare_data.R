# =============================================================================
# 01_prepare_data.R
# Clean the French MTPL data and build rating factors.
#
# Data: freMTPL2freq (677,991 policies) and freMTPL2sev
# (26,444 individual claim amounts) from the CASdatasets R package
# (Charpentier, Dutang). French motor third-party liability, one year.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

load("data/freMTPL2freq.rda")
load("data/freMTPL2sev.rda")
freq <- as.data.table(freMTPL2freq)
sev  <- as.data.table(freMTPL2sev)
freq[, IDpol := as.integer(as.character(IDpol))]
sev[,  IDpol := as.integer(as.character(IDpol))]

# --- Data quality fixes (standard for this dataset; see Noll, Salzmann & Wuthrich 2018)
# A handful of policies have implausible claim counts or exposure > 1 year.
freq[, ClaimNb  := pmin(ClaimNb, 4)]
freq[, Exposure := pmin(Exposure, 1)]

# Total claim cost per policy, joined from the severity file
cost <- sev[, .(ClaimTotal = sum(ClaimAmount), nSev = .N), by = IDpol]
freq <- merge(freq, cost, by = "IDpol", all.x = TRUE)
freq[is.na(ClaimTotal), `:=`(ClaimTotal = 0, nSev = 0L)]
# Sanity check: every policy with claims should have amounts in the severity file
cat("Policies with claims but no recorded amount:", freq[ClaimNb > 0 & nSev == 0, .N], "\n")

# --- Rating factors -----------------------------------------------------------
# GLMs are linear in the link scale, so continuous variables with non-linear
# effects are banded. Bands follow the shape of the one-way plots.
freq[, DrivAgeBand := cut(DrivAge, c(17, 20, 25, 30, 40, 50, 60, 70, 80, Inf),
                          right = TRUE)]
freq[, VehAgeBand  := cut(VehAge, c(-1, 0, 2, 5, 10, 15, Inf))]
freq[, VehPowerF   := factor(pmin(VehPower, 12))]
freq[, BonusMalusC := pmin(BonusMalus, 150)]           # cap the long tail
freq[, LogDensity  := log(Density)]
freq[, VehBrand    := relevel(VehBrand, ref = "B1")]
freq[, Region      := relevel(Region, ref = "Centre")]
freq[, VehGas      := factor(VehGas)]

# Choose base levels with the most exposure (makes relativities easier to read)
freq[, DrivAgeBand := relevel(DrivAgeBand, ref = "(40,50]")]
freq[, VehAgeBand  := relevel(VehAgeBand,  ref = "(5,10]")]
freq[, VehPowerF   := relevel(VehPowerF,   ref = "6")]

# --- Train / test split (80/20, fixed seed) ------------------------------------
set.seed(2026)
freq[, split := ifelse(runif(.N) < 0.8, "train", "test")]

saveRDS(freq, "data/model_data.rds")
cat("Rows:", nrow(freq), " Train:", freq[split == "train", .N],
    " Test:", freq[split == "test", .N], "\n")
cat("Overall claim frequency:", round(sum(freq$ClaimNb) / sum(freq$Exposure), 4), "per year\n")
