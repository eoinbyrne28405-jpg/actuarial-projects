# =============================================================================
# 01_life_table.R
# Build a full life table from Irish Life Table No. 17 (2015-2017) and
# analyse survival and life expectancy by sex.
#
# Data: CSO Irish Life Tables No. 17, 2015-2017, single-year qx and lx for
# males and females, ages 0-105. Source document: CSO, via the Human Life-Table
# Database (lifetable.de). See data/README.md.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

ilt <- read.csv("data/irish_life_table_17.csv")

# Check the published lx are consistent with qx: l_{x+1} = l_x (1 - q_x)
for (s in c("male", "female")) {
  l <- ilt[[paste0("lx_", s)]]; q <- ilt[[paste0("qx_", s)]]
  stopifnot(max(abs(head(l, -1) * (1 - head(q, -1)) - l[-1])) < 1.5)
}

# Build a complete life table for one sex from qx alone.
# Deaths are assumed uniform over each year of age (UDD), so a person dying at
# age x lives on average half a year, except at age 0 where infant deaths are
# concentrated in the first weeks (a0 = 0.1 is a standard assumption).
build_life_table <- function(age, qx, radix = 100000) {
  qx[length(qx)] <- 1                       # everyone dies in the final age group
  n  <- length(qx)
  lx <- radix * c(1, cumprod(1 - qx)[-n])
  dx <- lx * qx
  ax <- c(0.1, rep(0.5, n - 1))
  Lx <- lx - (1 - ax) * dx                  # person-years lived between x and x+1
  Lx[n] <- 0.5 * lx[n]                     # final age: assume death mid-year
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx                             # complete expectation of life
  data.frame(age, qx, px = 1 - qx, lx, dx, Lx, Tx, ex)
}

lt_m <- build_life_table(ilt$age, ilt$qx_male)
lt_f <- build_life_table(ilt$age, ilt$qx_female)
write.csv(round(lt_m, 6), "output/tables/life_table_male.csv",   row.names = FALSE)
write.csv(round(lt_f, 6), "output/tables/life_table_female.csv", row.names = FALSE)

key_ages <- c(0, 20, 40, 65, 80)
ex_tbl <- data.frame(age = key_ages,
                     male = round(lt_m$ex[match(key_ages, lt_m$age)], 1),
                     female = round(lt_f$ex[match(key_ages, lt_f$age)], 1))
ex_tbl$gap <- ex_tbl$female - ex_tbl$male
cat("Life expectancy (years):\n"); print(ex_tbl, row.names = FALSE)
write.csv(ex_tbl, "output/tables/life_expectancy_key_ages.csv", row.names = FALSE)

# Survival probabilities from 65 and from birth
tpx <- function(lt, x, t) lt$lx[lt$age == x + t] / lt$lx[lt$age == x]
surv_tbl <- data.frame(
  measure = c("Newborn reaches 65", "Newborn reaches 85", "65-year-old reaches 85",
              "65-year-old reaches 95", "Median age at death (from birth)"),
  male = c(tpx(lt_m, 0, 65), tpx(lt_m, 0, 85), tpx(lt_m, 65, 20), tpx(lt_m, 65, 30),
           lt_m$age[which(lt_m$lx < 50000)[1]] - 1),
  female = c(tpx(lt_f, 0, 65), tpx(lt_f, 0, 85), tpx(lt_f, 65, 20), tpx(lt_f, 65, 30),
             lt_f$age[which(lt_f$lx < 50000)[1]] - 1))
surv_tbl[1:4, 2:3] <- round(surv_tbl[1:4, 2:3], 3)
cat("\nSurvival:\n"); print(surv_tbl, row.names = FALSE)
write.csv(surv_tbl, "output/tables/survival_probabilities.csv", row.names = FALSE)

# Modal age at death (peak of the d_x curve, excluding infancy)
modal <- function(lt) lt$age[lt$age > 10][which.max(lt$dx[lt$age > 10])]
cat("\nModal age at death: male", modal(lt_m), " female", modal(lt_f), "\n")

# --- Plots ---------------------------------------------------------------------
both <- rbind(cbind(lt_m, sex = "Male"), cbind(lt_f, sex = "Female"))
pal  <- c(Male = "#2171b5", Female = "#cb181d")

p1 <- ggplot(both[both$age < 105, ], aes(age, qx, colour = sex)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = pal) +
  labs(title = "Irish mortality rates by age, 2015-2017 (log scale)",
       subtitle = "High infant mortality, a low point around age 10, a sharp rise in the late teens, then near-exponential growth",
       x = "Age", y = expression(q[x]~"(probability of death within a year)"), colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/mortality_rates_log.png", p1, width = 8, height = 5, dpi = 150)

p2 <- ggplot(both, aes(age, lx / 1000, colour = sex)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Survival curve: of 100,000 births, how many are still alive at each age",
       x = "Age", y = "Survivors (thousands)", colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/survival_curves.png", p2, width = 8, height = 5, dpi = 150)

p3 <- ggplot(both[both$age > 0, ], aes(age, dx, colour = sex)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Distribution of age at death (curve of deaths)",
       x = "Age", y = "Deaths per 100,000 births", colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/curve_of_deaths.png", p3, width = 8, height = 5, dpi = 150)

saveRDS(list(male = lt_m, female = lt_f), "output/life_tables.rds")
