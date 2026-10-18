# =============================================================================
# 03_life_insurance_values.R
# Expected present values (EPVs) of life insurance and annuity benefits using
# the Irish life table, plus net premiums and sensitivity tests.
#
# Notation (standard actuarial):
#   v = 1/(1+i),  d = i/(1+i)
#   k p_x   probability a life aged x survives k years
#   ä_x     whole life annuity-due of 1 per year         = sum_k v^k  k p_x
#   A_x     whole life assurance of 1 at end of year of death
#           = sum_k v^(k+1) k p_x q_{x+k}
#   A1_x:n  n-year term assurance; ä_x:n n-year temporary annuity-due
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
lt <- readRDS("output/life_tables.rds")

# Survival and death probabilities for a life aged x, from a vector of qx.
# improvement: annual % reduction in mortality for future years (0 = use the
# period table as is, which assumes no future improvement).
cohort_q <- function(qx_all, ages, x, improvement = 0) {
  q <- qx_all[ages >= x]
  q * (1 - improvement)^(seq_along(q) - 1)
}

epv <- function(qx_all, ages, x, i, n = Inf, improvement = 0) {
  q  <- cohort_q(qx_all, ages, x, improvement)
  q[length(q)] <- 1
  K  <- length(q)
  kp <- c(1, cumprod(1 - q))[1:K]          # k p_x for k = 0..K-1
  v  <- 1 / (1 + i)
  k  <- 0:(K - 1)
  in_term <- k < n
  list(
    annuity_due = sum(v^k * kp * in_term),                # ä_x (or ä_x:n)
    assurance   = sum(v^(k + 1) * kp * q * in_term),      # A_x (or term A1_x:n)
    pure_endow  = if (is.finite(n)) v^n * c(kp, 0)[n + 1] else 0
  )
}

ages <- lt$male$age
qm <- lt$male$qx; qf <- lt$female$qx

# --- 1. Check: A_x + d * ä_x = 1 (holds exactly for whole life) ---------------
i <- 0.03; d <- i / (1 + i)
e <- epv(qm, ages, 40, i)
cat(sprintf("Check at 40, male, 3%%: A + d*ä = %.10f\n", e$assurance + d * e$annuity_due))
stopifnot(abs(e$assurance + d * e$annuity_due - 1) < 1e-10)

# --- 2. Table of whole-life values by age ----------------------------------------
tab <- do.call(rbind, lapply(seq(20, 90, 5), function(x) {
  m <- epv(qm, ages, x, i); f <- epv(qf, ages, x, i)
  data.frame(age = x, annuity_male = m$annuity_due, annuity_female = f$annuity_due,
             assurance_male = m$assurance, assurance_female = f$assurance)
}))
cat("\nWhole life values at 3% interest:\n"); print(round(tab, 4), row.names = FALSE)
write.csv(round(tab, 5), "output/tables/whole_life_values_3pct.csv", row.names = FALSE)

# --- 3. Worked examples ------------------------------------------------------------
# (a) 20-year term life cover of EUR 250,000 bought at 30: level annual net premium
#     P = SA * A1_{30:20} / ä_{30:20}
SA <- 250000
term_prem <- sapply(list(Male = qm, Female = qf), function(q) {
  e <- epv(q, ages, 30, i, n = 20); SA * e$assurance / e$annuity_due })

# (b) Cost at 65 of a pension of EUR 10,000 a year for life, paid annually in advance
pension_cost <- sapply(list(Male = qm, Female = qf), function(q)
  10000 * epv(q, ages, 65, i)$annuity_due)

# (c) Same pension if mortality keeps improving by 1.5% a year (longevity risk)
pension_improved <- sapply(list(Male = qm, Female = qf), function(q)
  10000 * epv(q, ages, 65, i, improvement = 0.015)$annuity_due)

examples <- data.frame(
  example = c("20-yr term assurance, EUR 250k, age 30: annual net premium",
              "Life pension EUR 10k p.a. from 65: single premium (period table)",
              "Same pension with 1.5% p.a. mortality improvement",
              "Increase in pension cost from improvement"),
  male   = c(round(term_prem["Male"]), round(pension_cost["Male"]),
             round(pension_improved["Male"]),
             paste0(round(100 * (pension_improved["Male"] / pension_cost["Male"] - 1), 1), "%")),
  female = c(round(term_prem["Female"]), round(pension_cost["Female"]),
             round(pension_improved["Female"]),
             paste0(round(100 * (pension_improved["Female"] / pension_cost["Female"] - 1), 1), "%")))
cat("\nWorked examples (3% interest):\n"); print(examples, row.names = FALSE)
write.csv(examples, "output/tables/worked_examples.csv", row.names = FALSE)

# --- 4. Interest rate sensitivity ------------------------------------------------
rates <- c(0.01, 0.02, 0.03, 0.04, 0.05)
sens <- data.frame(interest = paste0(100 * rates, "%"),
  pension_cost_male   = sapply(rates, function(r) round(10000 * epv(qm, ages, 65, r)$annuity_due)),
  pension_cost_female = sapply(rates, function(r) round(10000 * epv(qf, ages, 65, r)$annuity_due)),
  term_premium_male   = sapply(rates, function(r) { e <- epv(qm, ages, 30, r, n = 20)
                                                    round(SA * e$assurance / e$annuity_due) }))
cat("\nInterest rate sensitivity:\n"); print(sens, row.names = FALSE)
write.csv(sens, "output/tables/interest_sensitivity.csv", row.names = FALSE)

# --- 5. Plots ------------------------------------------------------------------------
grid <- expand.grid(age = 20:95, rate = c(0.01, 0.03, 0.05), sex = c("Male", "Female"),
                    stringsAsFactors = FALSE)
grid$annuity <- mapply(function(x, r, s) epv(if (s == "Male") qm else qf, ages, x, r)$annuity_due,
                       grid$age, grid$rate, grid$sex)
grid$rate <- factor(paste0(100 * grid$rate, "% interest"),
                    levels = paste0(c(1, 3, 5), "% interest"))
p1 <- ggplot(grid, aes(age, annuity, colour = sex, linetype = rate)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(Male = "#2171b5", Female = "#cb181d")) +
  labs(title = "Value of a whole-life annuity of 1 per year (annuity-due)",
       subtitle = "Cost of a lifetime income falls with age and with higher interest rates",
       x = "Age at purchase", y = expression(ddot(a)[x]), colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/annuity_values.png", p1, width = 8, height = 5, dpi = 150)

p2 <- ggplot(data.frame(rate = 100 * rates, cost = sens$pension_cost_male, sex = "Male") |>
               rbind(data.frame(rate = 100 * rates, cost = sens$pension_cost_female, sex = "Female")),
             aes(rate, cost / 1000, colour = sex)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c(Male = "#2171b5", Female = "#cb181d")) +
  labs(title = "Single premium at 65 for a EUR 10,000 a year lifetime pension",
       x = "Interest rate (%)", y = "Cost (EUR thousands)", colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/pension_cost_vs_interest.png", p2, width = 8, height = 5, dpi = 150)
