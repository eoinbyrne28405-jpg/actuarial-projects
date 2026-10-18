# =============================================================================
# 02_gompertz.R
# Fit Gompertz and Gompertz-Makeham laws of mortality to adult Irish mortality.
#
#   Gompertz:  mu_x = B * c^x            (ageing: hazard grows exponentially)
#   Makeham:   mu_x = A + B * c^x        (adds an age-independent accident term)
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
lt <- readRDS("output/life_tables.rds")

# Force of mortality from qx (constant force within each year of age):
#   mu_{x+1/2} ~ -log(1 - q_x)
fit_ages <- 40:90

# The life table is a synthetic cohort of 100,000 births, so we can treat
# d_x deaths out of L_x person-years of exposure as Poisson data:
#   D_x ~ Poisson(E_x * mu_x)
# This is the standard way actuaries graduate mortality rates. Maximum
# likelihood automatically gives more weight to ages with more deaths.
fit_one <- function(tab, sex) {
  d   <- tab[tab$age %in% fit_ages, ]
  d$x <- d$age + 0.5
  d$mu <- d$dx / d$Lx                       # central death rate ~ mu_{x+1/2}

  # Gompertz: log mu = log B + x log c, a Poisson GLM with log(exposure) offset
  g  <- glm(round(dx) ~ x, family = poisson(), offset = log(Lx), data = d)
  B  <- exp(coef(g)[1]); cc <- exp(coef(g)[2])

  # Makeham: maximise the Poisson log-likelihood directly
  nll <- function(p) {
    mu <- exp(p[1]) + exp(p[2] + p[3] * d$x)
    -sum(d$dx * log(mu) - d$Lx * mu)
  }
  opt <- optim(unname(c(log(5e-4), coef(g)[1], coef(g)[2])), nll, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-12))
  pm  <- c(A = exp(opt$par[1]), B = exp(opt$par[2]), k = opt$par[3])

  d$gompertz <- B * cc^d$x
  d$makeham  <- pm["A"] + pm["B"] * exp(pm["k"] * d$x)
  d$sex <- sex

  rmse <- function(fit) sqrt(mean((log(fit) - log(d$mu))^2))
  params <- data.frame(
    sex = sex,
    gompertz_B = signif(B, 4), gompertz_c = round(cc, 5),
    mortality_doubling_years = round(log(2) / log(cc), 1),
    makeham_A = signif(pm["A"], 4), makeham_B = signif(pm["B"], 4),
    makeham_c = round(exp(pm["k"]), 5),
    rmse_log_gompertz = round(rmse(d$gompertz), 4),
    rmse_log_makeham  = round(rmse(d$makeham), 4))
  list(data = d, params = params)
}

fm <- fit_one(lt$male, "Male")
ff <- fit_one(lt$female, "Female")
params <- rbind(fm$params, ff$params)
print(params, row.names = FALSE)
write.csv(params, "output/tables/gompertz_makeham_parameters.csv", row.names = FALSE)

d <- rbind(fm$data, ff$data)
pal <- c(Male = "#2171b5", Female = "#cb181d")
p <- ggplot(d, aes(x, mu, colour = sex)) +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_line(aes(y = gompertz), linewidth = 0.8) +
  geom_line(aes(y = makeham), linewidth = 0.6, linetype = "dashed") +
  scale_y_log10() +
  scale_colour_manual(values = pal) +
  labs(title = "Gompertz (solid) and Makeham (dashed) fits to Irish adult mortality",
       subtitle = sprintf("Mortality doubles roughly every %.1f years (males) and %.1f years (females)",
                          fm$params$mortality_doubling_years, ff$params$mortality_doubling_years),
       x = "Age", y = expression("Force of mortality "*mu[x]*" (log scale)"), colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/gompertz_fit.png", p, width = 8, height = 5, dpi = 150)

resid <- rbind(transform(d, model = "Gompertz", r = log(mu) - log(gompertz)),
               transform(d, model = "Makeham",  r = log(mu) - log(makeham)))
res <- ggplot(resid, aes(x, r, colour = sex)) +
  facet_wrap(~ model) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = 1.4) +
  scale_colour_manual(values = pal) +
  labs(title = "Residuals (log observed - log fitted)",
       subtitle = "Gompertz misses the extra, age-independent mortality at younger adult ages",
       x = "Age", y = "Log residual", colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/gompertz_residuals.png", res, width = 8, height = 5, dpi = 150)
