# =============================================================================
# 02_mack_and_bootstrap.R
# Reserve uncertainty: Mack (1993) analytic standard errors and an
# over-dispersed Poisson (ODP) bootstrap of the full reserve distribution.
# =============================================================================

suppressPackageStartupMessages({
  library(ChainLadder)
  library(ggplot2)
})
set.seed(2026)

tri <- RAA

# --- 1. Mack chain ladder -----------------------------------------------------
# Distribution-free: assumes E[C_{i,j+1} | C_{i,j}] = f_j C_{i,j} and
# Var[C_{i,j+1} | C_{i,j}] = sigma_j^2 C_{i,j}, accident years independent.
mack <- MackChainLadder(tri, est.sigma = "Mack")
print(mack)

mack_tbl <- data.frame(
  accident_year = rownames(tri),
  latest   = round(getLatestCumulative(tri)),
  ultimate = round(mack$FullTriangle[, ncol(tri)]),
  reserve  = round(mack$FullTriangle[, ncol(tri)] - getLatestCumulative(tri)),
  mack_se  = round(mack$Mack.S.E[, ncol(tri)]),
  cv       = round(mack$Mack.S.E[, ncol(tri)] /
                   pmax(mack$FullTriangle[, ncol(tri)] - getLatestCumulative(tri), 1e-9), 3)
)
mack_tbl$cv[mack_tbl$reserve == 0] <- NA
total_reserve <- sum(mack_tbl$reserve)
total_se      <- as.numeric(mack$Total.Mack.S.E)
mack_tbl <- rbind(mack_tbl, data.frame(accident_year = "Total",
  latest = sum(mack_tbl$latest), ultimate = sum(mack_tbl$ultimate),
  reserve = total_reserve, mack_se = round(total_se),
  cv = round(total_se / total_reserve, 3)))
write.csv(mack_tbl, "output/tables/mack_summary.csv", row.names = FALSE, na = "")

# --- 2. Diagnostics: the Mack assumptions should look reasonable --------------
png("output/figures/mack_diagnostics.png", width = 1200, height = 900, res = 130)
plot(mack, lattice = FALSE)
dev.off()

# --- 3. ODP bootstrap ---------------------------------------------------------
# Resamples scaled Pearson residuals of an ODP GLM, refits the chain ladder on
# each pseudo-triangle and adds process variance (gamma), giving a full
# predictive distribution of the total reserve rather than just a mean and SE.
B <- 5000
boot <- BootChainLadder(tri, R = B, process.distr = "gamma")
print(boot)

sims <- as.numeric(boot$IBNR.Totals)
q    <- quantile(sims, c(0.5, 0.75, 0.9, 0.95, 0.995))
boot_tbl <- data.frame(
  statistic = c("Mean", "Std. deviation", "Median", "75th percentile",
                "90th percentile", "95th percentile", "99.5th percentile (Solvency II VaR)"),
  value = round(c(mean(sims), sd(sims), q))
)
write.csv(boot_tbl, "output/tables/bootstrap_reserve_distribution.csv", row.names = FALSE)
print(boot_tbl, row.names = FALSE)

p <- ggplot(data.frame(reserve = sims / 1000), aes(reserve)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#9ecae1", colour = "white") +
  geom_density(colour = "#08519c", linewidth = 0.8) +
  geom_vline(xintercept = total_reserve / 1000, linetype = "dashed") +
  geom_vline(xintercept = q["99.5%"] / 1000, colour = "#cb181d") +
  annotate("text", x = total_reserve / 1000, y = Inf, vjust = 1.5, hjust = -0.05,
           label = "Chain ladder\nbest estimate", size = 3.2) +
  annotate("text", x = q["99.5%"] / 1000, y = Inf, vjust = 1.5, hjust = -0.05,
           label = "99.5th\npercentile", size = 3.2, colour = "#cb181d") +
  labs(title = sprintf("Bootstrapped distribution of total RAA reserve (%s simulations)",
                       format(B, big.mark = ",")),
       x = "Total reserve ($m)", y = "Density") +
  theme_minimal(base_size = 11)
ggsave("output/figures/bootstrap_reserve_distribution.png", p, width = 8, height = 5, dpi = 150)

se_plot <- head(mack_tbl, -1)
se_plot <- se_plot[se_plot$reserve > 0, ]
p2 <- ggplot(se_plot, aes(accident_year, reserve / 1000)) +
  geom_col(fill = "#fd8d3c", width = 0.65) +
  geom_errorbar(aes(ymin = pmax(reserve - 1.96 * mack_se, 0) / 1000,
                    ymax = (reserve + 1.96 * mack_se) / 1000), width = 0.25) +
  labs(title = "Reserve by accident year with approximate 95% Mack intervals",
       x = "Accident year", y = "Reserve ($m)") +
  theme_minimal(base_size = 11)
ggsave("output/figures/mack_reserve_intervals.png", p2, width = 8, height = 5, dpi = 150)

cat(sprintf("\nTotal reserve %s, Mack SE %s (CV %.1f%%), bootstrap mean %s, 99.5%% %s\n",
            format(total_reserve, big.mark = ","), format(round(total_se), big.mark = ","),
            100 * total_se / total_reserve, format(round(mean(sims)), big.mark = ","),
            format(round(q["99.5%"]), big.mark = ",")))
