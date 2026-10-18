# =============================================================================
# 01_chain_ladder.R
# Deterministic chain ladder on the RAA cumulative claims triangle.
#
# Data: RAA - Reinsurance Association of America, general liability,
#       accident years 1981-1990, cumulative incurred claims ($000s).
#       Shipped with the ChainLadder package (Mack, 1993 uses this triangle).
# =============================================================================

suppressPackageStartupMessages({
  library(ChainLadder)
  library(ggplot2)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

tri <- RAA
n   <- nrow(tri)
cat("RAA cumulative triangle ($000s):\n")
print(tri)

# --- 1. Age-to-age (link) ratios ---------------------------------------------
link_ratios <- ata(tri)
cat("\nIndividual age-to-age factors:\n")
print(round(link_ratios, 3))

# --- 2. Volume-weighted development factors (implemented by hand) -------------
# f_j = sum_i C_{i,j+1} / sum_i C_{i,j}, over accident years where both exist.
dev_factors <- sapply(1:(n - 1), function(j) {
  rows <- 1:(n - j)
  sum(tri[rows, j + 1]) / sum(tri[rows, j])
})
names(dev_factors) <- paste0(1:(n - 1), "-", 2:n)

# Cross-check against the package's volume-weighted average
stopifnot(all.equal(unname(dev_factors), unname(attr(link_ratios, "vwtd"))))

cdf <- rev(cumprod(rev(dev_factors)))   # cumulative development factor to ultimate
cat("\nVolume-weighted development factors:\n"); print(round(dev_factors, 4))
cat("\nCumulative development factors (to ultimate):\n"); print(round(cdf, 4))

# --- 3. Project the lower triangle --------------------------------------------
full <- tri
for (j in 1:(n - 1)) {
  missing <- is.na(full[, j + 1])
  full[missing, j + 1] <- full[missing, j] * dev_factors[j]
}

latest   <- getLatestCumulative(tri)
ultimate <- full[, n]
reserve  <- ultimate - latest

# CDF applying to each accident year: AY 1981 is fully developed (CDF 1),
# AY 1990 has only its first development year observed (largest CDF).
ay_cdf <- c(1, rev(cdf))

summary_tbl <- data.frame(
  accident_year = as.integer(rownames(tri)),
  latest        = round(latest),
  dev_to_ult    = round(ay_cdf, 4),
  ultimate      = round(ultimate),
  reserve       = round(reserve)
)

summary_tbl$accident_year <- as.character(summary_tbl$accident_year)
total <- data.frame(accident_year = "Total", latest = sum(summary_tbl$latest),
                    dev_to_ult = NA, ultimate = sum(summary_tbl$ultimate),
                    reserve = sum(summary_tbl$reserve))
summary_out <- rbind(summary_tbl, total)

cat("\nChain ladder reserve summary:\n"); print(summary_out, row.names = FALSE)
write.csv(summary_out, "output/tables/chain_ladder_summary.csv", row.names = FALSE, na = "")
write.csv(round(as.matrix(full)), "output/tables/completed_triangle.csv")
write.csv(data.frame(period = names(dev_factors), factor = round(dev_factors, 4),
                     cdf = round(cdf, 4)),
          "output/tables/development_factors.csv", row.names = FALSE)

# --- 4. Plots -----------------------------------------------------------------
long <- do.call(rbind, lapply(1:n, function(i) data.frame(
  ay = rownames(tri)[i], dev = 1:n, value = as.numeric(full[i, ]),
  type = ifelse(is.na(tri[i, ]), "Projected", "Observed"))))

p1 <- ggplot(long, aes(dev, value / 1000, group = ay, colour = ay)) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = type), size = 1.8) +
  scale_shape_manual(values = c(Observed = 16, Projected = 1)) +
  scale_colour_viridis_d(end = 0.92) +
  scale_x_continuous(breaks = 1:n) +
  labs(title = "RAA claims development: observed and chain ladder projections",
       x = "Development year", y = "Cumulative claims ($m)",
       colour = "Accident year", shape = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/figures/development_curves.png", p1, width = 8, height = 5, dpi = 150)

heat <- do.call(rbind, lapply(1:(n - 1), function(i) data.frame(
  ay = rownames(link_ratios)[i], period = colnames(link_ratios),
  f = as.numeric(link_ratios[i, ]))))
heat <- heat[!is.na(heat$f), ]
heat$period <- factor(heat$period, levels = colnames(link_ratios))
p2 <- ggplot(heat, aes(period, ay, fill = log(f))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", f)), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#2171b5", guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(title = "Individual age-to-age factors", x = "Development period",
       y = "Accident year") +
  theme_minimal(base_size = 11)
ggsave("output/figures/link_ratio_heatmap.png", p2, width = 8, height = 5, dpi = 150)

bars <- rbind(
  data.frame(ay = summary_tbl$accident_year, part = "Emerged to date", value = summary_tbl$latest),
  data.frame(ay = summary_tbl$accident_year, part = "Reserve (still to emerge)", value = summary_tbl$reserve))
bars$part <- factor(bars$part, levels = c("Reserve (still to emerge)", "Emerged to date"))
p3 <- ggplot(bars, aes(factor(ay), value / 1000, fill = part)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Emerged to date" = "#2171b5",
                               "Reserve (still to emerge)" = "#fd8d3c")) +
  labs(title = "Chain ladder ultimate claims by accident year",
       x = "Accident year", y = "$m", fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("output/figures/ultimate_vs_latest.png", p3, width = 8, height = 5, dpi = 150)

cat("\nDone. Total chain ladder reserve: $", format(round(sum(reserve)), big.mark = ","), "k\n")
