#  Simulates the bighorn survey used in the Week 15 take-home WOD (Filter by Hand).
#  Writes data/BighornCounts.csv: 25 years of helicopter counts of ewes in one mountain range.

# Truth: Gompertz state-space model on the log scale -----------------------------------------
#   x[t] = a + b * x[t-1] + c * w[t] + process noise,   log(count[t]) = x[t] + observation noise

truth <- c(a = 0.2 * log(120), b = 0.8, c = 0.12, sigma_p = 0.10, sigma_o = 0.15)
n_yr <- 25
if (!exists("seed")) seed <- 16
if (!exists("write_csv")) write_csv <- TRUE

set.seed(seed)
precip_z <- round(as.numeric(scale(rnorm(n_yr))), 2)    # winter precipitation before the survey (SD units)
x <- numeric(n_yr)
x_prev <- rnorm(1, truth["a"] / (1 - truth["b"]), truth["sigma_p"] / sqrt(1 - truth["b"]^2))
for (t in seq_len(n_yr)) {
  x[t] <- truth["a"] + truth["b"] * x_prev + truth["c"] * precip_z[t] + rnorm(1, 0, truth["sigma_p"])
  x_prev <- x[t]
}
ewes <- round(exp(x + rnorm(n_yr, 0, truth["sigma_o"])))

bighorn <- data.frame(year = 2001:2025, ewes = ewes, precip_z = precip_z)
if (write_csv) write.csv(bighorn, "data/BighornCounts.csv", row.names = FALSE)
