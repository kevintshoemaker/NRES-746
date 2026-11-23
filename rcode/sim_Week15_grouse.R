#  Simulates the sage-grouse trend series used on the Week 15 page (dynamic models).
#  Writes data/SageGrouseTrend.csv: 30 years of spring male counts for one population.

# Truth: Gompertz state-space model on the log scale -----------------------------------------
#   x[t] = a + b * x[t-1] + c * w[t] + process noise,   log(count[t]) = x[t] + observation noise

truth <- c(a = 0.3 * log(300), b = 0.7, c = 0.15, sigma_p = 0.15, sigma_o = 0.20)
n_yr <- 30
if (!exists("seed")) seed <- 18
if (!exists("write_csv")) write_csv <- TRUE       # FALSE returns the data frame without writing

set.seed(seed)
precip_z <- round(as.numeric(scale(rnorm(n_yr))), 2)    # spring precipitation the year before (SD units)
x <- numeric(n_yr)
x_prev <- rnorm(1, truth["a"] / (1 - truth["b"]), truth["sigma_p"] / sqrt(1 - truth["b"]^2))
for (t in seq_len(n_yr)) {
  x[t] <- truth["a"] + truth["b"] * x_prev + truth["c"] * precip_z[t] + rnorm(1, 0, truth["sigma_p"])
  x_prev <- x[t]
}
males <- round(exp(x + rnorm(n_yr, 0, truth["sigma_o"])))

grouse <- data.frame(year = 1996:2025, males = males, precip_z = precip_z, log_true = round(x, 3))
if (write_csv) write.csv(grouse, "data/SageGrouseTrend.csv", row.names = FALSE)
