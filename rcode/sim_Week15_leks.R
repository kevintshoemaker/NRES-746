#  Simulates the lek survey used in the Week 15 demo (a spatiotemporal count model).
#  Writes data/LekCounts.csv: 60 leks counted each spring, 2011-2025.

# Truth -------------------------------------------------------------------------------------
#   log(mu[s,t]) = b0 + b_sage * sage[s] + b_precip * precip[s,t] + omega[s] + eps[s,t]
#   eps[,t] = rho * eps[,t-1] + spatially correlated innovations;  count ~ NegBin(mu, size)

truth <- c(b0 = log(25), b_sage = 0.4, b_precip = 0.25, sigma_omega = 0.5, ell_omega = 12,
           sigma_eps = 0.4, ell_eps = 15, rho = 0.7, size = 8)
years <- 2011:2025
n_lek <- 60
if (!exists("seed")) seed <- 1
if (!exists("write_csv")) write_csv <- TRUE

set.seed(seed)
leks <- data.frame(lek = sprintf("L%02d", seq_len(n_lek)),
                   x_km = round(runif(n_lek, 2, 98), 1), y_km = round(runif(n_lek, 2, 78), 1))
D <- as.matrix(dist(leks[, c("x_km", "y_km")]))
draw_field <- function(sd, ell) as.vector(t(chol(sd^2 * exp(-D / ell) + diag(1e-8, n_lek))) %*% rnorm(n_lek))

# sagebrush cover (%): higher in the north and east, plus local variation
leks$sagebrush_pct <- round(pmin(45, pmax(3, 12 + 0.10 * leks$x_km + 0.12 * leks$y_km + draw_field(5, 10))))

# precipitation the spring before each count (SD units): a regional year effect plus storm tracks
precip <- sapply(years, function(yr) rnorm(1, 0, 0.8) + draw_field(0.6, 30))
precip <- (precip - mean(precip)) / sd(precip)

# fields: a permanent one, and one that changes from year to year with AR(1) memory
omega <- draw_field(truth["sigma_omega"], truth["ell_omega"])
eps <- matrix(0, n_lek, length(years))
eps[, 1] <- draw_field(truth["sigma_eps"], truth["ell_eps"])
for (t in 2:length(years)) {
  eps[, t] <- truth["rho"] * eps[, t - 1] +
    draw_field(truth["sigma_eps"] * sqrt(1 - truth["rho"]^2), truth["ell_eps"])
}

lek_counts <- do.call(rbind, lapply(seq_along(years), function(t) {
  log_mu <- truth["b0"] + truth["b_sage"] * (leks$sagebrush_pct - 20) / 10 + truth["b_precip"] * precip[, t] +
    omega + eps[, t]
  data.frame(leks, year = years[t], precip_z = round(precip[, t], 2),
             males = rnbinom(n_lek, mu = exp(log_mu), size = truth["size"]), eps_true = round(eps[, t], 3))
}))
if (write_csv) write.csv(lek_counts, "data/LekCounts.csv", row.names = FALSE)
