#  NRES 746, Week 6: Likelihood, maximum likelihood, and a power analysis in code ----
#  Code from Week6_LikelihoodAndPowerAnalysis.html. Base R only.

# LIKELIHOOD AND MAXIMUM LIKELIHOOD ESTIMATION =================================

# 2. WOD: writing the likelihood -----------------------------------------------
#  4 nests, outcomes 1, 1, 0, 1; L4(p) = p^3 (1 - p)

p_grid <- seq(0, 1, by = 0.005)
lik_nests <- function(p) p^3 * (1 - p)
p_three <- c(0.25, 0.5, 0.75)

par(mar = c(4, 4.5, 0.5, 0.5))
plot(p_grid, lik_nests(p_grid), type = "l", lwd = 2,
     xlab = "Nest success probability, p", ylab = "Likelihood")
points(p_three, lik_nests(p_three), pch = 16, cex = 1.2)
lik_nests(p_three) * 256   # 3, 16, 27 (in 256ths)
#  [PAUSE] add up the three likelihoods. do they sum to 1? should they?

# 3. From likelihood to negative log-likelihood --------------------------------

prod(rep(0.5, 1852))       # underflows to exactly 0
sum(log(rep(0.5, 1852)))   # the sum of logs is an ordinary number

haypiles <- c(2, 5, 3, 0, 4, 4)

nll_pois <- function(lambda, y) {
  -sum(dpois(y, lambda = lambda, log = TRUE))
}

lambda_grid <- seq(0.5, 7, by = 0.01)
lik_vals <- sapply(lambda_grid, function(l) prod(dpois(haypiles, l)))
nll_vals <- sapply(lambda_grid, nll_pois, y = haypiles)

par(mfrow = c(1, 2), mar = c(4, 4.5, 0.5, 0.5))
plot(lambda_grid, lik_vals, type = "l", lwd = 2,
     xlab = "lambda (haypiles per plot)", ylab = "Likelihood")
abline(v = 3, lty = 3)
plot(lambda_grid, nll_vals, type = "l", lwd = 2,
     xlab = "lambda (haypiles per plot)", ylab = "Negative log-likelihood")
abline(v = 3, lty = 3)
par(mfrow = c(1, 1))

dpois(haypiles, 1.5)   # each plot's probability at lambda = 1.5 (as in the widget)
#  [PAUSE] which count is most to blame for the poor likelihood at lambda = 1.5? at 6?

# 4. The maximum likelihood estimate -------------------------------------------

prod(dpois(haypiles, 3))                   # likelihood at the MLE: tiny, and that's fine

lambda_grid[which.min(nll_vals)]           # grid search

fit_pika <- optim(par = 2, fn = nll_pois, y = haypiles,
                  method = "Brent", lower = 0.01, upper = 20)
fit_pika$par                               # optimizer
fit_pika$value                             # NLL at the MLE

coarse <- seq(1, 6, by = 0.5)
par(mar = c(4, 4.5, 0.5, 0.5))
plot(lambda_grid, nll_vals, type = "l", lwd = 2, col = "grey40",
     xlab = "lambda (haypiles per plot)", ylab = "Negative log-likelihood")
points(coarse, sapply(coarse, nll_pois, y = haypiles), pch = 16)
points(fit_pika$par, fit_pika$value, pch = 4, cex = 2, lwd = 3, col = "red")
legend("top", c("grid points", "optim() result"), pch = c(16, 4),
       col = c("black", "red"), bty = "n")
#  [PAUSE] the likelihood at the MLE is 7.1e-6. does that mean the model fits badly?

# Bonus: fitting by simulation -------------------------------------------------
#  black cherry trees: volume ~ Normal(a * D^b, sigma)

data(trees)
diam <- trees$Girth      # diameter in inches (mislabeled "Girth")
vol <- trees$Volume      # cubic feet

power_curve <- function(d, a, b) a * d^b

# 95% range of volumes simulated at each diameter, for parameters (a, b, sigma)
sim_band <- function(par, d_seq, n_sims = 1000) {
  sims <- replicate(n_sims, rnorm(length(d_seq), power_curve(d_seq, par[1], par[2]), par[3]))
  apply(sims, 1, quantile, probs = c(0.025, 0.975))
}

set.seed(6)
d_seq <- seq(8, 21, length.out = 50)
guesses <- list(c(0.1, 2.2, 3), c(0.5, 1.6, 3), c(0.5, 1.6, 15))

par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.5))
for (g in guesses) {
  band <- sim_band(g, d_seq)
  plot(diam, vol, pch = 16, ylim = c(-20, 110), xlab = "Diameter (in)",
       ylab = "Volume (cubic ft)",
       main = paste0("a = ", g[1], ", b = ", g[2], ", sigma = ", g[3]), cex.main = 1)
  polygon(c(d_seq, rev(d_seq)), c(band[1, ], rev(band[2, ])),
          col = adjustcolor("steelblue", 0.25), border = NA)
  lines(d_seq, power_curve(d_seq, g[1], g[2]), lwd = 2, col = "steelblue4")
}
par(mfrow = c(1, 1))
#  [PAUSE] why not just increase sigma until every point falls inside the band?

nll_power <- function(par, d, v) {
  -sum(dnorm(v, mean = power_curve(d, par[1], par[2]), sd = par[3], log = TRUE))
}

sapply(guesses, nll_power, d = diam, v = vol)

# FITTING A MODEL WITH SEVERAL PARAMETERS ======================================

# 1. Likelihood for continuous data --------------------------------------------

par(mar = c(4, 4.5, 0.5, 0.5))
plot(diam, vol, pch = 16, xlab = "Diameter (in)", ylab = "Volume (cubic ft)")

par_try <- c(a = 0.087, b = 2.24, sigma = 3.2)
mu_1 <- power_curve(diam[1], par_try["a"], par_try["b"])   # expected volume
dens_1 <- dnorm(vol[1], mean = mu_1, sd = par_try["sigma"])
c(expected = unname(mu_1), density = unname(dens_1))

par(mar = c(4, 4.5, 0.5, 0.5))
curve(dnorm(x, mu_1, par_try["sigma"]), from = 0, to = 22, lwd = 2,
      xlab = "Volume (cubic ft) of a tree with diameter 8.3 in", ylab = "Probability density")
abline(v = vol[1], col = "red", lwd = 2)

m3 <- 0.0283168   # cubic meters per cubic foot
dnorm(vol[1] * m3, mu_1 * m3, par_try["sigma"] * m3)   # same tree, density per cubic meter
prod(dnorm(vol, power_curve(diam, 0.087, 2.24), 3.2))              # product of 31 densities
sum(dnorm(vol, power_curve(diam, 0.087, 2.24), 3.2, log = TRUE))   # log-likelihood

# 2. Writing a negative log-likelihood function --------------------------------
#  nll_power() above: parameters in, expected values, log densities, minus the sum

nll_power(par_try, d = diam, v = vol)

# 3. Starting values and optim() -----------------------------------------------
#  b ~ 2 (cylinder); a ~ 10.3 / 8.3^2 from the first tree; sigma ~ 3 by eye

start <- c(a = 0.15, b = 2, sigma = 3)
par(mar = c(4, 4.5, 0.5, 0.5))
plot(diam, vol, pch = 16, xlab = "Diameter (in)", ylab = "Volume (cubic ft)")
curve(power_curve(x, start["a"], start["b"]), add = TRUE, lwd = 2, lty = 2)

fit_trees <- optim(par = start, fn = nll_power, d = diam, v = vol)
fit_trees$par           # MLEs of a, b, and sigma
fit_trees$value         # NLL at the MLE
fit_trees$convergence   # 0 means optim() thinks it converged

round(optim(c(a = 0.05, b = 2.5, sigma = 5), nll_power, d = diam, v = vol)$par, 3)

# 4. Checking the fit ----------------------------------------------------------

set.seed(7)
band_fit <- sim_band(fit_trees$par, d_seq)
par(mar = c(4, 4.5, 0.5, 0.5))
plot(diam, vol, pch = 16, xlab = "Diameter (in)", ylab = "Volume (cubic ft)",
     ylim = c(0, 85))
polygon(c(d_seq, rev(d_seq)), c(band_fit[1, ], rev(band_fit[2, ])),
        col = adjustcolor("steelblue", 0.25), border = NA)
curve(power_curve(x, fit_trees$par["a"], fit_trees$par["b"]), add = TRUE,
      lwd = 2, col = "steelblue4")
#  [PAUSE] does the same sigma for small and large trees look right?

nll_exp <- function(par, d, v) -sum(dnorm(v, par[1] * exp(par[2] * d), par[3], log = TRUE))
fit_exp <- optim(c(5, 0.15, 3), nll_exp, d = diam, v = vol)
c(exponential = fit_exp$value, power_law = fit_trees$value)
#  [PAUSE] is a difference that size enough to prefer the power law?

# Bonus: deviance and goodness of fit ------------------------------------------
#  deviance = 2 * (NLL_model - NLL_saturated) ~ chi-square, df = 6 - 1

nll_saturated <- -sum(dpois(haypiles, lambda = haypiles, log = TRUE))
deviance_pika <- 2 * (nll_pois(3, haypiles) - nll_saturated)
c(deviance = deviance_pika, df = 5,
  p_value = pchisq(deviance_pika, df = 5, lower.tail = FALSE))

# A POWER ANALYSIS IN CODE =====================================================

# 2. WOD: find the bugs (each function has one bug; don't run these) -----------

# detect_prob <- function(k, pobs = 0.02) {
#   1 - pobs^k
# }
#
# sim_detect <- function(R = 1000, visits = 3, p = 0.4) {
#   detected <- logical(R)
#   for (i in seq_len(R)) {
#     n_det <- rbinom(1, size = visits, prob = p)
#     detected <- n_det > 0
#   }
#   mean(detected)
# }
#  [PAUSE] find each bug. what one-line check would have caught it?

# 4. Building and testing the pieces -------------------------------------------
#  N_t = 1000 * lambda^t, lambda = 0.75^(1/25); detection per survey = 1 - 0.98^k

detect_prob <- function(k, pobs = 0.02) {
  1 - (1 - pobs)^k
}
detect_prob(1)   # one person-day should equal pobs
detect_prob(3)

true_abundance <- function(t, N0 = 1000, lambda = 0.75^(1/25)) {
  N0 * lambda^t
}
true_abundance(c(0, 25))   # should return 1000 and 750 together

simulate_counts <- function(t_survey, k, N0 = 1000, lambda = 0.75^(1/25),
                            pobs = 0.02) {
  N_t <- true_abundance(t_survey, N0, lambda)
  rbinom(length(t_survey), size = round(N_t), prob = detect_prob(k, pobs))
}
set.seed(11)
counts <- simulate_counts(1:25, k = 3)
length(counts)   # should be 25

par(mar = c(4, 4.5, 0.5, 0.5))
plot(1:25, counts, pch = 16, xlab = "Year", ylab = "Count",
     ylim = c(0, max(counts) * 1.1))

detect_decline <- function(t_survey, counts, alpha = 0.05) {
  s <- summary(lm(log(counts) ~ t_survey))$coefficients[2, ]
  unname(s[1] < 0 && s[4] <= alpha)
}
detect_decline(1:25, counts)

set.seed(12)
null_hits <- replicate(1000, {
  detect_decline(1:25, simulate_counts(1:25, 3, lambda = 1))
})
mean(null_hits)   # stable population: should be about 0.025

# 5. Power across designs ------------------------------------------------------

estimate_power <- function(k, s, R = 1000, lambda = 0.75^(1/25),
                           alpha = 0.05) {
  t_survey <- seq(1, 25, by = s)
  mean(replicate(R, {
    counts <- simulate_counts(t_survey, k, lambda = lambda)
    detect_decline(t_survey, counts, alpha)
  }))
}
set.seed(14)
estimate_power(3, 1, lambda = 1)   # false-decline rate
estimate_power(3, 1)               # power

set.seed(746)
grid <- expand.grid(k = 1:15, s = 1:5)
grid$power <- mapply(estimate_power, grid$k, grid$s,
                     MoreArgs = list(R = 1000))   # about 20-30 s

par(mar = c(4, 4.5, 0.5, 0.5))
plot(NA, xlim = c(1, 15), ylim = c(0, 1), xlab = "Person-days per survey, k",
     ylab = "Estimated power")
abline(h = 0.75, lty = 2, col = "grey40")
for (ss in 1:5) {
  lines(grid$k[grid$s == ss], grid$power[grid$s == ss], col = ss, lwd = 2)
}
legend("bottomright", legend = paste("s =", 1:5), col = 1:5, lwd = 2,
       bty = "n", cex = 0.8)

cost_of <- function(k, s) length(seq(1, 25, by = s)) * (2000 + 200 * k)
grid$cost <- mapply(cost_of, grid$k, grid$s)
ok <- grid[grid$power >= 0.75, ]
head(ok[order(ok$cost), ], 5)

sqrt(0.75 * 0.25 / 1000)   # Monte Carlo SE of a power near 0.75

top3 <- head(ok[order(ok$cost), ], 3)
set.seed(999)
top3$power_5000 <- mapply(estimate_power, top3$k, top3$s, MoreArgs = list(R = 5000))
top3$mc_se <- sqrt(top3$power_5000 * (1 - top3$power_5000) / 5000)
top3
#  [PAUSE] k = 10, s = 4 is cheapest at about 0.76. recommend it, or k = 11, s = 4?
