#  NRES 746, Week 13: Bayesian hierarchical models; splines, GAMs, and Gaussian processes =====
#  Sections match the Week 13 web page. Needs mvtnorm, lme4, cmdstanr, posterior, brms, splines, mgcv.

# BAYESIAN HIERARCHICAL MODELS ===========================================================

rm(list = ls())
library(mvtnorm)
library(lme4)

# 2. WOD: partial pooling by hand ---------------------------------------------------------

# a_j = mu + w_j (ybar_j - mu), w_j = sb^2 / (sb^2 + s^2 / n_j)
w_wod <- 0.5^2 / (0.5^2 + 1^2 / c(A = 4, B = 16))
3 + w_wod * (c(A = 4.2, B = 4.0) - 3)

# 3. Partial pooling ---------------------------------------------------------------------

# 12 pinyon stands, records of 5-30 years; simulated
# truth: mean intercept 0.9, mean slope 0.25, tau = (0.15, 0.10), rho = -0.6, sigma = 0.12
set.seed(12)
n_years <- rep(c(30, 24, 18, 12, 8, 5), each = 2)        # record length per stand
tau_true <- c(0.15, 0.10)                                 # among-stand SDs: intercept, slope
rho_true <- -0.6
G_true <- diag(tau_true) %*% matrix(c(1, rho_true, rho_true, 1), 2) %*% diag(tau_true)
ab_true <- rmvnorm(12, c(0.9, 0.25), G_true)              # each stand's intercept and slope
stand_names <- sprintf("S%02d", 1:12)
drought <- data.frame(stand = factor(rep(stand_names, n_years)),
                      year = unlist(lapply(n_years, function(n) seq(2021 - n, 2020))))
j <- as.integer(drought$stand)
drought$precip <- round(rgamma(nrow(drought), shape = 9, scale = 22 / 9), 1)
drought$precip10 <- (drought$precip - 20) / 10           # 10-cm units, 0 = 20 cm
drought$width <- round(ab_true[j, 1] + ab_true[j, 2] * drought$precip10 +
                         rnorm(nrow(drought), 0, 0.12), 2)

# no pooling: one least-squares line per stand; maximum likelihood random slope fit
nopool <- t(sapply(split(drought, drought$stand), function(s) coef(lm(width ~ precip10, data = s))))
fit_rs_lmer <- lmer(width ~ precip10 + (precip10 | stand), data = drought, REML = FALSE)
ml <- c(fixef(fit_rs_lmer), as.data.frame(VarCorr(fit_rs_lmer))$sdcor)
names(ml) <- c("mu_a", "mu_b", "tau_a", "tau_b", "rho", "sigma")
round(ml, 3)


# random intercept (no precipitation): a_j = mu + w_j (ybar_j - mu), w_j = sb^2 / (sb^2 + s^2 / n_j)
fit_ri0 <- lmer(width ~ 1 + (1 | stand), data = drought, REML = FALSE)
sd_ri0 <- as.data.frame(VarCorr(fit_ri0))$sdcor          # sigma_b, sigma
ybar <- tapply(drought$width, drought$stand, mean)
n_j <- as.vector(table(drought$stand))
w <- sd_ri0[1]^2 / (sd_ri0[1]^2 + sd_ri0[2]^2 / n_j)
round(cbind(own_mean = ybar, n_j, w, formula = fixef(fit_ri0) + w * (ybar - fixef(fit_ri0)),
            lmer = coef(fit_ri0)$stand[, 1]), 3)

#  [PAUSE] S10 and S12 have similar own means but 8 vs. 5 years. Which moves more toward the overall mean?

# 5. Priors for a covariance matrix ------------------------------------------------------

# LKJ(eta) marginal for one correlation in a K x K matrix: (rho + 1)/2 ~ Beta(eta - 1 + K/2, same)
lkj_marginal <- function(rho, eta, K) dbeta((rho + 1) / 2, eta - 1 + K / 2, eta - 1 + K / 2) / 2
curve(lkj_marginal(x, eta = 2, K = 2), -1, 1, lwd = 2, xlab = "rho", ylab = "Prior density",
      ylim = c(0, 2))
curve(lkj_marginal(x, eta = 1, K = 2), add = TRUE, lty = 2)
curve(lkj_marginal(x, eta = 1, K = 10), add = TRUE, lty = 3, lwd = 2)

#  [PAUSE] For a 10 x 10 matrix with eta = 1, how likely is any one correlation to exceed 0.5 in absolute value?

# 6. The random slope model in brms ---------------------------------------------------------

library(brms)
library(posterior)
fit_brms <- brm(width ~ precip10 + (precip10 | stand), data = drought,
                prior = c(prior(normal(1, 0.5), class = Intercept),
                          prior(normal(0, 0.5), class = b),
                          prior(normal(0, 0.3), class = sd),
                          prior(lkj(2), class = cor),
                          prior(normal(0, 0.3), class = sigma)),
                backend = "cmdstanr", chains = 4, cores = 4, seed = 746, refresh = 0)

# check the sampler, then the estimates
summary(fit_brms)
prior_summary(fit_brms)

# posterior for rho against its LKJ(2) prior and the MLE
draws <- as_draws_df(fit_brms)
hist(draws$cor_stand__Intercept__precip10, breaks = 40, freq = FALSE, xlim = c(-1, 1),
     main = "", xlab = "rho")
curve(lkj_marginal(x, eta = 2, K = 2), add = TRUE, lwd = 2, lty = 2)
abline(v = ml["rho"], col = "firebrick", lwd = 2)

# shrinkage: each stand's own fit -> its posterior mean
ab_post <- coef(fit_brms)$stand[, "Estimate", ]
plot(nopool, xlab = "Intercept: width at 20 cm (mm)", ylab = "Slope (mm per 10 cm)")
arrows(nopool[, 1], nopool[, 2], ab_post[, 1], ab_post[, 2], length = 0.07, col = "steelblue")
points(ab_post, pch = 19, col = "steelblue")

#  [PAUSE] Why is the posterior for rho so much wider than the others, and its median closer to 0 than the MLE?

# try it: change lkj(2) to lkj(1), then lkj(10), and refit

# Bonus: the same model written in Stan ------------------------------------------------

library(cmdstanr)
slopes_code <- "
data {
  int<lower=1> N;                       // stand-years
  int<lower=1> J;                       // stands
  array[N] int<lower=1, upper=J> stand;
  vector[N] x;                          // precipitation (10 cm, centered at 20 cm)
  vector[N] y;                          // ring width (mm)
}
parameters {
  vector[2] mu;                         // mean intercept and slope
  vector<lower=0>[2] tau;               // among-stand SDs of intercepts and slopes
  cholesky_factor_corr[2] L;            // Cholesky factor of the correlation matrix
  matrix[2, J] z;                       // standardized stand effects
  real<lower=0> sigma;                  // residual SD
}
transformed parameters {
  // row j: stand j's intercept and slope
  matrix[J, 2] ab = rep_matrix(mu', J) + (diag_pre_multiply(tau, L) * z)';
}
model {
  mu[1] ~ normal(1, 0.5);
  mu[2] ~ normal(0, 0.5);
  tau ~ normal(0, 0.3);                 // half-normal, because tau > 0
  L ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();
  sigma ~ normal(0, 0.3);
  y ~ normal(ab[stand, 1] + ab[stand, 2] .* x, sigma);
}
generated quantities {
  real rho = multiply_lower_tri_self_transpose(L)[1, 2];
}
"
slopes_model <- cmdstan_model(write_stan_file(slopes_code))   # compiles (slow, once)
stan_data <- list(N = nrow(drought), J = 12, stand = as.integer(drought$stand),
                  x = drought$precip10, y = drought$width)
fit_stan <- slopes_model$sample(data = stan_data, chains = 4, parallel_chains = 4,
                                seed = 746, refresh = 0)

# check the sampler, then the estimates
fit_stan$diagnostic_summary()
fit_stan$summary(c("mu", "tau", "rho", "sigma"))

# try it: change lkj_corr_cholesky(2) to (1), then (10), recompile, and refit

stancode(fit_brms)          # the program brms wrote: same Cholesky factor, non-centered effects

# SPLINES, GAMS, AND GAUSSIAN PROCESSES ==================================================

library(splines)
library(mgcv)

# 1. A thermal performance curve ----------------------------------------------------------

# Bd growth rate (per day) for 48 cultures on a thermal gradient plate; simulated
# truth: a two-piece curve (Gaussian rise, quadratic fall), peak 0.3 at 20 C, zero at 28 C
tpc_true <- function(temp) 0.3 * ifelse(temp <= 20, exp(-((temp - 20) / 8)^2), 1 - ((temp - 20) / 8)^2)
set.seed(13)
bd <- data.frame(temp = round(seq(4, 29, length.out = 48) + runif(48, -0.2, 0.2), 1))
bd$r <- round(tpc_true(bd$temp) + rnorm(48, 0, 0.03), 3)
head(bd)

# Gaussian (bell) curve with Normal noise, by maximum likelihood
bell <- function(temp, rmax, topt, width) rmax * exp(-(temp - topt)^2 / (2 * width^2))
nll_bell <- function(par) {
  -sum(dnorm(bd$r, bell(bd$temp, par[1], par[2], par[3]), exp(par[4]), log = TRUE))
}
fit_bell <- optim(c(0.3, 20, 5, log(0.05)), nll_bell)
round(setNames(c(fit_bell$par[1:3], exp(fit_bell$par[4])), c("rmax", "topt", "width", "sigma")), 3)

plot(r ~ temp, data = bd, pch = 16, xlab = "Temperature (C)", ylab = "Growth rate r (per day)")
curve(bell(x, fit_bell$par[1], fit_bell$par[2], fit_bell$par[3]), add = TRUE, lwd = 2, col = "firebrick")

#  [PAUSE] Which features of the scatter does the bell curve get wrong? Which matter for management?

# 2. Splines are linear models ------------------------------------------------------------

# f(x) = sum_k beta_k B_k(x): each basis function is a column of the model matrix
X <- bs(bd$temp, df = 6)          # 6 cubic B-spline columns (3 interior knots)
round(head(X), 3)
attr(X, "knots")
fit_bs6 <- lm(r ~ bs(temp, df = 6), data = bd)
fit_bs25 <- lm(r ~ bs(temp, df = 25), data = bd)

temp_grid <- seq(4, 29, by = 0.05)
matplot(temp_grid, predict(X, temp_grid), type = "l", lty = 1, lwd = 2,
        xlab = "Temperature (C)", ylab = "Basis function value")
abline(v = attr(X, "knots"), lty = 3)
plot(r ~ temp, data = bd, pch = 16, xlab = "Temperature (C)", ylab = "Growth rate r (per day)")
lines(temp_grid, predict(fit_bs6, data.frame(temp = temp_grid)), lwd = 2, col = "steelblue")
lines(temp_grid, predict(fit_bs25, data.frame(temp = temp_grid)), lty = 2)

#  [PAUSE] Random-walk weights with a small step SD: what would you assume about the weights to get a smooth curve?

# 3. Penalizing wiggliness ----------------------------------------------------------------

# penalized NLL = NLL + lambda * sum of squared second differences of the weights
# a Normal(0, tau^2) prior on each difference gives the same penalty, lambda = 1 / (2 tau^2)

#  [PAUSE] As lambda grows, why does the curve become a sloped line rather than a flat one?

# 4. Fitting a GAM -----------------------------------------------------------------------

fit_gam <- gam(r ~ s(temp, k = 20), data = bd, method = "REML")
summary(fit_gam)
gam.vcomp(fit_gam)          # the smooth's SD, printed like a variance component

# AIC: bell vs. unpenalized spline vs. GAM (the GAM's parameter count is its edf)
data.frame(model = c("bell", "bs, 6 basis functions", "gam"),
           AIC = c(2 * fit_bell$value + 8, AIC(fit_bs6), AIC(fit_gam)))

# GAM curve, its optimum, and where it crosses zero
pg <- predict(fit_gam, data.frame(temp = temp_grid), se.fit = TRUE)
plot(r ~ temp, data = bd, pch = 16, xlab = "Temperature (C)", ylab = "Growth rate r (per day)")
polygon(c(temp_grid, rev(temp_grid)), c(pg$fit + 1.96 * pg$se.fit, rev(pg$fit - 1.96 * pg$se.fit)),
        col = adjustcolor("steelblue", 0.2), border = NA)
lines(temp_grid, pg$fit, lwd = 2, col = "steelblue")
curve(bell(x, fit_bell$par[1], fit_bell$par[2], fit_bell$par[3]), add = TRUE, lwd = 2, lty = 2, col = "firebrick")
topt_gam <- temp_grid[which.max(pg$fit)]
ctmax_gam <- temp_grid[temp_grid > topt_gam & pg$fit < 0][1]
c(topt = topt_gam, ctmax = ctmax_gam)

# 5. Gaussian processes -------------------------------------------------------------------

# y ~ MVN(m 1, K + sigma^2 I), K_ij = alpha^2 exp(-(x_i - x_j)^2 / (2 ell^2))
dist_temp <- abs(outer(bd$temp, bd$temp, "-"))
nll_gp <- function(par) {
  m <- par[1]
  alpha <- exp(par[2])
  ell <- exp(par[3])
  sigma <- exp(par[4])
  V <- alpha^2 * exp(-dist_temp^2 / (2 * ell^2)) + sigma^2 * diag(nrow(bd))
  -dmvnorm(bd$r, rep(m, nrow(bd)), V, log = TRUE)
}
fit_gp <- optim(c(0.15, log(0.15), log(4), log(0.03)), nll_gp, method = "BFGS")
gp_par <- c(m = fit_gp$par[1], exp(fit_gp$par[2:4]))
names(gp_par) <- c("m", "alpha", "ell", "sigma")
round(gp_par, 3)

# prediction: f(x*) = m + k*' (K + sigma^2 I)^-1 (y - m)
gp_predict <- function(new_temp, p = gp_par) {
  K <- p["alpha"]^2 * exp(-dist_temp^2 / (2 * p["ell"]^2)) + p["sigma"]^2 * diag(nrow(bd))
  k_star <- p["alpha"]^2 * exp(-outer(new_temp, bd$temp, "-")^2 / (2 * p["ell"]^2))
  as.vector(p["m"] + k_star %*% solve(K, bd$r - p["m"]))
}
gp_fit_grid <- gp_predict(temp_grid)
max(abs(gp_fit_grid - pg$fit))          # largest gap between the GP and GAM curves
lines(temp_grid, gp_fit_grid, lty = 3, lwd = 2)

#  [PAUSE] Exponential vs. squared exponential kernel: what changes about the curves between observations?

# Discussion: predictions outside the data, at 32 C
c(bell = bell(32, fit_bell$par[1], fit_bell$par[2], fit_bell$par[3]),
  gam = as.vector(predict(fit_gam, data.frame(temp = 32))),
  gp = gp_predict(32))

#  [PAUSE] Why do the three predictions at 32 C disagree? Which would you use?

# Bonus: uncertainty in derived quantities ------------------------------------------------

Xp <- predict(fit_gam, data.frame(temp = temp_grid), type = "lpmatrix")   # basis at each grid point
set.seed(746)
beta_sims <- rmvnorm(2000, coef(fit_gam), vcov(fit_gam))
curves <- Xp %*% t(beta_sims)                                              # one column per curve
topt_sims <- temp_grid[apply(curves, 2, which.max)]
ctmax_sims <- apply(curves, 2, function(f) {
  above <- temp_grid > 20 & f < 0
  if (any(above)) temp_grid[which(above)[1]] else NA
})
round(rbind(topt = quantile(topt_sims, c(0.5, 0.025, 0.975)),
            ctmax = quantile(ctmax_sims, c(0.5, 0.025, 0.975), na.rm = TRUE)), 1)

# Bonus: a smooth in brms -----------------------------------------------------------------

library(brms)
fit_brms_gam <- brm(r ~ s(temp, k = 20), data = bd, backend = "cmdstanr",
                    chains = 4, cores = 4, seed = 746, refresh = 0,
                    control = list(adapt_delta = 0.95))
summary(fit_brms_gam)          # sds(stemp_1): SD of the smooth's wiggly part
stancode(fit_brms_gam)
