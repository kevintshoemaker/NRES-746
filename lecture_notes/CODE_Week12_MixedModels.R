#  NRES 746, Week 12: autocorrelation, mixed models in R, Bayesian hierarchical models =====
#  Sections match the Week 12 web page. Needs mvtnorm, lme4, nlme, glmmTMB, cmdstanr, posterior, brms.

# AUTOCORRELATION AND MIXED MODELS IN R ===================================================

rm(list = ls())
library(mvtnorm)

# 1. Big picture: the pinyon ring-width data ---------------------------------------------

# pinyon ring widths (mm): 4 stands x 5 trees x 25 years; simulated
set.seed(21)
years <- 1996:2020
n_yr <- length(years)
precip <- round(rgamma(n_yr, shape = 9, scale = 22 / 9), 1)   # winter precipitation (cm)
stand_eff <- rnorm(4, 0, 0.12)
ar_dev <- matrix(0, 4, n_yr)            # shared year-to-year deviations within each stand
for (j in 1:4) {
  ar_dev[j, 1] <- rnorm(1, 0, 0.08 / sqrt(1 - 0.6^2))
  for (t in 2:n_yr) ar_dev[j, t] <- 0.6 * ar_dev[j, t - 1] + rnorm(1, 0, 0.08)
}
rings <- expand.grid(year = years, tree = 1:5, stand = LETTERS[1:4])
rings$tree_id <- paste0(rings$stand, rings$tree)
tree_eff <- setNames(rnorm(20, 0, 0.1), unique(rings$tree_id))
j <- match(rings$stand, LETTERS)
t <- match(rings$year, years)
rings$precip <- precip[t]
rings$width <- round(0.35 + stand_eff[j] + tree_eff[rings$tree_id] + 0.025 * rings$precip +
                       ar_dev[cbind(j, t)] + rnorm(nrow(rings), 0, 0.1), 2)

# stand chronologies (mean of 5 trees per stand-year) and tree means (mean over 25 years)
chron <- aggregate(width ~ stand + year + precip, data = rings, FUN = mean)
chron <- chron[order(chron$stand, chron$year), ]
chron$stand <- factor(chron$stand)
tree_means <- aggregate(width ~ stand + tree_id, data = rings, FUN = mean)

# 2. Autocorrelation: AR(1) --------------------------------------------------------------

chron_a <- chron[chron$stand == "A", ]
chron_a$year_c <- chron_a$year - 2008          # centered year: 0 = 2008
fit_lm_a <- lm(width ~ precip + year_c, data = chron_a)
summary(fit_lm_a)$coefficients
par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 1))
plot(chron_a$year, resid(fit_lm_a), type = "b", pch = 19, xlab = "Year", ylab = "Residual")
abline(h = 0, lty = 2)
acf(resid(fit_lm_a), main = "", lag.max = 10)
par(mfrow = c(1, 1))

# lag-1 residual correlation in each stand (stand A's is inside the bounds, +/-0.39)
lag1 <- function(dat) {
  dat$year_c <- dat$year - 2008
  acf(resid(lm(width ~ precip + year_c, data = dat)), plot = FALSE)$acf[2]
}
round(sapply(split(chron, chron$stand), lag1), 2)

# AR(1) covariance: V_st = sigma^2 / (1 - rho^2) * rho^|s - t|
lags <- abs(outer(1:25, 1:25, "-"))
V_ar1 <- function(rho, sigma) sigma^2 / (1 - rho^2) * rho^lags
round(V_ar1(0.6, 0.08)[1:5, 1:5], 4)

# mean: a + b * precip + c * year_c; AR(1) errors with rho = tanh(par[4]), sigma = exp(par[5])
nll_ar1 <- function(par) {
  mu <- par[1] + par[2] * chron_a$precip + par[3] * chron_a$year_c
  V <- V_ar1(tanh(par[4]), exp(par[5]))
  -dmvnorm(chron_a$width, mu, V, log = TRUE)
}
fit_ar1 <- optim(c(coef(fit_lm_a), 0, log(0.1)), nll_ar1, method = "BFGS", hessian = TRUE)
c(a = fit_ar1$par[1], b = fit_ar1$par[2], c = fit_ar1$par[3],
  rho = tanh(fit_ar1$par[4]), sigma = exp(fit_ar1$par[5]))

library(nlme)
fit_gls <- gls(width ~ precip + year_c, data = chron_a,
               correlation = corAR1(form = ~ year), method = "ML")
coef(fit_gls)
coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)   # rho (nlme calls it Phi)

# trend: estimate, SE, p under independence vs. AR(1) (Hessian SE)
se_ar1 <- sqrt(diag(solve(fit_ar1$hessian)))
summary(fit_lm_a)$coefficients["year_c", ]
c(estimate = fit_ar1$par[3], se = se_ar1[3], p = 2 * pnorm(-abs(fit_ar1$par[3] / se_ar1[3])))

#  [PAUSE] With rho near 0.9, what does one AR(1) series look like? How would you tell it from a slow trend?

# 3. Why these models are harder to fit --------------------------------------------------

# random intercept for the 20 tree means: V = sigma_b^2 Z Z^T + sigma^2 I
y_tree <- tree_means$width
Z <- model.matrix(~ 0 + stand, data = tree_means)
V_ri <- function(sigma_b, sigma) sigma_b^2 * Z %*% t(Z) + sigma^2 * diag(20)
nll_ri <- function(par) -dmvnorm(y_tree, rep(par[1], 20), V_ri(exp(par[2]), exp(par[3])), log = TRUE)
fit_ri <- optim(c(mean(y_tree), log(0.1), log(0.1)), nll_ri, method = "BFGS")
c(mu = fit_ri$par[1], sigma_b = exp(fit_ri$par[2]), sigma = exp(fit_ri$par[3]))

# milliseconds per NLL evaluation, averaged over 'reps' evaluations
time_nll <- function(n, reps = 10) {
  y <- rnorm(n)
  V <- 0.5^abs(outer(1:n, 1:n, "-"))
  ms <- function(expr) 1000 * system.time(for (i in 1:reps) expr)[["elapsed"]] / reps
  c(n = n, independent = ms(sum(dnorm(y, 0, 1, log = TRUE))),
    mvn = ms(dmvnorm(y, rep(0, n), V, log = TRUE)))
}
set.seed(1)
timings <- t(sapply(c(250, 500, 1000, 2000), time_nll))
timings

# 4. Random slopes -----------------------------------------------------------------------

# drought sensitivity: 12 stands, records of 5-30 years; simulated
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
head(drought)

# no pooling: one least-squares line per stand
nopool <- t(sapply(split(drought, drought$stand), function(s) coef(lm(width ~ precip10, data = s))))
plot(nopool, pch = 19, xlab = "Intercept: width at 20 cm (mm)", ylab = "Slope (mm per 10 cm)")
text(nopool, labels = stand_names, pos = 4, cex = 0.7)

# 5. Mixed models in R -------------------------------------------------------------------

library(lme4)
fit_ri_lmer <- lmer(width ~ precip10 + (1 | stand), data = drought, REML = FALSE)
fit_rs_lmer <- lmer(width ~ precip10 + (precip10 | stand), data = drought, REML = FALSE)
fixef(fit_rs_lmer)
VarCorr(fit_rs_lmer)
AIC(fit_ri_lmer, fit_rs_lmer)
vc_rs <- as.data.frame(VarCorr(fit_rs_lmer))$sdcor        # tau_a, tau_b, rho, sigma

fit_rs_lme <- lme(width ~ precip10, random = ~ precip10 | stand, data = drought, method = "ML")

library(glmmTMB)
fit_rs_tmb <- glmmTMB(width ~ precip10 + (precip10 | stand), data = drought)   # ML by default

# by hand: one MVN block per stand, V_j = Z_j G Z_j^T + sigma^2 I
rows_by_stand <- split(seq_len(nrow(drought)), drought$stand)
nll_rs <- function(par) {
  tau <- exp(par[3:4])
  rho <- tanh(par[5])
  sigma <- exp(par[6])
  G <- matrix(c(tau[1]^2, rho * tau[1] * tau[2], rho * tau[1] * tau[2], tau[2]^2), 2)
  nll <- 0
  for (rows in rows_by_stand) {
    Zj <- cbind(1, drought$precip10[rows])
    Vj <- Zj %*% G %*% t(Zj) + sigma^2 * diag(length(rows))
    nll <- nll - dmvnorm(drought$width[rows], as.vector(Zj %*% par[1:2]), Vj, log = TRUE)
  }
  nll
}
fit_rs_optim <- optim(c(0.9, 0.25, log(0.1), log(0.1), 0, log(0.1)), nll_rs, method = "BFGS")

tmb_sd <- attr(VarCorr(fit_rs_tmb)$cond$stand, "stddev")
lme_vc <- VarCorr(fit_rs_lme)
compare <- rbind(
  optim = c(fit_rs_optim$par[1:2], exp(fit_rs_optim$par[3:4]), tanh(fit_rs_optim$par[5]),
            exp(fit_rs_optim$par[6]), -fit_rs_optim$value),
  lmer = c(fixef(fit_rs_lmer), vc_rs, logLik(fit_rs_lmer)),
  lme = c(fixef(fit_rs_lme), as.numeric(lme_vc[1:2, "StdDev"]), as.numeric(lme_vc[2, "Corr"]),
          fit_rs_lme$sigma, logLik(fit_rs_lme)),
  glmmTMB = c(fixef(fit_rs_tmb)$cond, tmb_sd,
              attr(VarCorr(fit_rs_tmb)$cond$stand, "correlation")[1, 2],
              sigma(fit_rs_tmb), logLik(fit_rs_tmb)))
colnames(compare) <- c("mu_a", "mu_b", "tau_a", "tau_b", "rho", "sigma", "logLik")
round(compare, 4)

#  [PAUSE] The intercept-slope correlation is about -0.57. What does that say about these stands, ecologically?

# Discussion: the four short-record stands alone
short <- droplevels(drought[drought$stand %in% c("S09", "S10", "S11", "S12"), ])
fit_short <- lmer(width ~ precip10 + (precip10 | stand), data = short, REML = FALSE)
VarCorr(fit_short)

# Bonus: a random intercept and AR(1) together -------------------------------------------

Z_chron <- model.matrix(~ 0 + stand, data = chron)
V_combined <- function(sigma_b, rho, sigma) {
  sigma_b^2 * Z_chron %*% t(Z_chron) + kronecker(diag(4), V_ar1(rho, sigma))
}
nll_combined <- function(par) {
  mu <- par[1] + par[2] * chron$precip
  V <- V_combined(exp(par[3]), tanh(par[4]), exp(par[5]))
  -dmvnorm(chron$width, mu, V, log = TRUE)
}
fit_comb <- optim(c(0.4, 0.02, log(0.1), 0, log(0.1)), nll_combined, method = "BFGS",
                  control = list(maxit = 1000))
c(a = fit_comb$par[1], b = fit_comb$par[2], sigma_b = exp(fit_comb$par[3]),
  rho = tanh(fit_comb$par[4]), sigma = exp(fit_comb$par[5]))

fit_lme <- lme(width ~ precip, random = ~ 1 | stand, data = chron,
               correlation = corAR1(form = ~ year | stand), method = "ML")
fixef(fit_lme)
coef(fit_lme$modelStruct$corStruct, unconstrained = FALSE)
c(-fit_comb$value, logLik(fit_lme))       # same maximized log-likelihood

# Bonus: a Poisson GLMM by hand ----------------------------------------------------------

# cones on 8 trees per stand; log(lambda) = b0 + b1 * vigor + b_j; simulated
set.seed(12)
cones <- data.frame(stand = factor(rep(stand_names, each = 8)),
                    vigor = round(rnorm(96), 2))
b_stand <- rnorm(12, 0, 0.5)
cones$count <- rpois(96, exp(1.5 + 0.4 * cones$vigor + b_stand[as.integer(cones$stand)]))

# each stand's likelihood integrates over its random intercept b
nll_glmm <- function(par) {
  sigma_b <- exp(par[3])
  nll <- 0
  for (g in split(cones, cones$stand)) {
    stand_lik <- function(b) {
      sapply(b, function(bi) exp(sum(dpois(g$count, exp(par[1] + par[2] * g$vigor + bi), log = TRUE)))) *
        dnorm(b, 0, sigma_b)
    }
    nll <- nll - log(integrate(stand_lik, -Inf, Inf)$value)
  }
  nll
}
fit_glmm <- optim(c(1, 0, log(0.3)), nll_glmm, method = "BFGS")
fit_glmm_tmb <- glmmTMB(count ~ vigor + (1 | stand), family = poisson, data = cones)
rbind(integrate = c(fit_glmm$par[1:2], exp(fit_glmm$par[3]), -fit_glmm$value),
      glmmTMB = c(fixef(fit_glmm_tmb)$cond, attr(VarCorr(fit_glmm_tmb)$cond$stand, "stddev"),
                  logLik(fit_glmm_tmb)))

# BAYESIAN HIERARCHICAL MODELS ===========================================================

# 2. Partial pooling ---------------------------------------------------------------------

# random intercept (no precipitation): a_j = mu + w_j (ybar_j - mu), w_j = sb^2 / (sb^2 + s^2 / n_j)
fit_ri0 <- lmer(width ~ 1 + (1 | stand), data = drought, REML = FALSE)
sd_ri0 <- as.data.frame(VarCorr(fit_ri0))$sdcor          # sigma_b, sigma
ybar <- tapply(drought$width, drought$stand, mean)
n_j <- as.vector(table(drought$stand))
w <- sd_ri0[1]^2 / (sd_ri0[1]^2 + sd_ri0[2]^2 / n_j)
round(cbind(own_mean = ybar, n_j, w, formula = fixef(fit_ri0) + w * (ybar - fixef(fit_ri0)),
            lmer = coef(fit_ri0)$stand[, 1]), 3)

#  [PAUSE] S10 and S12 have similar own means but 8 vs. 5 years. Which moves more toward the overall mean?

# 3. Priors for a covariance matrix ------------------------------------------------------

# LKJ(eta) marginal for one correlation in a K x K matrix: (rho + 1)/2 ~ Beta(eta - 1 + K/2, same)
lkj_marginal <- function(rho, eta, K) dbeta((rho + 1) / 2, eta - 1 + K / 2, eta - 1 + K / 2) / 2
curve(lkj_marginal(x, eta = 2, K = 2), -1, 1, lwd = 2, xlab = "rho", ylab = "Prior density",
      ylim = c(0, 2))
curve(lkj_marginal(x, eta = 1, K = 2), add = TRUE, lty = 2)
curve(lkj_marginal(x, eta = 1, K = 10), add = TRUE, lty = 3, lwd = 2)

#  [PAUSE] For a 10 x 10 matrix with eta = 1, how likely is any one correlation to exceed 0.5 in absolute value?

# 4. Random slopes in Stan ---------------------------------------------------------------

library(cmdstanr)
library(posterior)
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

# posterior for rho against its LKJ(2) prior and the MLE
draws <- as_draws_df(fit_stan$draws(c("rho", "ab")))
hist(draws$rho, breaks = 40, freq = FALSE, xlim = c(-1, 1), main = "", xlab = "rho")
curve(lkj_marginal(x, eta = 2, K = 2), add = TRUE, lwd = 2, lty = 2)
abline(v = compare["lmer", "rho"], col = "firebrick", lwd = 2)

# shrinkage: each stand's own fit -> its posterior mean
ab_post <- cbind(sapply(1:12, function(s) mean(draws[[sprintf("ab[%d,1]", s)]])),
                 sapply(1:12, function(s) mean(draws[[sprintf("ab[%d,2]", s)]])))
plot(nopool, xlab = "Intercept: width at 20 cm (mm)", ylab = "Slope (mm per 10 cm)")
arrows(nopool[, 1], nopool[, 2], ab_post[, 1], ab_post[, 2], length = 0.07, col = "steelblue")
points(ab_post, pch = 19, col = "steelblue")

#  [PAUSE] Why is the posterior for rho so much wider than the others, and its median closer to 0 than the MLE?

# try it: change lkj_corr_cholesky(2) to (1), then (10), recompile, and refit

# 5. The same model in brms --------------------------------------------------------------

library(brms)
fit_brms <- brm(width ~ precip10 + (precip10 | stand), data = drought,
                prior = c(prior(normal(0, 0.3), class = sd),
                          prior(lkj(2), class = cor),
                          prior(normal(0, 0.3), class = sigma)),
                backend = "cmdstanr", chains = 4, cores = 4, seed = 746, refresh = 0)
summary(fit_brms)
prior_summary(fit_brms)
stancode(fit_brms)          # the Stan program brms wrote
