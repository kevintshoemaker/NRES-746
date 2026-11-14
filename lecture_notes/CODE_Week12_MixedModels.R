#  NRES 746, Week 12: autocorrelation and mixed models in R =====
#  Sections match the Week 12 web page. Needs mvtnorm, lme4, nlme, glmmTMB.
#  Package syntax: lecture_notes/NOTES_Reference_MixedModelSyntax.pdf (nlme, lme4, glmmTMB side by side).

# AUTOCORRELATION AND MIXED MODELS IN R ===================================================

rm(list = ls())
library(mvtnorm)

# Pinyon ring-width data ----------------------------------------------------------------

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

# 1. Autocorrelation: AR(1) ---------------------------------------------------------------

chron_a <- chron[chron$stand == "A", ]
chron_a$year_c <- chron_a$year - 2008          # centered year: 0 = 2008
fit_lm_a <- lm(width ~ precip + year_c, data = chron_a)
summary(fit_lm_a)$coefficients
par(mfrow = c(3, 1), mar = c(4, 4.5, 2, 1))
plot(chron_a$year, chron_a$width, type = "b", pch = 19, xlab = "Year",
     ylab = "Mean ring width (mm)", main = "Stand A chronology")
trend <- predict(fit_lm_a, newdata = data.frame(precip = mean(chron_a$precip),
                                                year_c = chron_a$year_c))
lines(chron_a$year, trend, lwd = 2, col = "red")   # year trend at mean precipitation
plot(chron_a$year, resid(fit_lm_a), type = "b", pch = 19, xlab = "Year",
     ylab = "Residual", main = "Residuals from lm()")
abline(h = 0, lty = 2)
acf(resid(fit_lm_a), lag.max = 10, main = "")
title("Autocorrelation of the residuals", line = 0.8)
par(mfrow = c(1, 1))

# lag-1 residual correlation in each stand (stand A's is inside the bounds, +/-0.39)
lag1 <- function(dat) {
  dat$year_c <- dat$year - 2008
  acf(resid(lm(width ~ precip + year_c, data = dat)), plot = FALSE)$acf[2]
}
round(sapply(split(chron, chron$stand), lag1), 2)

# AR(1) covariance: V_st = sigma_e^2 / (1 - phi^2) * phi^|s - t|
lags <- abs(outer(1:25, 1:25, "-"))
V_ar1 <- function(phi, sigma_e) sigma_e^2 / (1 - phi^2) * phi^lags
round(V_ar1(0.6, 0.08)[1:5, 1:5], 4)

# mean: a + b * precip + c * year_c; AR(1) errors with phi = tanh(par[4]), sigma_e = exp(par[5])
nll_ar1 <- function(par) {
  mu <- par[1] + par[2] * chron_a$precip + par[3] * chron_a$year_c
  V <- V_ar1(tanh(par[4]), exp(par[5]))
  -dmvnorm(chron_a$width, mu, V, log = TRUE)
}
fit_ar1 <- optim(c(coef(fit_lm_a), 0, log(0.1)), nll_ar1, method = "BFGS", hessian = TRUE)
c(a = fit_ar1$par[1], b = fit_ar1$par[2], c = fit_ar1$par[3],
  phi = tanh(fit_ar1$par[4]), sigma_e = exp(fit_ar1$par[5]))

library(nlme)
fit_gls <- gls(width ~ precip + year_c, data = chron_a,
               correlation = corAR1(form = ~ year), method = "ML")
coef(fit_gls)
coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)   # phi (nlme calls it Phi)

# trend: estimate, SE, p under independence vs. AR(1) (Hessian SE)
se_ar1 <- sqrt(diag(solve(fit_ar1$hessian)))
summary(fit_lm_a)$coefficients["year_c", ]
c(estimate = fit_ar1$par[3], se = se_ar1[3], p = 2 * pnorm(-abs(fit_ar1$par[3] / se_ar1[3])))

#  [PAUSE] With phi near 0.9, what does one AR(1) series look like? How would you tell it from a slow trend?

# AR errors or a lagged response?
# AR(1) errors rewritten: y_t = a(1 - phi) + phi y_{t-1} + b x_t - phi b x_{t-1} + e_t
chron_a$width_lag <- c(NA, head(chron_a$width, -1))      # last year's width
chron_a$precip_lag <- c(NA, head(chron_a$precip, -1))    # last year's precipitation
fit_lag <- lm(width ~ width_lag + precip + year_c, data = chron_a)
fit_lag_x <- lm(width ~ width_lag + precip + precip_lag + year_c, data = chron_a)
round(coef(fit_lag), 4)          # gamma near 0
round(coef(fit_lag_x), 4)        # gamma near phi; precip_lag near -gamma * beta

# 2. Why these models are harder to fit ---------------------------------------------------

# random intercept for the 20 tree means: V = sigma_a^2 Z Z^T + sigma_e^2 I
y_tree <- tree_means$width
Z <- model.matrix(~ 0 + stand, data = tree_means)
V_ri <- function(sigma_a, sigma_e) sigma_a^2 * Z %*% t(Z) + sigma_e^2 * diag(20)
nll_ri <- function(par) -dmvnorm(y_tree, rep(par[1], 20), V_ri(exp(par[2]), exp(par[3])), log = TRUE)
fit_ri <- optim(c(mean(y_tree), log(0.1), log(0.1)), nll_ri, method = "BFGS")
c(mu = fit_ri$par[1], sigma_a = exp(fit_ri$par[2]), sigma_e = exp(fit_ri$par[3]))

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

# 3. Random slopes ------------------------------------------------------------------------

# drought sensitivity: 12 stands, records of 5-30 years; simulated
# truth: mean intercept 0.9, mean slope 0.25, sigma_a, sigma_b = (0.15, 0.10), rho_ab = -0.6, sigma_e = 0.12
set.seed(12)
n_years <- rep(c(30, 24, 18, 12, 8, 5), each = 2)        # record length per stand
sd_true <- c(0.15, 0.10)                                   # among-stand SDs: intercept, slope
rho_true <- -0.6
G_true <- diag(sd_true) %*% matrix(c(1, rho_true, rho_true, 1), 2) %*% diag(sd_true)
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

# 4. Mixed models in R --------------------------------------------------------------------

library(lme4)
fit_ri_lmer <- lmer(width ~ precip10 + (1 | stand), data = drought, REML = FALSE)
fit_rs_lmer <- lmer(width ~ precip10 + (precip10 | stand), data = drought, REML = FALSE)
fixef(fit_rs_lmer)
VarCorr(fit_rs_lmer)
AIC(fit_ri_lmer, fit_rs_lmer)
vc_rs <- as.data.frame(VarCorr(fit_rs_lmer))$sdcor        # sigma_a, sigma_b, rho_ab, sigma_e

fit_rs_lme <- lme(width ~ precip10, random = ~ precip10 | stand, data = drought, method = "ML")

library(glmmTMB)
fit_rs_tmb <- glmmTMB(width ~ precip10 + (precip10 | stand), data = drought)   # ML by default

# by hand: one MVN block per stand, V_j = Z_j G Z_j^T + sigma_e^2 I
rows_by_stand <- split(seq_len(nrow(drought)), drought$stand)
nll_rs <- function(par) {
  sds <- exp(par[3:4])
  rho_ab <- tanh(par[5])
  sigma_e <- exp(par[6])
  G <- matrix(c(sds[1]^2, rho_ab * sds[1] * sds[2], rho_ab * sds[1] * sds[2], sds[2]^2), 2)
  nll <- 0
  for (rows in rows_by_stand) {
    Zj <- cbind(1, drought$precip10[rows])
    Vj <- Zj %*% G %*% t(Zj) + sigma_e^2 * diag(length(rows))
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
colnames(compare) <- c("mu_a", "mu_b", "sigma_a", "sigma_b", "rho_ab", "sigma_e", "logLik")
round(compare, 4)

# Discussion: the four short-record stands alone
short <- droplevels(drought[drought$stand %in% c("S09", "S10", "S11", "S12"), ])
fit_short <- lmer(width ~ precip10 + (precip10 | stand), data = short, REML = FALSE)
VarCorr(fit_short)

# Bonus: AR(1) by generalized least squares ----------------------------------------------

# for fixed phi: beta = (X'R^-1 X)^-1 X'R^-1 y, sigma_e^2 = r'R^-1 r / n; then search over phi alone
X_a <- cbind(1, chron_a$precip, chron_a$year_c)
y_a <- chron_a$width
n_a <- length(y_a)

gls_profile <- function(phi) {
  R <- phi^lags / (1 - phi^2)
  R_inv <- solve(R)
  beta <- solve(t(X_a) %*% R_inv %*% X_a, t(X_a) %*% R_inv %*% y_a)
  r <- as.vector(y_a - X_a %*% beta)
  sigma2 <- as.numeric(t(r) %*% R_inv %*% r) / n_a
  list(beta = as.vector(beta), sigma = sqrt(sigma2),
       loglik = dmvnorm(y_a, as.vector(X_a %*% beta), sigma2 * R, log = TRUE))
}
phi_gls <- optimize(function(phi) -gls_profile(phi)$loglik, c(-0.95, 0.95))$minimum
phi_gls
gls_profile(phi_gls)          # matches optim() and gls(); gls() reports sigma_e / sqrt(1 - phi^2)

# profile likelihood and approximate 95% interval for phi
phi_seq <- seq(-0.5, 0.9, by = 0.005)
ll_seq <- sapply(phi_seq, function(r) gls_profile(r)$loglik)
range(phi_seq[ll_seq > max(ll_seq) - 1.92])
plot(phi_seq, ll_seq, type = "l", lwd = 2, xlab = "phi", ylab = "Profile log-likelihood")
abline(h = max(ll_seq) - 1.92, lty = 2)

# Bonus: a random intercept and AR(1) together -------------------------------------------

Z_chron <- model.matrix(~ 0 + stand, data = chron)
V_combined <- function(sigma_a, phi, sigma_e) {
  sigma_a^2 * Z_chron %*% t(Z_chron) + kronecker(diag(4), V_ar1(phi, sigma_e))
}
nll_combined <- function(par) {
  mu <- par[1] + par[2] * chron$precip
  V <- V_combined(exp(par[3]), tanh(par[4]), exp(par[5]))
  -dmvnorm(chron$width, mu, V, log = TRUE)
}
fit_comb <- optim(c(0.4, 0.02, log(0.1), 0, log(0.1)), nll_combined, method = "BFGS",
                  control = list(maxit = 1000))
c(mu = fit_comb$par[1], beta = fit_comb$par[2], sigma_a = exp(fit_comb$par[3]),
  phi = tanh(fit_comb$par[4]), sigma_e = exp(fit_comb$par[5]))

fit_lme <- lme(width ~ precip, random = ~ 1 | stand, data = chron,
               correlation = corAR1(form = ~ year | stand), method = "ML")
fixef(fit_lme)
coef(fit_lme$modelStruct$corStruct, unconstrained = FALSE)
c(-fit_comb$value, logLik(fit_lme))       # same maximized log-likelihood

# Bonus: a Poisson GLMM by hand ----------------------------------------------------------

# marginal moments: mu = exp(b0 + b1 x + sigma_a^2 / 2); Var = mu + mu^2 (exp(sigma_a^2) - 1)
mu_glmm <- exp(1.5 + 0.5^2 / 2)               # vigor 0, beta0 = 1.5, sigma_a = 0.5
moments <- c(mean = mu_glmm, var = mu_glmm + mu_glmm^2 * (exp(0.5^2) - 1), cov = mu_glmm^2 * (exp(0.5^2) - 1))

set.seed(746)
a_sim <- rnorm(1e5, 0, 0.5)                  # 100,000 stands, two trees each
y_tree1 <- rpois(1e5, exp(1.5 + a_sim))
y_tree2 <- rpois(1e5, exp(1.5 + a_sim))
rbind(formula = moments, simulation = c(mean(y_tree1), var(y_tree1), cov(y_tree1, y_tree2)))
moments[3] / moments[2]                  # within-stand correlation

# cones on 8 trees per stand; log(lambda) = b0 + b1 * vigor + b_j; simulated
set.seed(12)
cones <- data.frame(stand = factor(rep(stand_names, each = 8)),
                    vigor = round(rnorm(96), 2))
a_stand <- rnorm(12, 0, 0.5)
cones$count <- rpois(96, exp(1.5 + 0.4 * cones$vigor + a_stand[as.integer(cones$stand)]))

# each stand's likelihood integrates over its random intercept a
nll_glmm <- function(par) {
  sigma_a <- exp(par[3])
  nll <- 0
  for (g in split(cones, cones$stand)) {
    stand_lik <- function(a) {
      sapply(a, function(ai) exp(sum(dpois(g$count, exp(par[1] + par[2] * g$vigor + ai), log = TRUE)))) *
        dnorm(a, 0, sigma_a)
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
# moved to lecture_notes/CODE_Week13_HierarchicalModelsAndGAMs.R (Week 13 page)
