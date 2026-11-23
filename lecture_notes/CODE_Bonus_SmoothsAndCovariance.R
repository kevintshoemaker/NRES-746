# Bonus: Smooths and correlated errors, two views of one model (NRES 746)
# Packages: mvtnorm, lme4, nlme, mgcv

library(mvtnorm)
library(lme4)
library(nlme)
library(mgcv)

##########
# Two ways to write the same model
##########

# 2. A random intercept, both ways -------------------------------------

set.seed(7)
pond <- factor(rep(LETTERS[1:4], each = 5))
a_true <- rnorm(4, 0, 0.8)
load <- round(3 + a_true[pond] + rnorm(20, 0, 0.5), 1)   # log Bd load
tapply(load, pond, mean)

# marginal view: y ~ MVN(mu, V), V = sigma_a^2 Z Z' + sigma_e^2 I
Z <- model.matrix(~ pond - 1)
nll_ri <- function(p) {
  V <- exp(2 * p[2]) * Z %*% t(Z) + exp(2 * p[3]) * diag(20)
  -dmvnorm(load, rep(p[1], 20), V, log = TRUE)
}
fit_ri <- optim(c(3, 0, 0), nll_ri, method = "BFGS")
mu_hat <- fit_ri$par[1]
sig_a <- exp(fit_ri$par[2])
sig_e <- exp(fit_ri$par[3])
c(mu = mu_hat, sigma_a = sig_a, sigma_e = sig_e, logLik = -fit_ri$value)

# conditional view: ridge on the pond effects, lambda = sigma_e^2 / sigma_a^2
V_ri <- sig_a^2 * Z %*% t(Z) + sig_e^2 * diag(20)
lambda_ri <- sig_e^2 / sig_a^2
C <- cbind(1, Z)
P <- diag(c(0, rep(lambda_ri, 4)))
est <- solve(t(C) %*% C + P, t(C) %*% load)
blup <- sig_a^2 * t(Z) %*% solve(V_ri, load - mu_hat)   # G Z' V^-1 (y - X beta)
m_ri <- lmer(load ~ 1 + (1 | pond), REML = FALSE)
cbind(conditional = est[-1], marginal = drop(blup), lmer = ranef(m_ri)$pond[, 1])

# 3. Smooths are random effects ----------------------------------------

set.seed(12)
n <- 40
year <- 1985:2024
yr <- year - 1985
temp <- round(4 + 0.06 * yr + rnorm(n, 0, 0.7), 1)   # March mean temperature
D <- abs(outer(yr, yr, "-"))
field <- drop(t(chol(25 * exp(-D / 6))) %*% rnorm(n))
day <- round(120 - 3 * temp + field + rnorm(n, 0, 2))   # first flowering, day of year
phen <- data.frame(year, yr, temp, day)

m_lm <- lm(day ~ temp, data = phen)
summary(m_lm)$coef

g_ml <- gam(day ~ temp + s(yr, bs = "ps", k = 20), data = phen, method = "ML")
c(temp = coef(g_ml)[["temp"]], lambda = g_ml$sp, edf = sum(g_ml$edf[-(1:2)]))

# split the smooth into unpenalized (fixed) and penalized (iid random) parts
sm <- smoothCon(s(yr, bs = "ps", k = 20), data = phen, absorb.cons = TRUE)[[1]]
re <- smooth2random(sm, names(phen), type = 2)
Xf <- re$Xf
Xr <- re$rand[[1]]
c(ncol(Xf), ncol(Xr))

# marginal fit: V = sigma_b^2 Xr Xr' + sigma_e^2 I, beta by GLS
X_gam <- cbind(1, phen$temp, Xf)
nll_gam <- function(q) {
  V <- exp(2 * q[1]) * Xr %*% t(Xr) + exp(2 * q[2]) * diag(n)
  b <- solve(t(X_gam) %*% solve(V, X_gam), t(X_gam) %*% solve(V, phen$day))
  -dmvnorm(phen$day, drop(X_gam %*% b), V, log = TRUE)
}
fit_gm <- optim(c(0, log(2)), nll_gam, method = "BFGS")
sd_gm <- exp(fit_gm$par)
V_gm <- sd_gm[1]^2 * Xr %*% t(Xr) + sd_gm[2]^2 * diag(n)
b_gm <- solve(t(X_gam) %*% solve(V_gm, X_gam), t(X_gam) %*% solve(V_gm, phen$day))
c(temp = b_gm[2], lambda = sd_gm[2]^2 / sd_gm[1]^2)   # matches g_ml
#  [PAUSE] With a first-difference penalty, which part of the smooth moves between the mean and the covariance?

# 4. Covariance functions are smooths ----------------------------------

# V = sigma_s^2 exp(-|s - t| / ell) + sigma_e^2 I
X <- cbind(1, phen$temp)
nll_exp <- function(p) {
  V <- exp(2 * p[3]) * exp(-D / exp(p[4])) + exp(2 * p[5]) * diag(n)
  -dmvnorm(phen$day, drop(X %*% p[1:2]), V, log = TRUE)
}
fit_exp <- optim(c(coef(m_lm), log(4), log(5), log(2)), nll_exp,
                 method = "BFGS", hessian = TRUE)
pe <- fit_exp$par
sig_s <- exp(pe[3]); ell <- exp(pe[4]); sig_e2 <- exp(pe[5])
se_exp <- sqrt(diag(solve(fit_exp$hessian)))
c(temp = pe[[2]], SE = se_exp[[2]], sigma_s = sig_s, ell = ell, sigma_e = sig_e2,
  logLik = -fit_exp$value)
exp(-1 / ell)   # AR(1) phi at annual spacing

m_gls <- gls(day ~ temp, data = phen, method = "ML",
             correlation = corExp(form = ~ yr, nugget = TRUE))
logLik(m_gls)

# K = U Lambda U' = B B', B = U Lambda^(1/2)
K <- sig_s^2 * exp(-D / ell)
eg <- eigen(K, symmetric = TRUE)
B <- eg$vectors %*% diag(sqrt(eg$values))
share <- cumsum(eg$values) / sum(eg$values)
share[c(5, 10)]

par(mfrow = c(2, 3))
for (k in 1:6) {
  plot(year, B[, k], type = "l", lwd = 2, ylim = range(B[, 1:6]),
       xlab = "Year", ylab = "", main = paste("Function", k))
  abline(h = 0, col = "grey80")
}
par(mfrow = c(1, 1))

# conditional fit with gamma ~ N(0, I): penalty sigma_e^2 gamma'gamma
C <- cbind(X, B)
P <- diag(c(0, 0, rep(sig_e2^2, n)))
est_c <- solve(t(C) %*% C + P, t(C) %*% phen$day)
V_exp <- K + sig_e2^2 * diag(n)
b_gls <- solve(t(X) %*% solve(V_exp, X), t(X) %*% solve(V_exp, phen$day))
krig <- K %*% solve(V_exp, phen$day - X %*% b_gls)
c(temp_conditional = est_c[2], temp_marginal = b_gls[2],
  max_curve_difference = max(abs(B %*% est_c[-(1:2)] - krig)))

plot(year, phen$day - drop(X %*% b_gls), pch = 16, col = "grey55",
     xlab = "Year", ylab = "Flowering date minus temperature effect (days)")
lines(year, krig, lwd = 4, col = "grey75")
lines(year, B %*% est_c[-(1:2)], lwd = 1.5, lty = 2)

# low-rank versions: keep the first k basis functions
low_rank <- t(sapply(c(3, 5, 10, 20, 40), function(k) {
  Bk <- B[, seq_len(k), drop = FALSE]
  Vk <- Bk %*% t(Bk) + sig_e2^2 * diag(n)
  bk <- solve(t(X) %*% solve(Vk, X), t(X) %*% solve(Vk, phen$day))
  c(k = k, variance_share = share[k], temp = bk[2],
    logLik = dmvnorm(phen$day, drop(X %*% bk), Vk, log = TRUE))
}))
round(low_rank, 3)
#  [PAUSE] Which kernel needs more basis functions to reproduce its covariance, and why?

##########
# Where the two views part ways
##########

# 5. The unpenalized part of a smooth ----------------------------------

X3 <- cbind(1, phen$temp, phen$yr)
nll_exp3 <- function(p) {
  V <- exp(2 * p[4]) * exp(-D / exp(p[5])) + exp(2 * p[6]) * diag(n)
  -dmvnorm(phen$day, drop(X3 %*% p[1:3]), V, log = TRUE)
}
fit_exp3 <- optim(c(coef(lm(day ~ temp + yr, data = phen)), log(4), log(5), log(2)),
                  nll_exp3, method = "BFGS", hessian = TRUE)
se_exp3 <- sqrt(diag(solve(fit_exp3$hessian)))
c(temp = fit_exp3$par[[2]], trend = fit_exp3$par[[3]], trend_SE = se_exp3[[3]],
  logLik = -fit_exp3$value)

g_reml <- gam(day ~ temp + s(yr, bs = "ps", k = 20), data = phen, method = "REML")
g_fx <- gam(day ~ temp + s(yr, k = 10, fx = TRUE), data = phen)   # unpenalized
rbind(lm = c(coef(m_lm)[2], summary(m_lm)$coef[2, 2]),
      exponential = c(pe[2], se_exp[2]),
      exponential_trend = c(fit_exp3$par[2], se_exp3[2]),
      gam_reml = c(coef(g_reml)[2], summary(g_reml)$se[2]),
      gam_fixed = c(coef(g_fx)[2], summary(g_fx)$se[2]))

# 6. Non-normal responses ----------------------------------------------

# E[exp(eta + s)] = exp(eta) exp(sigma_s^2 / 2)
set.seed(746)
s_draw <- rnorm(1e6, 0, 0.8)
c(at_s_zero = exp(log(5)), marginal_simulated = mean(exp(log(5) + s_draw)),
  marginal_formula = 5 * exp(0.8^2 / 2))
