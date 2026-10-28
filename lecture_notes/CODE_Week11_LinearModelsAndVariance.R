
#  NRES 746, Week 11: linear models as likelihoods, linear algebra, modeling variance =========
#  R demos only; sections match the Week 11 web page. Needs mvtnorm, lme4, nlme.

# STANDARD STATISTICS AS LIKELIHOOD MODELS ============================================

rm(list = ls())

# 3. One model behind lm() ------------------------------------------------------------

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

# stand chronologies: mean of the 5 trees in each stand and year
chron <- aggregate(width ~ stand + year + precip, data = rings, FUN = mean)
chron <- chron[order(chron$stand, chron$year), ]
head(chron)

plot(width ~ precip, data = chron, col = as.integer(factor(chron$stand)), pch = 19,
     xlab = "Winter precipitation (cm)", ylab = "Ring width (mm)")

# ANCOVA: one intercept per stand (a[1:4]), shared slope b; sigma on the log scale
nll_ancova <- function(par, dat) {
  a <- par[1:4]
  b <- par[5]
  sigma <- exp(par[6])
  mu <- a[as.integer(dat$stand)] + b * dat$precip
  -sum(dnorm(dat$width, mu, sigma, log = TRUE))
}
chron$stand <- factor(chron$stand)
fit_ancova <- optim(c(rep(0.5, 4), 0, log(0.1)), nll_ancova, dat = chron,
                    method = "BFGS", control = list(maxit = 1000))
fit_lm <- lm(width ~ 0 + stand + precip, data = chron)   # 0 + : one intercept per stand

cbind(optim = fit_ancova$par[1:5], lm = coef(fit_lm))
c(optim = exp(fit_ancova$par[6]), lm = summary(fit_lm)$sigma)   # MLE divides by n; lm by n - 5

#  [PAUSE] Each precipitation value is shared by all four stands in a year. Does that matter for independence?

# 4. The model matrix -----------------------------------------------------------------

rows_1996 <- chron$year == 1996
model.matrix(~ 0 + stand + precip, data = chron)[rows_1996, ]   # one column per stand
model.matrix(~ stand + precip, data = chron)[rows_1996, ]       # R's default coding
coef(lm(width ~ stand + precip, data = chron))

#  [PAUSE] In a model with an intercept, why can't R include indicator columns for all four stands?

# 5. Generalized linear models --------------------------------------------------------

tree_means <- aggregate(width ~ stand + tree_id, data = rings, FUN = mean)
set.seed(5)
tree_means$cones <- rpois(20, exp(-0.5 + 3 * tree_means$width))

# Poisson GLM by hand: log(mu) = a + b * width
nll_pois <- function(par) {
  -sum(dpois(tree_means$cones, exp(par[1] + par[2] * tree_means$width), log = TRUE))
}
fit_pois <- optim(c(0, 1), nll_pois, method = "BFGS")
fit_glm <- glm(cones ~ width, family = poisson, data = tree_means)
cbind(optim = fit_pois$par, glm = coef(fit_glm))

# Discussion: lag-1 correlation of ANCOVA residuals within stands
mean(sapply(split(resid(fit_lm), chron$stand), function(r) cor(r[-1], r[-length(r)])))

# Bonus: nonlinear least squares ------------------------------------------------------

fit_nls <- nls(width ~ a * precip / (b + precip), data = chron,
               start = list(a = 1.5, b = 20))
coef(fit_nls)


# LINEAR ALGEBRA FOR STATISTICS =======================================================

# 1. Vectors and matrices -------------------------------------------------------------

A <- matrix(c(2, 1, 0, 3, -1, 4), nrow = 2)
A
dim(A)
A[2, 3]     # row 2, column 3

# 2. Matrix multiplication ------------------------------------------------------------

X <- model.matrix(~ stand + precip, data = chron)
beta_hat <- coef(lm(width ~ stand + precip, data = chron))
mu_hat <- X %*% beta_hat               # 100 x 5 times 5 x 1 = 100 x 1
all.equal(as.vector(mu_hat), unname(fitted(lm(width ~ stand + precip, data = chron))))

# 3. Transpose, identity, and inverse -------------------------------------------------

# normal equations: t(X) X beta = t(X) y
XtX <- t(X) %*% X
XtX                                      # 5 x 5
beta_ne <- solve(XtX, t(X) %*% chron$width)   # solves XtX beta = Xty
cbind(normal_eq = as.vector(beta_ne), lm = beta_hat)

#  [PAUSE] Suppose X also had an indicator column for stand A. What would go wrong with solve(XtX)?

# 4. Determinant and quadratic form ---------------------------------------------------

# with V = sigma^2 I, the MVN is a product of univariate normals
library(mvtnorm)
r <- resid(fit_lm)
sigma_mle <- sqrt(mean(r^2))
c(sum_dnorm = sum(dnorm(r, 0, sigma_mle, log = TRUE)),
  dmvnorm = dmvnorm(r, rep(0, 100), sigma_mle^2 * diag(100), log = TRUE))

# 5. Valid covariance matrices --------------------------------------------------------

V_ok <- matrix(c(4, 1, 1, 2), 2)       # variances 4 and 2, covariance 1
eigen(V_ok)$values                     # all positive: positive definite


# MODELING VARIANCE: RANDOM EFFECTS AND AUTOCORRELATION ===============================

# 3. A random intercept ---------------------------------------------------------------

# V = sigma_b^2 Z Z^T + sigma^2 I: block diagonal, one 5 x 5 block per stand
y_tree <- tree_means$width
Z <- model.matrix(~ 0 + stand, data = tree_means)
ZZt <- Z %*% t(Z)
V_ri <- function(sigma_b, sigma) sigma_b^2 * ZZt + sigma^2 * diag(20)
V_ri(0.12, 0.1)[1:7, 1:7]     # first stand's block, and the start of the next

nll_ri <- function(par) {
  V <- V_ri(exp(par[2]), exp(par[3]))
  -dmvnorm(y_tree, rep(par[1], 20), V, log = TRUE)
}
fit_ri <- optim(c(mean(y_tree), log(0.1), log(0.1)), nll_ri, method = "BFGS")
c(mu = fit_ri$par[1], sigma_b = exp(fit_ri$par[2]), sigma = exp(fit_ri$par[3]))

library(lme4)
fit_lmer <- lmer(width ~ 1 + (1 | stand), data = tree_means, REML = FALSE)
fixef(fit_lmer)
as.data.frame(VarCorr(fit_lmer))[, c("grp", "sdcor")]

sb_hat <- exp(fit_ri$par[2])
s_hat <- exp(fit_ri$par[3])
sb_hat^2 / (sb_hat^2 + s_hat^2)       # intraclass correlation

#  [PAUSE] The among-stand SD comes from only four stands. How much would you trust it?

# 4. Autocorrelation: AR(1) -----------------------------------------------------------

chron_a <- chron[chron$stand == "A", ]
chron_a$year_c <- chron_a$year - 2008          # centered year: 0 = 2008
fit_lm_a <- lm(width ~ precip + year_c, data = chron_a)
summary(fit_lm_a)$coefficients
par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 1))
plot(chron_a$year, resid(fit_lm_a), type = "b", pch = 19, xlab = "Year", ylab = "Residual")
abline(h = 0, lty = 2)
acf(resid(fit_lm_a), main = "", lag.max = 10)
par(mfrow = c(1, 1))

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
coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)   # rho

# trend: estimate, SE, p under independence vs. AR(1) (Hessian SE)
se_ar1 <- sqrt(diag(solve(fit_ar1$hessian)))
summary(fit_lm_a)$coefficients["year_c", ]
c(estimate = fit_ar1$par[3], se = se_ar1[3], p = 2 * pnorm(-abs(fit_ar1$par[3] / se_ar1[3])))

#  [PAUSE] With rho near 0.9, what does one AR(1) series look like? How would you tell it from a slow trend?

# 5. Why these models are harder to fit -----------------------------------------------

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

# Bonus: a random intercept and AR(1) together ----------------------------------------

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
