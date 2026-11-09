
#  NRES 746, Week 11: linear models as likelihoods, linear algebra, modeling variance =========
#  R demos only; sections match the Week 11 web page. Needs mvtnorm, lme4.

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

# classic ANCOVA tests: F-tests between nested models
fit_precip <- lm(width ~ precip, data = chron)
fit_slopes <- lm(width ~ 0 + stand * precip, data = chron)   # one slope per stand
anova(fit_precip, fit_lm)    # do stands differ, given precipitation?
anova(fit_lm, fit_slopes)    # same precipitation slope in every stand?

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

library(lme4)          # syntax for nlme, lme4, glmmTMB: lecture_notes/NOTES_Reference_MixedModelSyntax.pdf
fit_lmer <- lmer(width ~ 1 + (1 | stand), data = tree_means, REML = FALSE)
fixef(fit_lmer)
as.data.frame(VarCorr(fit_lmer))[, c("grp", "sdcor")]

sb_hat <- exp(fit_ri$par[2])
s_hat <- exp(fit_ri$par[3])
sb_hat^2 / (sb_hat^2 + s_hat^2)       # intraclass correlation

#  [PAUSE] The among-stand SD comes from only four stands. How much would you trust it?

# 4. Autocorrelation and fitting --------------------------------------------------------

# AR(1), why structured-variance models are harder to fit, and the combined
# random intercept + AR(1) bonus moved to CODE_Week12_MixedModels.R
