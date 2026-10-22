
#  NRES 746, Week 10: LRT, AIC, penalized likelihood (Mon.); Bayesian model comparison (Wed.) =
#  R demos only; section numbers match the NOTES files. Needs emdbook, MASS, glmnet.

# MONDAY, Oct. 26: Likelihood ratio tests, AIC, ridge and lasso ===============

rm(list=ls())

# 3. Likelihood ratio test (~8 min) -----------------------------------------------

library(emdbook)
fir <- na.omit(FirDBHFec[, c("TOTCONES", "DBH", "WAVE_NON")])   # balsam fir, 242 trees
fir$TOTCONES <- round(fir$TOTCONES)
wave <- as.numeric(fir$WAVE_NON)   # 1 = non-wave, 2 = wave

plot(TOTCONES ~ DBH, data = fir, col = wave, pch = 19,
     xlab = "DBH (cm)", ylab = "Total cones")
legend("topleft", c("non-wave", "wave"), col = 1:2, pch = 19, bty = "n")

# cones ~ NegBin(mu = a * DBH^b, size = k); a, b, k can be vectors (one value per tree)
nll_nb <- function(a, b, k) {
  -sum(dnbinom(fir$TOTCONES, mu = a * fir$DBH^b, size = k, log = TRUE))
}

# full model: separate a, b, k for wave and non-wave sites (6 parameters)
nll_full <- function(p) nll_nb(p[1:2][wave], p[3:4][wave], p[5:6][wave])
# no-wave model: one a, b, k for all trees (3 parameters)
nll_nowave <- function(p) nll_nb(p[1], p[2], p[3])

# Nelder-Mead to get close, then BFGS to polish
fit_model <- function(nll, start) {
  fit <- optim(start, nll, control = list(maxit = 10000))
  optim(fit$par, nll, method = "BFGS", control = list(maxit = 10000))
}
fit_full <- fit_model(nll_full, c(1, 1, 2, 2, 1, 1))
fit_nowave <- fit_model(nll_nowave, c(1, 2, 1))
c(full = fit_full$value, nowave = fit_nowave$value)   # min NLL

# deviance = 2 * (NLL_reduced - NLL_full) ~ chi-square, df = difference in parameters
deviance <- 2 * (fit_nowave$value - fit_full$value)
deviance
qchisq(0.95, df = 3)                            # critical value
pchisq(deviance, df = 3, lower.tail = FALSE)    # p-value

curve(dchisq(x, df = 3), 0, 12, lwd = 2, xlab = "Deviance", ylab = "Density",
      main = "Deviance under the reduced model (df = 3)")
abline(v = qchisq(0.95, 3), lty = 2)
abline(v = deviance, col = "red", lwd = 3)

# Poisson, no wave effect (2 parameters)
nll_pois <- function(p) -sum(dpois(fir$TOTCONES, p[1] * fir$DBH^p[2], log = TRUE))
fit_pois <- fit_model(nll_pois, c(1, 2))
2 * (fit_pois$value - fit_nowave$value)
#  [PAUSE] is the Poisson model nested in the negative binomial? can we use an LRT?

# 4. AIC (~10 min) -------------------------------------------------------------

# two intermediate models: wave changes only a and k (5 parameters), or only b (4)
nll_constb <- function(p) nll_nb(p[1:2][wave], p[3], p[4:5][wave])
nll_constak <- function(p) nll_nb(p[1], p[2:3][wave], p[4])
fit_constb <- fit_model(nll_constb, c(1, 1, 2, 1, 1))
fit_constak <- fit_model(nll_constak, c(1, 2, 2, 1))

fits <- list(full = fit_full, const_b = fit_constb, const_ak = fit_constak,
             no_wave = fit_nowave, poisson = fit_pois)
n <- nrow(fir)
aic_tab <- data.frame(model = names(fits),
                      k = sapply(fits, function(f) length(f$par)),
                      nll = sapply(fits, function(f) f$value))
aic_tab$AIC <- 2 * aic_tab$nll + 2 * aic_tab$k                  # -2 log L + 2k
aic_tab$dAIC <- aic_tab$AIC - min(aic_tab$AIC)
aic_tab$weight <- exp(-aic_tab$dAIC / 2) / sum(exp(-aic_tab$dAIC / 2))
aic_tab$AICc <- aic_tab$AIC + 2 * aic_tab$k * (aic_tab$k + 1) / (n - aic_tab$k - 1)
aic_tab$BIC <- 2 * aic_tab$nll + log(n) * aic_tab$k             # -2 log L + k log(n)
aic_tab <- aic_tab[order(aic_tab$AIC), ]
rownames(aic_tab) <- NULL
cbind(aic_tab[1:2], round(aic_tab[-(1:2)], 2))
#  [PAUSE] the no-wave model has weight 0.52. is there a 52% chance it's the true model?

# 5. Penalized likelihood: ridge and lasso (~10 min) ----------------------------

library(MASS)
library(glmnet)

# brown creeper counts at 50 point-count stations, 10 standardized habitat covariates.
# Simulated, so we know the truth: only big_trees, snags, and road_dist matter.
covars <- c("canopy", "big_trees", "snags", "shrub", "elev", "slope",
            "north", "road_dist", "stream_dist", "grazing")
beta_true <- c(0, 0.5, 0.4, 0, 0, 0, 0, 0.3, 0, 0)
cov_mat <- diag(10)
cov_mat[1, 2] <- cov_mat[2, 1] <- 0.6    # canopy ~ big trees
cov_mat[1, 4] <- cov_mat[4, 1] <- -0.4   # canopy ~ shrub
cov_mat[5, 6] <- cov_mat[6, 5] <- 0.5    # elevation ~ slope
set.seed(19)
X <- mvrnorm(50, rep(0, 10), cov_mat)
colnames(X) <- covars
creepers <- rpois(50, exp(log(2) + X %*% beta_true))   # log(mu) = log(2) + X beta

# ordinary MLE with every covariate: no estimate is exactly 0
fit_glm <- glm(creepers ~ X, family = poisson)
round(coef(fit_glm)[-1], 2)

# all 2^10 = 1024 subsets, ranked by AIC
subsets <- expand.grid(rep(list(c(FALSE, TRUE)), 10))
subset_aic <- apply(subsets, 1, function(use) {
  if (!any(use)) return(AIC(glm(creepers ~ 1, family = poisson)))
  AIC(glm(creepers ~ X[, use, drop = FALSE], family = poisson))
})
covars[unlist(subsets[which.min(subset_aic), ])]   # best subset
sum(subset_aic - min(subset_aic) < 2)              # models within 2 AIC units

# penalized NLL: NLL(beta) + lambda * sum(beta^2) (ridge) or lambda * sum(|beta|) (lasso)
fit_ridge <- glmnet(X, creepers, family = "poisson", alpha = 0)
fit_lasso <- glmnet(X, creepers, family = "poisson", alpha = 1)
par(mfrow = c(1, 2))
plot(fit_ridge, xvar = "lambda", label = TRUE, main = "Ridge")
plot(fit_lasso, xvar = "lambda", label = TRUE, main = "Lasso")
par(mfrow = c(1, 1))
#  [PAUSE] read the lasso panel right to left: which covariates enter first?

# choose lambda by 10-fold cross-validation: predict each held-out tenth
set.seed(746)
cv_lasso <- cv.glmnet(X, creepers, family = "poisson", alpha = 1)
plot(cv_lasso)   # left dashed line: lambda.min; right: lambda.1se (simplest within 1 SE)
coef_lasso <- coef(cv_lasso, s = "lambda.1se")[-1]

round(data.frame(truth = beta_true,
                 MLE = coef(fit_glm)[-1],
                 AIC_best = ifelse(unlist(subsets[which.min(subset_aic), ]), 1, 0),
                 lasso = coef_lasso,
                 row.names = covars), 2)
#  [PAUSE] the lasso kept road_dist but estimates 0.07; the truth is 0.3. is the lasso broken?





# WEDNESDAY, Oct. 28: Bayesian model comparison and regularization ==============
#  (placeholder: to be built)
