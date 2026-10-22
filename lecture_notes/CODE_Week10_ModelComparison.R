
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





# WEDNESDAY, Oct. 28: Bayesian model comparison and shrinkage priors ===========
#  Uses fir, wave, aic_tab, X, creepers, covars, beta_true, coef_lasso from Monday's section.

# 3. WAIC and LOO (~12 min) -------------------------------------------------------

options(cmdstanr_no_ver_check = TRUE)
library(cmdstanr)
library(loo)

# fir NegBin power law; G = 1: one a, b, k for all trees; G = 2: separate for wave/non-wave
fir_code <- "
data {
  int<lower=1> N;
  int<lower=1> G;
  array[N] int<lower=0> cones;
  vector<lower=0>[N] dbh;
  array[N] int<lower=1, upper=G> grp;
}
parameters {
  vector[G] log_a;
  vector[G] b;
  vector<lower=0>[G] k;
}
model {
  log_a ~ normal(0, 5);
  b ~ normal(0, 5);
  k ~ exponential(0.1);
  cones ~ neg_binomial_2(exp(log_a[grp]) .* pow(dbh, b[grp]), k[grp]);
}
generated quantities {
  vector[N] log_lik;    // one value per tree, for WAIC and LOO
  for (i in 1:N)
    log_lik[i] = neg_binomial_2_lpmf(cones[i] | exp(log_a[grp[i]]) * dbh[i]^b[grp[i]], k[grp[i]]);
}
"
fir_model <- cmdstan_model(write_stan_file(fir_code))
fit_fir <- function(G) {
  grp <- if (G == 1) rep(1L, nrow(fir)) else wave
  fir_model$sample(data = list(N = nrow(fir), G = G, cones = fir$TOTCONES, dbh = fir$DBH,
                               grp = grp),
                   chains = 4, parallel_chains = 4, seed = 746, refresh = 0)
}
fit_nowave_b <- fit_fir(1)
fit_full_b <- fit_fir(2)

# DIC (Bolker 6.6.3): pD = mean deviance - deviance at posterior means; DIC = Dhat + 2 pD
dic <- function(fit, G) {
  grp <- if (G == 1) rep(1L, nrow(fir)) else wave
  d_bar <- mean(-2 * rowSums(fit$draws("log_lik", format = "matrix")))
  pm <- function(v) colMeans(fit$draws(v, format = "matrix"))
  d_hat <- -2 * sum(dnbinom(fir$TOTCONES, mu = exp(pm("log_a")[grp]) * fir$DBH^pm("b")[grp],
                            size = pm("k")[grp], log = TRUE))
  c(DIC = d_hat + 2 * (d_bar - d_hat), pD = d_bar - d_hat)
}
waic_nowave <- waic(fit_nowave_b$draws("log_lik", format = "matrix"))
waic_full <- waic(fit_full_b$draws("log_lik", format = "matrix"))
loo_nowave <- fit_nowave_b$loo()
loo_full <- fit_full_b$loo()
loo_nowave    # elpd_loo, p_loo, looic, and the Pareto k check

# the same two models, four ways (lower is better; penalties: k, pD, p_waic, p_loo)
ic_tab <- rbind(
  AIC = aic_tab$AIC[match(c("no_wave", "full"), aic_tab$model)],
  DIC = c(dic(fit_nowave_b, 1)["DIC"], dic(fit_full_b, 2)["DIC"]),
  WAIC = c(waic_nowave$estimates["waic", 1], waic_full$estimates["waic", 1]),
  LOOIC = c(loo_nowave$estimates["looic", 1], loo_full$estimates["looic", 1]))
colnames(ic_tab) <- c("no_wave", "full")
round(ic_tab, 1)
loo_compare(list(no_wave = loo_nowave, full = loo_full))
#  [PAUSE] the full model is 2.1 elpd behind, SE 1.8. what can you conclude?

# 4. Shrinkage priors (~12 min) ----------------------------------------------------

# creepers again: Poisson regression; prior on the slopes is Normal(0, s) or Laplace(0, s)
pois_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=0> y;
  int<lower=1, upper=2> prior;    // 1 = Normal, 2 = Laplace
  real<lower=0> s;
}
parameters {
  real b0;
  vector[K] beta;
}
model {
  b0 ~ normal(0, 5);
  if (prior == 1) beta ~ normal(0, s);
  else beta ~ double_exponential(0, s);
  y ~ poisson_log(b0 + X * beta);
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = poisson_log_lpmf(y[i] | b0 + X[i] * beta);
}
"
pois_model <- cmdstan_model(write_stan_file(pois_code))
fit_creepers <- function(prior, s) {
  pois_model$sample(data = list(N = 50, K = 10, X = X, y = creepers, prior = prior, s = s),
                    chains = 4, parallel_chains = 4, seed = 746, refresh = 0)
}
fits_shrink <- list(vague = fit_creepers(1, 10), normal = fit_creepers(1, 0.2),
                    laplace = fit_creepers(2, 0.15))

# posterior medians and 90% intervals for each slope, under each prior; x = truth
post_sum <- lapply(fits_shrink, function(f) {
  as.matrix(f$summary("beta", "median", ~quantile(.x, c(0.05, 0.95)))[, -1])
})
par(mar = c(4, 7, 1, 1))
plot(NA, xlim = c(-0.5, 1), ylim = c(0.5, 10.5), yaxt = "n", xlab = "Slope (per SD)", ylab = "")
axis(2, at = 10:1, labels = covars, las = 1, cex.axis = 0.8)
abline(v = 0, col = "grey70")
offsets <- c(vague = 0.25, normal = 0, laplace = -0.25)
for (p in names(post_sum)) {
  yy <- (10:1) + offsets[p]
  segments(post_sum[[p]][, 2], yy, post_sum[[p]][, 3], yy, lwd = 2, col = match(p, names(post_sum)))
  points(post_sum[[p]][, 1], yy, pch = 19, col = match(p, names(post_sum)))
}
points(beta_true, 10:1, pch = 4, cex = 1.3, lwd = 2)
points(coef_lasso, 10:1, pch = 2)
legend("bottomright", c("vague prior", "Normal(0, 0.2)", "Laplace(0, 0.15)", "truth", "glmnet lasso"),
       col = c(1:3, 1, 1), pch = c(19, 19, 19, 4, 2), bty = "n", cex = 0.8)
par(mar = c(5, 4, 4, 2) + 0.1)
#  [PAUSE] lasso set northness to exactly 0; the Laplace posterior says 0.08 (-0.07 to 0.28). which is more honest?

sapply(fits_shrink, function(f) f$loo()$estimates["p_loo", 1])   # effective parameters shrink
loo_compare(lapply(fits_shrink, function(f) f$loo()))

# prior shapes: Normal, Laplace, horseshoe (simulated: Normal with half-Cauchy scale)
set.seed(746)
hs <- rnorm(1e6, 0, 0.1 * abs(rcauchy(1e6)))
curve(dnorm(x, 0, 0.2), -1, 1, n = 501, lwd = 2, ylim = c(0, 6), xlab = "Slope", ylab = "Prior density")
curve(dexp(abs(x), 1 / 0.15) / 2, -1, 1, n = 501, add = TRUE, lwd = 2, col = 3)
lines(density(hs[abs(hs) < 1], bw = 0.01, n = 1024), lwd = 2, col = 4)
legend("topright", c("Normal (ridge)", "Laplace (lasso)", "horseshoe"), col = c(1, 3, 4),
       lwd = 2, bty = "n")
