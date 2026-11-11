# NRES 746, Bonus: Evaluating model performance (code from the web page)
# Packages: emdbook, cmdstanr

# WHAT MAKES A MODEL GOOD? =================================================

# 2. Example: myxomatosis titers ---------------------------------------------

library(emdbook)
myx <- subset(MyxoTiter_sum, grade == 1)   # 27 rabbits

ricker <- function(a, b, t) a * t * exp(-b * t)

# Gamma NLL with Ricker mean; 'd' lets us refit to subsets later
nll_myx <- function(p, d = myx) {
  if (any(p <= 0)) return(Inf)
  mu <- ricker(p[1], p[2], d$day)
  -sum(dgamma(d$titer, shape = p[3], scale = mu / p[3], log = TRUE))
}

fit <- optim(c(a = 1, b = 0.2, shape = 50), nll_myx)
mle <- fit$par
round(mle, 3)

# GOODNESS OF FIT ===========================================================

# 1. Plug-in prediction bounds -------------------------------------------------

days <- seq(1, 10, by = 0.1)
mu_days <- ricker(mle["a"], mle["b"], days)
upper <- qgamma(0.975, shape = mle["shape"], scale = mu_days / mle["shape"])
lower <- qgamma(0.025, shape = mle["shape"], scale = mu_days / mle["shape"])

plot(titer ~ day, data = myx, pch = 16, xlim = c(1, 10), ylim = c(0, 12),
     xlab = "Days since infection", ylab = "Titer")
lines(days, mu_days, lwd = 2)
lines(days, upper, lty = 2)
lines(days, lower, lty = 2)

# 2. Parametric bootstrap ------------------------------------------------------

#  [PAUSE] How does a parametric bootstrap differ from the ordinary bootstrap?

mu_obs <- ricker(mle["a"], mle["b"], myx$day)
rmse <- function(y, mu) sqrt(mean((y - mu)^2))

set.seed(746)
sims <- replicate(1000, rgamma(nrow(myx), shape = mle["shape"], scale = mu_obs / mle["shape"]))
rmse_sim <- apply(sims, 2, rmse, mu = mu_obs)
rmse_obs <- rmse(myx$titer, mu_obs)
p_boot <- mean(rmse_sim >= rmse_obs)

par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5))
boxplot(c(sims) ~ rep(myx$day, times = 1000), at = sort(unique(myx$day)),
        xlim = c(1.5, 9.5), ylim = c(0, 12), col = "grey90", range = 0,
        xlab = "Days since infection", ylab = "Titer")
points(titer ~ day, data = myx, pch = 16, col = "firebrick")
hist(rmse_sim, breaks = 30, col = "grey80", border = "white", main = "",
     xlab = "RMSE of a simulated data set")
abline(v = rmse_obs, lwd = 3, col = "firebrick")
par(mfrow = c(1, 1))
p_boot   # share of simulated data sets with a larger RMSE

# 3. Posterior predictive checks -----------------------------------------------

library(cmdstanr)
ricker_code <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] day;
  vector<lower=0>[N] titer;
}
parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> shape;
}
model {
  vector[N] mu = a * day .* exp(-b * day);
  titer ~ gamma(shape, shape ./ mu);   // rate = shape / mean
}
"
ricker_model <- cmdstan_model(write_stan_file(ricker_code))
fit_stan <- ricker_model$sample(
  data = list(N = nrow(myx), day = myx$day, titer = myx$titer),
  chains = 4, parallel_chains = 4, seed = 746
)
post <- as.data.frame(posterior::as_draws_df(fit_stan$draws(c("a", "b", "shape"))))

# RMSE of real and simulated titers around each posterior draw's mean
set.seed(746)
use <- sample(nrow(post), 1000)
ppc <- t(sapply(use, function(i) {
  mu <- ricker(post$a[i], post$b[i], myx$day)
  y_rep <- rgamma(nrow(myx), shape = post$shape[i], scale = mu / post$shape[i])
  c(obs = rmse(myx$titer, mu), rep = rmse(y_rep, mu))
}))
p_bayes <- mean(ppc[, "rep"] >= ppc[, "obs"])   # posterior predictive p-value

# posterior predictive 95% band across days: 10 simulated titers per draw and day
band <- sapply(days, function(t) {
  mu <- rep(ricker(post$a, post$b, t), 10)
  shape <- rep(post$shape, 10)
  quantile(rgamma(length(mu), shape = shape, scale = mu / shape), c(0.025, 0.975))
})

par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5))
plot(titer ~ day, data = myx, pch = 16, xlim = c(1, 10), ylim = c(0, 12),
     xlab = "Days since infection", ylab = "Titer")
polygon(c(days, rev(days)), c(band[1, ], rev(band[2, ])), col = adjustcolor("steelblue", 0.25),
        border = NA)
lines(days, ricker(median(post$a), median(post$b), days), lwd = 2)
plot(ppc[, "obs"], ppc[, "rep"], pch = 16, cex = 0.5, col = adjustcolor("black", 0.4),
     xlim = range(ppc), ylim = range(ppc), xlab = "RMSE, real data", ylab = "RMSE, simulated data")
abline(0, 1, col = "firebrick", lwd = 2)
par(mfrow = c(1, 1))
p_bayes

# PREDICTIVE PERFORMANCE =====================================================

# 1. Measures of skill ---------------------------------------------------------

nll_null <- function(p) {
  if (any(p <= 0)) return(Inf)
  -sum(dgamma(myx$titer, shape = p[2], scale = p[1] / p[2], log = TRUE))
}
fit_null <- optim(c(mean = 7, shape = 10), nll_null)

r2 <- 1 - sum((myx$titer - mu_obs)^2) / sum((myx$titer - mean(myx$titer))^2)
mcfadden <- 1 - fit$value / fit_null$value     # NLLs: the minus signs cancel
round(c(R2 = r2, McFadden = mcfadden, RMSE = rmse_obs), 3)

#  [PAUSE] Would an overfit model pass the goodness-of-fit checks? How would
#          its R2 and RMSE look? How well would it predict new rabbits?

# 2. Validation with new data --------------------------------------------------

new1 <- data.frame(day = 2:8, titer = c(4.4, 7.2, 6.8, 5.9, 9.1, 8.3, 8.8))

# R2, McFadden, RMSE for new data, using the parameters fitted to the original rabbits
skill <- function(d) {
  mu <- ricker(mle["a"], mle["b"], d$day)
  ll_model <- sum(dgamma(d$titer, shape = mle["shape"], scale = mu / mle["shape"], log = TRUE))
  ll_null <- sum(dgamma(d$titer, shape = fit_null$par["shape"],
                        scale = fit_null$par["mean"] / fit_null$par["shape"], log = TRUE))
  c(R2 = 1 - sum((d$titer - mu)^2) / sum((d$titer - mean(d$titer))^2),
    McFadden = 1 - ll_model / ll_null,
    RMSE = rmse(d$titer, mu))
}
round(skill(new1), 3)

# 3. Extrapolation -------------------------------------------------------------

new2 <- data.frame(day = 10:16, titer = c(6.8, 8.0, 4.5, 3.1, 2.7, 1.2, 0.04))
round(skill(new2), 3)

# 4. Cross-validation ----------------------------------------------------------

pred_cv <- numeric(nrow(myx))
for (i in seq_len(nrow(myx))) {
  fit_i <- optim(mle, nll_myx, d = myx[-i, ])        # fit without rabbit i
  pred_cv[i] <- ricker(fit_i$par["a"], fit_i$par["b"], myx$day[i])   # predict rabbit i
}

var_explained <- function(pred) 1 - mean((myx$titer - pred)^2) / mean((myx$titer - mean(myx$titer))^2)
round(rbind(all_data = c(RMSE = rmse(myx$titer, mu_obs), var_explained = var_explained(mu_obs)),
            cross_validation = c(rmse(myx$titer, pred_cv), var_explained(pred_cv))), 3)

#  [PAUSE] (widget) Which polynomial degree gives the lowest cross-validation error?

# BINARY RESPONSES: ROC CURVES AND AUC ======================================

# 1. Example: pika occupancy (simulated) -----------------------------------------

set.seed(746)
n <- 120
patches <- data.frame(elev = round(runif(n, 2400, 3400)),
                      area = round(exp(rnorm(n, log(0.8), 0.9)), 2))
for (j in 1:8) patches[[paste0("x", j)]] <- rnorm(n)    # irrelevant covariates
eta <- -0.3 + 1.3 * (patches$elev - 2900) / 300 + 0.8 * (log(patches$area) - log(0.8)) / 0.9
patches$occupied <- rbinom(n, 1, plogis(eta))

f_habitat <- occupied ~ scale(elev) + scale(log(area))
f_big <- update(f_habitat, paste(". ~ . +", paste0("x", 1:8, collapse = " + ")))
m_habitat <- glm(f_habitat, family = binomial, data = patches)
m_big <- glm(f_big, family = binomial, data = patches)

# 10-fold cross-validation: each patch predicted by a model that never saw it
folds <- sample(rep(1:10, length.out = n))
cv_predict <- function(formula) {
  p <- numeric(n)
  for (k in 1:10) {
    m <- glm(formula, family = binomial, data = patches[folds != k, ])
    p[folds == k] <- predict(m, patches[folds == k, ], type = "response")
  }
  p
}
p_cv_habitat <- cv_predict(f_habitat)
p_cv_big <- cv_predict(f_big)

# 2. Thresholds and the ROC curve ----------------------------------------------

#  [PAUSE] (widget) What does the ROC curve look like for a perfect classifier?
#          For a coin flip? For one worse than a coin?

# 3. Area under the curve ------------------------------------------------------

# AUC = P(random occupied patch scores higher than random unoccupied one), from ranks
auc <- function(p, y) {
  r <- rank(p)
  n1 <- sum(y)
  n0 <- sum(1 - y)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ROC curve: true and false positive rates over every threshold
roc_points <- function(p, y) {
  thr <- c(Inf, sort(unique(p), decreasing = TRUE))
  cbind(fpr = sapply(thr, function(t) mean(p[y == 0] >= t)),
        tpr = sapply(thr, function(t) mean(p[y == 1] >= t)))
}

y <- patches$occupied
dev <- function(p) -2 * sum(dbinom(y, 1, p, log = TRUE))
dev_null <- dev(rep(mean(y), n))
perf <- rbind(
  habitat_all_data = c(AUC = auc(fitted(m_habitat), y), McFadden = 1 - dev(fitted(m_habitat)) / dev_null),
  habitat_cv       = c(auc(p_cv_habitat, y), 1 - dev(p_cv_habitat) / dev_null),
  big_all_data     = c(auc(fitted(m_big), y), 1 - dev(fitted(m_big)) / dev_null),
  big_cv           = c(auc(p_cv_big, y), 1 - dev(p_cv_big) / dev_null)
)
round(perf, 3)

plot(roc_points(fitted(m_big), y), type = "s", lwd = 2, col = "grey60",
     xlab = "False positive rate", ylab = "True positive rate")
lines(roc_points(p_cv_big, y), type = "s", lwd = 2, lty = 2, col = "grey60")
lines(roc_points(p_cv_habitat, y), type = "s", lwd = 2, col = "steelblue")
abline(0, 1, lty = 3)
legend("bottomright", c("big model, all data", "big model, cross-validated",
                        "habitat model, cross-validated"),
       lwd = 2, lty = c(1, 2, 1), col = c("grey60", "grey60", "steelblue"), bty = "n")
