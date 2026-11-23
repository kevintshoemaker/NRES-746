#  NRES 746, Week 15: Dynamic models ==========================================================
#  Sections match the Week 15 web page. Needs mvtnorm, RTMB, sdmTMB, unmarked.

# STATE-SPACE MODELS =========================================================================

rm(list = ls())
site <- "https://kevintshoemaker.github.io/NRES-746/data/"

# 1. Two kinds of noise ----------------------------------------------------------------------

# one sage-grouse population, spring male counts 1996-2025 (simulated)
# truth: a = 0.3 * log(300), b = 0.7, c = 0.15, sigma_p = 0.15, sigma_o = 0.20
grouse <- read.csv(paste0(site, "SageGrouseTrend.csv"))
head(grouse, 4)

# Gompertz model on the log scale, x = log(N), w = last spring's precipitation (SD units):
#   state:        x[t] = a + b * x[t-1] + c * w[t] + process noise (SD sigma_p)
#   observation:  y[t] = x[t] + observation noise (SD sigma_o)
plot(grouse$year, grouse$males, type = "b", pch = 16, xlab = "Year", ylab = "Males counted")

#  [PAUSE] A very wet spring adds 0.3 to x. With b = 0.7, how much is left one year later? Two years?
#          What is left one year after an observer overcounts by the same amount?

# 2. Fitting with one kind of noise ----------------------------------------------------------

y <- log(grouse$males)
w <- grouse$precip_z
n <- length(y)

# process noise only: regress this year's log count on last year's (one-step-ahead fitting)
fit_proc <- lm(y[-1] ~ y[-n] + w[-1])
round(summary(fit_proc)$coefficients, 3)

# observation noise only: one deterministic trajectory from x[1] (trajectory matching)
nll_traj <- function(par) {
  x <- numeric(n)
  x[1] <- par[4]
  for (t in 2:n) x[t] <- par[1] + par[2] * x[t - 1] + par[3] * w[t]
  -sum(dnorm(y, x, exp(par[5]), log = TRUE))
}
fit_traj <- optim(c(coef(fit_proc), y[1], log(0.2)), nll_traj, method = "BFGS", control = list(maxit = 500))
round(c(a = fit_traj$par[1], b = fit_traj$par[2], c = fit_traj$par[3], sigma_o = exp(fit_traj$par[5])), 3)

#  [PAUSE] The 2016 count was well below the counts on either side. What does each fit say happened,
#          and what does each predict for 2017 as a result?

# 3. The state-space model -------------------------------------------------------------------

# (a) one multivariate Normal: V[s,t] = sigma_p^2 / (1 - b^2) * b^|s-t| + sigma_o^2 * [s = t]
nll_mvn <- function(par) {
  a <- par[1]
  b <- plogis(par[2])          # keeps 0 < b < 1
  c <- par[3]
  sigma_p <- exp(par[4])
  sigma_o <- exp(par[5])
  mu <- numeric(n)
  m <- a / (1 - b)
  for (t in seq_len(n)) mu[t] <- m <- a + b * m + c * w[t]
  V <- sigma_p^2 / (1 - b^2) * b^abs(outer(1:n, 1:n, "-")) + sigma_o^2 * diag(n)
  -mvtnorm::dmvnorm(y, mu, V, log = TRUE)
}
start <- c(1.7, qlogis(0.7), 0.1, log(0.15), log(0.15))
fit_mvn <- optim(start, nll_mvn, method = "BFGS")
fit_mvn$value

# (b) the Kalman filter: predict, score, update, one year at a time
#  d: covariate effect on the observation (0 until section 4)
kalman <- function(a, b, c, sigma_p, sigma_o, d = 0, y = log(grouse$males), w = grouse$precip_z) {
  n <- length(y)
  x_pred <- P_pred <- x_filt <- P_filt <- numeric(n)
  x_hat <- a / (1 - b)                   # belief before year 1: at equilibrium
  P <- sigma_p^2 / (1 - b^2)
  nll <- 0
  for (t in seq_len(n)) {
    x_hat <- a + b * x_hat + c * w[t]    # predict
    P <- b^2 * P + sigma_p^2
    x_pred[t] <- x_hat
    P_pred[t] <- P
    nll <- nll - dnorm(y[t], x_hat + d * w[t], sqrt(P + sigma_o^2), log = TRUE)    # score
    K <- P / (P + sigma_o^2)             # update
    x_hat <- x_hat + K * (y[t] - d * w[t] - x_hat)
    P <- (1 - K) * P
    x_filt[t] <- x_hat
    P_filt[t] <- P
  }
  list(nll = nll, x_pred = x_pred, P_pred = P_pred, x_filt = x_filt, P_filt = P_filt, K = K)
}
nll_ss <- function(par, ...) kalman(par[1], plogis(par[2]), par[3], exp(par[4]), exp(par[5]), ...)$nll
c(mvn = nll_mvn(start), kalman = nll_ss(start))       # same likelihood, two ways
fit_ss <- optim(start, nll_ss, method = "BFGS", hessian = TRUE)

est_ss <- c(a = fit_ss$par[1], b = plogis(fit_ss$par[2]), c = fit_ss$par[3],
            sigma_p = exp(fit_ss$par[4]), sigma_o = exp(fit_ss$par[5]))
se_ss <- sqrt(diag(solve(fit_ss$hessian)))            # on the scale optim() worked on
round(est_ss, 3)
plogis(fit_ss$par[2] + c(-1.96, 1.96) * se_ss[2])     # 95% interval for b
exp(est_ss["a"] / (1 - est_ss["b"]))                  # equilibrium, males
kf <- kalman(est_ss["a"], est_ss["b"], est_ss["c"], est_ss["sigma_p"], est_ss["sigma_o"])
kf$K                                                  # weight on each new count

#  [PAUSE] In the widget, raise the observation noise from 0.05 to 0.4 (process noise 0.15).
#          What happens to the Kalman gain, and to how the filter treats one unusually high count?

# (c) hidden states as random effects: write the joint likelihood, let RTMB integrate the states out
library(RTMB)
nll_joint <- function(p) {
  b <- plogis(p$logit_b)
  sigma_p <- exp(p$log_sigma_p)
  sigma_o <- exp(p$log_sigma_o)
  x <- p$x                                                                              # the hidden states
  nll <- -dnorm(x[1], p$a / (1 - b) + p$c * w[1], sigma_p / sqrt(1 - b^2), log = TRUE)  # first state
  nll <- nll - sum(dnorm(x[-1], p$a + b * x[-n] + p$c * w[-1], sigma_p, log = TRUE))    # state equation
  nll - sum(dnorm(y, x, sigma_o, log = TRUE))                                           # observation equation
}
par_list <- list(a = 1.7, logit_b = qlogis(0.7), c = 0.1, log_sigma_p = log(0.15), log_sigma_o = log(0.15), x = y)
obj <- MakeADFun(nll_joint, par_list, random = "x", silent = TRUE)
fit_tmb <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
c(kalman = fit_ss$value, rtmb = fit_tmb$objective)

# 4. States, covariates, and forecasts -------------------------------------------------------

# smoothed states: the random-effect estimates, using the whole series
smooth <- as.data.frame(summary(sdr, "random"))
names(smooth) <- c("x_hat", "se")
head(round(smooth, 3), 3)

plot(grouse$year, grouse$males, pch = 16, log = "y", xlab = "Year", ylab = "Males (log scale)")
lines(grouse$year, exp(smooth$x_hat), col = "steelblue", lwd = 2)
lines(grouse$year, exp(smooth$x_hat - 1.96 * smooth$se), col = "steelblue", lty = 2)
lines(grouse$year, exp(smooth$x_hat + 1.96 * smooth$se), col = "steelblue", lty = 2)
lines(grouse$year, exp(grouse$log_true), lty = 3)     # the truth (simulated data only)

# the covariate in the observation equation instead: y[t] = x[t] + d * w[t] + noise
nll_ss_obs <- function(par) kalman(par[1], plogis(par[2]), 0, exp(par[4]), exp(par[5]), d = par[3])$nll
fit_ss_obs <- optim(start, nll_ss_obs, method = "BFGS")
c(precip_in_state = -fit_ss$value, precip_in_observation = -fit_ss_obs$value)      # log-likelihoods
est_ss["c"] * est_ss["b"]^(0:5)                       # effect of one wet spring, 0 to 5 years later
est_ss["c"] / (1 - est_ss["b"])                       # long-run effect of a lasting 1-SD change

#  [PAUSE] Which equation does each belong in: (a) acres burned last summer; (b) number of mornings
#          each lek was visited; (c) snow depth on the day of the count?

# forecast: the predict step repeated, with no update
forecast_ss <- function(w_future) {
  x_hat <- kf$x_filt[n]
  P <- kf$P_filt[n]
  out <- data.frame(year = max(grouse$year) + seq_along(w_future), x_hat = NA, sd = NA)
  for (h in seq_along(w_future)) {
    x_hat <- est_ss["a"] + est_ss["b"] * x_hat + est_ss["c"] * w_future[h]
    P <- est_ss["b"]^2 * P + est_ss["sigma_p"]^2
    out$x_hat[h] <- x_hat
    out$sd[h] <- sqrt(P)
  }
  out
}
fc_avg <- forecast_ss(rep(0, 5))       # average precipitation
fc_dry <- forecast_ss(rep(-1, 5))      # five dry springs in a row
round(cbind(year = fc_avg$year, average = exp(fc_avg$x_hat), dry = exp(fc_dry$x_hat)))
round(exp(cbind(lower = fc_avg$x_hat - 1.96 * fc_avg$sd, upper = fc_avg$x_hat + 1.96 * fc_avg$sd)))

# 5. What goes wrong -------------------------------------------------------------------------

# 500 series from the true model, fit both ways (takes a minute or two)
sim_gompertz <- function(n, a = 0.3 * log(300), b = 0.7, c = 0.15, sigma_p = 0.15, sigma_o = 0.20) {
  w <- rnorm(n)
  x <- numeric(n)
  x_prev <- rnorm(1, a / (1 - b), sigma_p / sqrt(1 - b^2))
  for (t in seq_len(n)) x_prev <- x[t] <- a + b * x_prev + c * w[t] + rnorm(1, 0, sigma_p)
  list(y = x + rnorm(n, 0, sigma_o), w = w)
}
fit_both <- function(n) {
  s <- sim_gompertz(n)
  b_lm <- coef(lm(s$y[-1] ~ s$y[-n] + s$w[-1]))[2]
  f <- try(optim(c(1.7, qlogis(0.7), 0.15, log(0.15), log(0.20)), nll_ss, y = s$y, w = s$w, method = "BFGS"),
           silent = TRUE)
  if (inherits(f, "try-error")) return(c(b_lm, NA, NA, NA))
  c(b_lm, plogis(f$par[2]), exp(f$par[4]), exp(f$par[5]))
}
set.seed(746)
sims30 <- t(replicate(500, fit_both(30)))
sims100 <- t(replicate(500, fit_both(100)))
rbind(n30 = colMeans(sims30, na.rm = TRUE), n100 = colMeans(sims100, na.rm = TRUE))   # b (lm), b, sigma_p, sigma_o
c(n30 = mean(sims30[, 3] < 0.01 | sims30[, 4] < 0.01, na.rm = TRUE),                  # share with an SD at zero
  n100 = mean(sims100[, 3] < 0.01 | sims100[, 4] < 0.01, na.rm = TRUE))

# profile log-likelihood for the two SDs (each point maximized over a, b, c)
sd_grid <- seq(0.01, 0.36, length.out = 30)
profile_nll <- outer(sd_grid, sd_grid, Vectorize(function(sp, so) {
  f <- function(p) kalman(p[1], plogis(p[2]), p[3], sp, so)$nll
  optim(fit_ss$par[1:3], f, method = "BFGS")$value
}))
contour(sd_grid, sd_grid, profile_nll - fit_ss$value, levels = c(0.25, 0.5, 1, 2, 3, 5, 8),
        xlab = "Process SD", ylab = "Observation SD")
points(est_ss["sigma_p"], est_ss["sigma_o"], pch = 16)

# if replicate counts had supplied sigma_o = 0.20
nll_known <- function(p) kalman(p[1], plogis(p[2]), p[3], exp(p[4]), 0.20)$nll
fit_known <- optim(fit_ss$par[1:4], nll_known, method = "BFGS", hessian = TRUE)
plogis(fit_known$par[2] + c(-1.96, 1.96) * sqrt(diag(solve(fit_known$hessian)))[2])   # 95% interval for b

#  [PAUSE] Three fits give b = 0.60, 0.91, and 0.75 (95% interval 0.35 to 0.94). The agency asks how
#          fast the population recovers from a bad year. What do you tell them?

# Bonus: reading ACF and PACF plots ----------------------------------------------------------

set.seed(746)
ar1_series <- arima.sim(list(ar = 0.7), n = 500)
ma1_series <- arima.sim(list(ma = 0.7), n = 500)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 0.5))
acf(ar1_series, lag.max = 12, main = "AR(1): ACF decays gradually")
pacf(ar1_series, lag.max = 12, main = "AR(1): PACF stops after lag 1")
acf(ma1_series, lag.max = 12, main = "MA(1): ACF stops after lag 1")
pacf(ma1_series, lag.max = 12, main = "MA(1): PACF decays gradually")
par(mfrow = c(1, 2))
acf(y, lag.max = 10, main = "Grouse log counts: ACF")
pacf(y, lag.max = 10, main = "Grouse log counts: PACF")
par(mfrow = c(1, 1))

# Bonus: counts instead of log counts --------------------------------------------------------

# only the observation line changes: negative binomial with mean exp(x)
males <- grouse$males
nll_joint_nb <- function(p) {
  b <- plogis(p$logit_b)
  sigma_p <- exp(p$log_sigma_p)
  x <- p$x
  mu <- exp(x)
  nll <- -dnorm(x[1], p$a / (1 - b) + p$c * w[1], sigma_p / sqrt(1 - b^2), log = TRUE)
  nll <- nll - sum(dnorm(x[-1], p$a + b * x[-n] + p$c * w[-1], sigma_p, log = TRUE))
  nll - sum(dnbinom2(males, mu, mu + mu^2 / exp(p$log_size), log = TRUE))      # variance = mu + mu^2 / size
}
obj_nb <- MakeADFun(nll_joint_nb, list(a = 1.7, logit_b = qlogis(0.7), c = 0.1, log_sigma_p = log(0.15),
                                       log_size = log(20), x = y), random = "x", silent = TRUE)
fit_nb <- nlminb(obj_nb$par, obj_nb$fn, obj_nb$gr)
round(c(a = unname(fit_nb$par["a"]), b = plogis(unname(fit_nb$par["logit_b"])), c = unname(fit_nb$par["c"]),
        sigma_p = exp(unname(fit_nb$par["log_sigma_p"])), size = exp(unname(fit_nb$par["log_size"]))), 3)

# A SPATIOTEMPORAL ANALYSIS FROM START TO FINISH =============================================

# 1. Problem and data ------------------------------------------------------------------------

# 60 leks counted each spring, 2011-2025 (simulated)
# truth: sagebrush 0.4, precipitation 0.25, sigma_s 0.5, sigma_st 0.4, rho 0.7, size 8
lk <- read.csv(paste0(site, "LekCounts.csv"))
lk$sage_c <- (lk$sagebrush_pct - 20) / 10        # sagebrush cover: 10-point units above 20%
head(lk, 3)
tapply(lk$males, lk$year, sum)                   # district total by year
round(tapply(lk$precip_z, lk$year, mean), 2)     # mean precipitation by year

# 2. Candidate models ------------------------------------------------------------------------

# log(mu[s,t]) = b0 + b1 * sage[s] + b2 * precip[s,t] + omega[s] + eps[s,t],  males ~ NegBin(mu, size)
# M0 no fields; M1 omega only; M2 + independent yearly fields; M3 + AR(1) fields; M4 + random walk

#  [PAUSE] Picture five years of maps of eps under M2 and under M3. How do they differ?
#          Which model has a better chance of forecasting next year?

# 3. Fit and check ---------------------------------------------------------------------------

library(sdmTMB)
#  one call for every candidate with fields; `...` passes extra arguments to sdmTMB()
fit_st <- function(spatiotemporal, data = lk, ...) {
  sdmTMB(males ~ sage_c + precip_z, data = data, mesh = make_mesh(data, c("x_km", "y_km"), cutoff = 5),
         family = nbinom2(), spatial = "on", time = "year", spatiotemporal = spatiotemporal, ...)
}
m0 <- sdmTMB(males ~ sage_c + precip_z, data = lk, family = nbinom2(), spatial = "off")
m1 <- fit_st("off")
m2 <- fit_st("iid")
m3 <- fit_st("ar1")
m4 <- fit_st("rw")

sanity(m3)
tidy(m3, conf.int = TRUE)
tidy(m3, "ran_pars", conf.int = TRUE)       # sigma_O = sigma_s, sigma_E = sigma_st, rho = AR(1), phi = NB size
tidy(m0)                                    # compare the sagebrush SE

# does a lek's residual this year predict its residual next year?
set.seed(746)
resid_lag1 <- function(fit) {
  r <- residuals(fit, type = "mle-mvn")       # randomized quantile residuals, Normal if the model is right
  o <- order(lk$lek, lk$year)
  same_lek <- lk$lek[o][-1] == lk$lek[o][-nrow(lk)]
  cor(r[o][-1][same_lek], r[o][-nrow(lk)][same_lek])
}
lag1 <- c(M1 = resid_lag1(m1), M2 = resid_lag1(m2), M3 = resid_lag1(m3))
round(lag1, 2)                                # lag-1 correlation of residuals within leks

# 4. Compare ---------------------------------------------------------------------------------

aic <- AIC(m0, m1, m2, m3, m4)
aic$dAIC <- aic$AIC - min(aic$AIC)
round(aic, 1)

# one-year-ahead forecasts from seven origins (takes a few minutes)
forecast_one <- function(origin) {
  train <- lk[lk$year <= origin, ]
  test <- lk[lk$year == origin + 1, ]
  fits <- list(M0 = sdmTMB(males ~ sage_c + precip_z, data = train, family = nbinom2(), spatial = "off"),
               M1 = fit_st("off", data = train, extra_time = origin + 1),
               M2 = fit_st("iid", data = train, extra_time = origin + 1),
               M3 = fit_st("ar1", data = train, extra_time = origin + 1),
               M4 = fit_st("rw", data = train, extra_time = origin + 1))
  sapply(fits, function(f) mean(abs(test$males - exp(predict(f, newdata = test)$est))))
}
forecast_err <- t(sapply(2018:2024, forecast_one))
rownames(forecast_err) <- 2019:2025
round(rbind(forecast_err, mean = colMeans(forecast_err)), 1)      # mean absolute error, males per lek

#  [PAUSE] M3 says a local departure fades (rho < 1); M4 says it persists (rho = 1). They forecast one
#          year about equally well. Would that hold for ten years? What does each say about a lek that
#          is far below expectation today?

# 5. Fields and forecasts --------------------------------------------------------------------

m3f <- fit_st("ar1", extra_time = 2026:2027)
grid <- expand.grid(x_km = seq(1, 99, 2), y_km = seq(1, 79, 2), year = c(2019, 2022, 2025, 2027),
                    sage_c = 0, precip_z = 0)
pred_grid <- predict(m3f, newdata = grid)
gx <- sort(unique(grid$x_km))
gy <- sort(unique(grid$y_km))
zl <- c(-1, 1) * max(abs(pred_grid$epsilon_st))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
for (yr in c(2019, 2022, 2025, 2027)) {
  image(gx, gy, matrix(pred_grid$epsilon_st[pred_grid$year == yr], length(gx)), zlim = zl, asp = 1,
        col = hcl.colors(65, "Blue-Red 3"), xlab = "Easting (km)", ylab = "Northing (km)", main = yr)
}
par(mfrow = c(1, 1))

# how much of the drop since the peak year is precipitation, and how much is the field?
yr_peak <- as.numeric(names(which.max(tapply(lk$males, lk$year, sum))))
p_obs <- predict(m3)
p_obs$precip_part <- tidy(m3)$estimate[3] * p_obs$precip_z
parts <- aggregate(cbind(precip_part, epsilon_st) ~ year, data = p_obs, FUN = mean)
round(parts[parts$year %in% c(yr_peak, 2025), ], 2)

# forecast the district total under two precipitation scenarios
leks <- unique(lk[, c("lek", "x_km", "y_km", "sage_c")])
new_avg <- merge(leks, data.frame(year = 2026:2027, precip_z = 0))      # every lek in both years
new_dry <- transform(new_avg, precip_z = -1)
set.seed(746)
sim_avg <- predict(m3f, newdata = new_avg, nsim = 500)                  # rows = lek-years, columns = draws
sim_dry <- predict(m3f, newdata = new_dry, nsim = 500)
total_by_year <- function(sims) {
  sapply(2026:2027, function(yr) quantile(colSums(exp(sims[new_avg$year == yr, ])), c(0.025, 0.5, 0.975)))
}
fc_tab <- rbind(average = as.vector(total_by_year(sim_avg)), dry = as.vector(total_by_year(sim_dry)))
colnames(fc_tab) <- c("2026 low", "2026 median", "2026 high", "2027 low", "2027 median", "2027 high")
round(fc_tab)                                                           # expected district total, males

#  [PAUSE] The 2027 intervals are wider than 2026's. Name two separate reasons.

# DYNAMIC OCCUPANCY ==========================================================================

# 1. A hidden state that is 0 or 1 -----------------------------------------------------------

# 100 talus patches, 8 summers, 3 visits per summer (simulated)
# truth: psi logit 0.5 + 1.0 * elev; extinction logit -1.6 + 0.8 * temp; colonization 0.12; detection 0.6
pikas <- read.csv(paste0(site, "PikaPatches.csv"))
head(pikas, 3)

years <- sort(unique(pikas$year))
n_yr <- length(years)
n_patch <- length(unique(pikas$patch))
pikas <- pikas[order(pikas$year, pikas$patch), ]
det <- matrix(pikas$visit1 + pikas$visit2 + pikas$visit3, n_patch, n_yr)     # detections (0 to 3): patch by year
elev_c <- (pikas$elev_m[pikas$year == years[1]] - 2800) / 500                # 500-m units above 2800 m
temp <- pikas$temp_z[match(years, pikas$year)]
rbind(temp = temp, detected = colMeans(det > 0))                             # share of patches with a detection

# 2. The likelihood by matrix multiplication -------------------------------------------------

# forward algorithm: alpha = (P(occupied and data so far), P(empty and data so far))
patch_lik <- function(d, psi, gam, eps, p) {
  alpha <- c(psi, 1 - psi)
  for (t in seq_along(d)) {
    if (t > 1) alpha <- alpha %*% matrix(c(1 - eps[t - 1], eps[t - 1], gam, 1 - gam), 2, byrow = TRUE)   # predict
    alpha <- alpha * c(p^d[t] * (1 - p)^(3 - d[t]), d[t] == 0)                                        # update
  }
  sum(alpha)
}
patch_lik(c(2, 0, 1), psi = 0.5, gam = 0.1, eps = c(0.2, 0.2), p = 0.6)      # the three-year example

#  [PAUSE] After a year with no detections, alpha is (0.00369, 0.0144). What is the probability the
#          patch was occupied that year, given the data so far?

nll_occ <- function(par) {
  psi <- plogis(par[1] + par[2] * elev_c)
  gam <- plogis(par[3])
  eps <- plogis(par[4] + par[5] * temp)        # temp[t] acts on the transition from year t to t + 1
  p <- plogis(par[6])
  -sum(sapply(seq_len(n_patch), function(i) log(patch_lik(det[i, ], psi[i], gam, eps, p))))
}
fit_occ <- optim(c(0, 0, -2, -1, 0, 0), nll_occ, method = "BFGS", hessian = TRUE)
round(cbind(estimate = fit_occ$par, se = sqrt(diag(solve(fit_occ$hessian)))), 2)

# the same model in unmarked
library(unmarked)
y_wide <- matrix(NA, n_patch, 3 * n_yr)
for (t in seq_len(n_yr)) {
  y_wide[, 3 * (t - 1) + 1:3] <- as.matrix(pikas[pikas$year == years[t], c("visit1", "visit2", "visit3")])
}
umf <- unmarkedMultFrame(y = y_wide, siteCovs = data.frame(elev_c = elev_c),
                         yearlySiteCovs = list(temp = matrix(temp, n_patch, n_yr, byrow = TRUE)),
                         numPrimary = n_yr)
fit_colext <- colext(psiformula = ~ elev_c, gammaformula = ~ 1, epsilonformula = ~ temp, pformula = ~ 1,
                     data = umf)
c(by_hand = -fit_occ$value, colext = as.numeric(logLik(fit_colext)))        # log-likelihoods

# extinction probability after a cool, average, and hot summer; detection; equilibrium occupancy
eps_hat <- function(tz) plogis(fit_occ$par[4] + fit_occ$par[5] * tz)
gam_hat <- plogis(fit_occ$par[3])
p_hat <- plogis(fit_occ$par[6])
round(c(cool = eps_hat(-1), average = eps_hat(0), hot = eps_hat(1), p = p_hat, miss3 = (1 - p_hat)^3), 2)
round(c(average = gam_hat / (gam_hat + eps_hat(0)), hot = gam_hat / (gam_hat + eps_hat(1))), 2)

# occupancy trajectory: psi[t+1] = psi[t] * (1 - eps[t]) + (1 - psi[t]) * gamma
psi_hat <- matrix(NA, n_patch, n_yr)
psi_hat[, 1] <- plogis(fit_occ$par[1] + fit_occ$par[2] * elev_c)
for (t in 2:n_yr) psi_hat[, t] <- psi_hat[, t - 1] * (1 - eps_hat(temp[t - 1])) + (1 - psi_hat[, t - 1]) * gam_hat
plot(years, colMeans(psi_hat), type = "b", pch = 16, ylim = c(0.2, 0.8), xlab = "Year",
     ylab = "Share of patches occupied")
lines(years, colMeans(det > 0), type = "b", lty = 2)                        # naive
lines(years, tapply(pikas$occupied_true, pikas$year, mean), type = "b", pch = 4, lty = 3)   # truth

#  [PAUSE] A patch had detections in 2019 and 2021 and none in three visits in 2020. Missed three times,
#          or extinct and recolonized? Which is more probable, and by how much?
