#  NRES 746, Week 14: Spatial regression =====================================================
#  Sections match the Week 14 web page. Needs mvtnorm, gstat, sp, glmmTMB, sdmTMB, fmesher, INLA, inlabru, sf, mgcv, DHARMa.

# SPATIAL REGRESSION =========================================================================

rm(list = ls())
library(mvtnorm)

# 1. Nearby plots are alike ------------------------------------------------------------------

# 200 cheatgrass plots in a 40 x 30 km area (simulated); coordinates in km
# truth: fire effect 0.5, elevation -2 per km, sigma_s 0.7, ell 6 km, sigma_e 0.3
site <- "https://kevintshoemaker.github.io/NRES-746/data/"
cg <- read.csv(paste0(site, "CheatgrassSurvey.csv"))
grid <- read.csv(paste0(site, "CheatgrassGrid.csv"))        # 0.5-km prediction grid
cg$elev_km <- (cg$elev_m - 1800) / 1000                     # km above 1800 m
grid$elev_km <- (grid$elev_m - 1800) / 1000
head(cg)

fit_lm <- lm(log_biomass ~ burned + elev_km, data = cg)
round(summary(fit_lm)$coefficients, 3)

# map of the residuals
cg$resid_lm <- resid(fit_lm)
plot(cg$x_km, cg$y_km, asp = 1, pch = 21, cex = 0.5 + 1.6 * abs(cg$resid_lm),
     bg = ifelse(cg$resid_lm > 0, "firebrick", "steelblue"), xlab = "Easting (km)", ylab = "Northing (km)")

#  [PAUSE] lm() assumes independent residuals. What on this map says otherwise?

# 2000 new surveys from the true model: how variable is lm's fire estimate, and what SE does it report?
X <- model.matrix(~ burned + elev_km, data = cg)
D <- as.matrix(dist(cg[, c("x_km", "y_km")]))            # distances between plots (km)
V_true <- 0.7^2 * exp(-D / 6) + 0.3^2 * diag(nrow(cg))
set.seed(746)
Y <- rmvnorm(2000, X %*% c(3, 0.5, -2), V_true)
sims <- t(apply(Y, 1, function(y) summary(lm(y ~ X - 1))$coefficients[2, 1:2]))
c(sd_of_estimates = sd(sims[, 1]), mean_reported_se = mean(sims[, 2]),
  coverage_95 = mean(abs(sims[, 1] - 0.5) < qt(0.975, 197) * sims[, 2]))

# 2. The variogram ---------------------------------------------------------------------------

# semivariance: half the mean squared difference between residuals of plots about d apart
library(gstat)
library(sp)
cg_sp <- cg
coordinates(cg_sp) <- ~ x_km + y_km
vg_lm <- variogram(resid_lm ~ 1, cg_sp, cutoff = 20, width = 2)     # 2-km distance bins
vg_lm[, c("np", "dist", "gamma")]
plot(vg_lm$dist, vg_lm$gamma, pch = 16, xlim = c(0, 20), ylim = c(0, 0.7),
     xlab = "Distance between plots (km)", ylab = "Semivariance of residuals")

#  [PAUSE] Where are the nugget, the sill, and the range on this plot?

# 3. A Gaussian process over the map ---------------------------------------------------------

# y ~ MVN(X beta, V), V = sigma_s^2 exp(-D / ell) + sigma_e^2 I
y <- cg$log_biomass
nll_spatial <- function(par) {
  beta <- par[1:3]
  sigma_s <- exp(par[4])
  ell <- exp(par[5])
  sigma_e <- exp(par[6])
  V <- sigma_s^2 * exp(-D / ell) + sigma_e^2 * diag(length(y))
  -dmvnorm(y, as.vector(X %*% beta), V, log = TRUE)
}
fit_hand <- optim(c(coef(fit_lm), log(0.5), log(5), log(0.5)), nll_spatial, method = "BFGS", hessian = TRUE)
se_hand <- sqrt(diag(solve(fit_hand$hessian)))
est_hand <- c(fit_hand$par[1:3], exp(fit_hand$par[4:6]))
names(est_hand) <- c("intercept", "burned", "elev_km", "sigma_s", "ell", "sigma_e")
round(est_hand, 3)
round(rbind(lm = c(coef(fit_lm)[2:3], sqrt(diag(vcov(fit_lm)))[2:3]),
            spatial = c(fit_hand$par[2:3], se_hand[2:3])), 3)        # estimates, then SEs
c(logLik_lm = as.numeric(logLik(fit_lm)), logLik_spatial = -fit_hand$value)

# the fitted model's variogram on top of the lm residual variogram
curve(est_hand["sigma_e"]^2 + est_hand["sigma_s"]^2 * (1 - exp(-x / est_hand["ell"])), add = TRUE, lwd = 2)

#  [PAUSE] Why did the fire effect change, and why did the elevation SE triple?

# 4. A mesh makes it fast --------------------------------------------------------------------

library(sdmTMB)
mesh <- make_mesh(cg, xy_cols = c("x_km", "y_km"), cutoff = 2)     # no two vertices closer than 2 km
mesh$mesh$n                                                        # number of vertices
plot(mesh)

# each plot's field value is a weighted average of the three vertices of its triangle
A <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(cg[, c("x_km", "y_km")]))
dim(A)
i <- which(Matrix::rowSums(A != 0) == 3)[1]
A[i, A[i, ] != 0]

# 5. One model, two packages -----------------------------------------------------------------

# glmmTMB: dense covariance matrix, exponential correlation (the model fit by hand)
library(glmmTMB)
cg$pos <- numFactor(cg$x_km, cg$y_km)
cg$all <- factor(1)
fit_glmmtmb <- glmmTMB(log_biomass ~ burned + elev_km + exp(pos + 0 | all), data = cg)
c(logLik(fit_glmmtmb), by_hand = -fit_hand$value)
summary(fit_glmmtmb)$coefficients$cond
exp(getME(fit_glmmtmb, "theta"))          # sigma_s, ell
sigma(fit_glmmtmb)                        # sigma_e

# sdmTMB: mesh, Matern (smoothness 1), maximum likelihood
fit_sdmtmb <- sdmTMB(log_biomass ~ burned + elev_km, data = cg, mesh = mesh, spatial = "on")
fit_sdmtmb
tidy(fit_sdmtmb, "ran_pars")              # range, phi = sigma_e, sigma_O = sigma_s

# 6. The GAM alternative: a smooth of the coordinates ----------------------------------------

library(mgcv)
fit_gam <- gam(log_biomass ~ burned + elev_km + s(x_km, y_km, k = 100), data = cg, method = "REML")
summary(fit_gam)$p.table
summary(fit_gam)$s.table

fit_gam_gp <- gam(log_biomass ~ burned + elev_km + s(x_km, y_km, bs = "gp", m = 2, k = 100),
                  data = cg, method = "REML")       # m = 2: exponential correlation
summary(fit_gam_gp)$p.table

# elevation is close to a plane in the coordinates, and both smooths leave planes unpenalized
summary(lm(elev_km ~ x_km + y_km, data = cg))$r.squared

# the fitted surfaces against the true field
pred_sd <- predict(fit_sdmtmb, newdata = grid)
gam_surface <- predict(fit_gam, newdata = grid, type = "terms")[, "s(x_km,y_km)"]
c(sdmTMB = cor(pred_sd$omega_s, grid$field_true), gam = cor(gam_surface, grid$field_true))

#  [PAUSE] Elevation SE: 0.19 (lm), 0.59 (exponential field), 3.9 (thin plate smooth). Which would you report?

# Bonus: a Bayesian fit with R-INLA ---------------------------------------------------------

# same mesh as sdmTMB, through inlabru; penalized complexity priors
library(INLA)
library(inlabru)
library(sf)
cg_sf <- st_as_sf(cg, coords = c("x_km", "y_km"), remove = FALSE)   # plots as spatial points
spde <- inla.spde2.pcmatern(mesh$mesh,
                            prior.range = c(5, 0.5),     # P(range < 5 km) = 0.5
                            prior.sigma = c(1, 0.05))    # P(sigma_s > 1) = 0.05
fit_inla <- bru(log_biomass ~ Intercept(1) + burned + elev_km + field(geometry, model = spde),
                data = cg_sf, family = "gaussian")
round(fit_inla$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")], 3)
round(fit_inla$summary.hyperpar[, c("mean", "0.025quant", "0.5quant", "0.975quant")], 3)

# Bonus: the same model in gls() --------------------------------------------------------------

library(nlme)
fit_gls <- gls(log_biomass ~ burned + elev_km, data = cg, method = "ML",
               correlation = corExp(form = ~ x_km + y_km, nugget = TRUE))
round(summary(fit_gls)$tTable, 3)
c(logLik_gls = as.numeric(logLik(fit_gls)), logLik_by_hand = -fit_hand$value)

# gls's range = ell; nugget = sigma_e^2 / (sigma_s^2 + sigma_e^2); residual SE = sqrt(sigma_s^2 + sigma_e^2)
gls_cor <- coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)
total_sd <- sigma(fit_gls)
round(c(sigma_s = total_sd * sqrt(1 - gls_cor[["nugget"]]), ell = gls_cor[["range"]],
        sigma_e = total_sd * sqrt(gls_cor[["nugget"]])), 3)
intervals(fit_gls, which = "var-cov")

# A SPATIAL ANALYSIS FROM START TO FINISH ====================================================

# (run section 1 and the mesh in section 4 first)

# 2. The mesh --------------------------------------------------------------------------------

mesh_check <- t(sapply(c(1, 2, 3, 5, 8), function(co) {
  m <- make_mesh(cg, c("x_km", "y_km"), cutoff = co)
  f <- sdmTMB(log_biomass ~ burned + elev_km, data = cg, mesh = m, spatial = "on")
  rp <- tidy(f, "ran_pars")
  c(cutoff = co, vertices = m$mesh$n, fire = tidy(f)$estimate[2], fire_se = tidy(f)$std.error[2],
    range = rp$estimate[rp$term == "range"], sigma_s = rp$estimate[rp$term == "sigma_O"],
    logLik = as.numeric(logLik(f)))
}))
round(mesh_check, 2)

# 3. Fit and check ---------------------------------------------------------------------------

m0 <- sdmTMB(log_biomass ~ burned + elev_km, data = cg, spatial = "off")
m1 <- sdmTMB(log_biomass ~ burned + elev_km, data = cg, mesh = mesh, spatial = "on")
m2 <- sdmTMB(log_biomass ~ burned * elev_km, data = cg, mesh = mesh, spatial = "on")
sanity(m1)

# DHARMa: where does each observation fall among simulations from the fitted model?
library(DHARMa)
set.seed(746)
dh_m1 <- dharma_residuals(simulate(m1, nsim = 250, type = "mle-mvn"), m1, return_DHARMa = TRUE, plot = FALSE)
plot(dh_m1)                               # uniform QQ plot (KS, dispersion, outlier tests); residuals vs. predicted
testDispersion(dh_m1, plot = FALSE)       # flags underdispersion: an artifact of keeping the fitted field

# spatial autocorrelation left in the residuals: Moran's I, M0 vs. M1
dh_m0 <- dharma_residuals(simulate(m0, nsim = 250, type = "mle-mvn"), m0, return_DHARMa = TRUE, plot = FALSE)
testSpatialAutocorrelation(dh_m0, x = cg$x_km, y = cg$y_km, plot = FALSE)
testSpatialAutocorrelation(dh_m1, x = cg$x_km, y = cg$y_km, plot = FALSE)

# residual variograms: flat near 1 if no spatial pattern is left
cg$resid_m1 <- residuals(m1, type = "mle-mvn")
cg_sp$resid_m1 <- cg$resid_m1
cg_sp$resid_m0 <- as.numeric(scale(cg$resid_lm))
vg_m0 <- variogram(resid_m0 ~ 1, cg_sp, cutoff = 20, width = 2)
vg_m1 <- variogram(resid_m1 ~ 1, cg_sp, cutoff = 20, width = 2)
plot(vg_m0$dist, vg_m0$gamma, ylim = c(0, 1.5), xlab = "Distance between plots (km)",
     ylab = "Semivariance of standardized residuals")
points(vg_m1$dist, vg_m1$gamma, pch = 16)

# 4. Compare ---------------------------------------------------------------------------------

aic <- AIC(m0, m1, m2)
aic$dAIC <- aic$AIC - min(aic$AIC)
aic
tidy(m2)

# cross-validation: random folds vs. whole 10-km blocks held out
set.seed(746)
fold_random <- sample(rep(1:10, length.out = nrow(cg)))
fold_block <- as.integer(factor(paste(floor(cg$x_km / 10), floor(cg$y_km / 10))))
cv_rmse <- function(formula, spatial, folds) {
  cv <- sdmTMB_cv(formula, data = cg, mesh = mesh, spatial = spatial, fold_ids = folds)
  sqrt(mean((cv$data$cv_predicted - cv$data$log_biomass)^2))
}
cv_tab <- rbind(
  random_folds = c(M0 = cv_rmse(log_biomass ~ burned + elev_km, "off", fold_random),
                   M1 = cv_rmse(log_biomass ~ burned + elev_km, "on", fold_random)),
  block_folds = c(M0 = cv_rmse(log_biomass ~ burned + elev_km, "off", fold_block),
                  M1 = cv_rmse(log_biomass ~ burned + elev_km, "on", fold_block)))
round(cv_tab, 3)                  # root mean squared prediction error, log biomass

#  [PAUSE] Which cross-validation result describes how well a map of the surveyed area will do?

# 5. Predict ---------------------------------------------------------------------------------

pred <- predict(m1, newdata = grid)
head(pred[, c("x_km", "y_km", "est", "est_non_rf", "omega_s")], 3)
set.seed(746)
sims <- predict(m1, newdata = grid, nsim = 500)      # 500 draws of the predicted surface
pred$se <- apply(sims, 1, sd)
pred$p_exceed <- rowMeans(sims > log(60))            # P(expected log biomass is above the threshold)

# maps: prediction, its standard error, probability of exceeding 60 g per square meter
ux <- sort(unique(grid$x_km))
uy <- sort(unique(grid$y_km))
show_map <- function(v, main, pal) {
  image(ux, uy, matrix(v, length(ux)), asp = 1, col = hcl.colors(64, pal, rev = TRUE), main = main,
        xlab = "Easting (km)", ylab = "Northing (km)")
  points(cg$x_km, cg$y_km, pch = ".")
}
par(mfrow = c(2, 2))
show_map(pred$est_non_rf, "Regression part", "YlOrRd")
show_map(pred$est, "Prediction: regression + field", "YlOrRd")
show_map(pred$se, "Standard error", "Purples")
show_map(pred$p_exceed, "P(biomass above 60 g per square meter)", "Reds")
par(mfrow = c(1, 1))

# against the true surface (known because the data are simulated)
pred0 <- predict(m0, newdata = grid)
c(rmse_m0 = sqrt(mean((pred0$est - grid$log_biomass_true)^2)),
  rmse_m1 = sqrt(mean((pred$est - grid$log_biomass_true)^2)))

# 6. What we tell the district ---------------------------------------------------------------

tidy(m1, conf.int = TRUE)
tidy(m1, "ran_pars", conf.int = TRUE)
exp(unlist(tidy(m1, conf.int = TRUE)[2, c("estimate", "conf.low", "conf.high")]))      # fire, as a multiplier

# Bonus: counts and other responses ----------------------------------------------------------

fit_pois_glm <- glm(seedlings ~ burned + elev_km, data = cg, family = poisson)
fit_pois_sp <- sdmTMB(seedlings ~ burned + elev_km, data = cg, mesh = mesh, spatial = "on", family = poisson())
round(rbind(glm = summary(fit_pois_glm)$coefficients[2, 1:2],
            spatial = unlist(tidy(fit_pois_sp)[2, c("estimate", "std.error")])), 3)     # fire effect
