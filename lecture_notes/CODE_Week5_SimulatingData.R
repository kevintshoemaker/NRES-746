#  NRES 746, Week 5: Simulating data from a model ------------------------------
#  Code from Week5_SimulatingData.html. Base R only.

# SIMULATING DATA FROM A MODEL ==================================================

# 2. WOD: choosing a data-generating model -------------------------------------
#  lek counts: mu_t = a exp(-b t), Y_t ~ NegBin(mu_t, k), with a = 55, b = 0.06, k = 8

set.seed(746)
yr <- 0:19
mu_lek <- 55 * exp(-0.06 * yr)                    # deterministic: expected count
lek <- rnbinom(20, mu = mu_lek, size = 8)         # stochastic: noise around it
lek
#  [PAUSE] data sets 2 (nest fate vs. cover) and 3 (chick mass vs. forbs): which function, which noise?

# 3. Random number generators --------------------------------------------------

set.seed(3)
rnorm(5, mean = 10, sd = 4.1)
rpois(5, lambda = 5.4)
round(rbeta(5, shape1 = 0.1, shape2 = 0.1), 3)

par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.5))
hist(rnorm(1000, 10, 4.1), breaks = 30, col = "grey85", border = "white",
     main = "Normal(10, 4.1)", xlab = "Value")
barplot(table(rpois(1000, 5.4)), col = "grey60", border = NA,
        main = "Poisson(5.4)", xlab = "Count", ylab = "Frequency")
hist(rbeta(1000, 0.1, 0.1), breaks = 30, col = "grey85", border = "white",
     main = "Beta(0.1, 0.1)", xlab = "Value")
par(mfrow = c(1, 1))

# 4. Signal and noise -----------------------------------------------------------
#  linear: mu = a + b x, Y ~ Normal(mu, sigma)

x_lin <- 1:20
a <- 2
b <- 1
y_det <- a + b * x_lin                            # deterministic part
set.seed(1001)
y_lin <- rnorm(20, mean = y_det, sd = 2)          # stochastic part

#  hyperbolic: mu = a b / (b + x), Y ~ NegBin(mu, k); variance mu + mu^2 / k

a <- 20
b <- 1
k <- 5
set.seed(1001)
x_hyp <- runif(50, min = 0, max = 5)              # densities sampled from the field
mu_hyp <- a * b / (b + x_hyp)
y_hyp <- rnbinom(50, mu = mu_hyp, size = k)

par(mfrow = c(1, 2), mar = c(4, 4.5, 0.5, 0.5))
plot(x_lin, y_lin, pch = 16, xlab = "x", ylab = "y")
lines(x_lin, y_det, lwd = 2, col = "steelblue4")
plot(x_hyp, y_hyp, pch = 16, xlab = "Population density", ylab = "Offspring per female")
curve(a * b / (b + x), add = TRUE, lwd = 2, col = "steelblue4")
par(mfrow = c(1, 1))
#  [PAUSE] exponential decline with normal, Poisson, lognormal noise: which spreads shrink with the mean?

# 5. Many replicates: could this model have produced our data? -----------------

tree_data <- data.frame(girth = 1:10,
                        volume = c(0.4, 2.7, 5.3, 13.3, 42.4,
                                   63.5, 45.8, 233.4, 213.7, 383.1))
plot(volume ~ girth, data = tree_data, pch = 16, xlab = "Girth", ylab = "Volume")

#  straight line: Volume ~ Normal(10 + 40 girth, 31)
set.seed(31)
n_sim <- 1000
sim_vol <- replicate(n_sim, rnorm(10, mean = 10 + 40 * tree_data$girth, sd = 31))
boxplot(t(sim_vol), names = tree_data$girth, col = "grey90",
        xlab = "Girth", ylab = "Volume (simulated)")
points(1:10, tree_data$volume, pch = 16, col = "red", cex = 1.5)

#  log(Volume) ~ Normal(log 0.3 + 3 log girth, 0.4), and Volume ~ Normal(100, 400)
sim_cube <- replicate(n_sim, exp(rnorm(10, log(0.3) + 3 * log(tree_data$girth), 0.4)))
sim_null <- replicate(n_sim, rnorm(10, mean = 100, sd = 400))

par(mfrow = c(1, 2), mar = c(4, 4.5, 2, 0.5))
boxplot(t(sim_cube), names = tree_data$girth, col = "grey90", log = "y",
        xlab = "Girth", ylab = "Volume (log axis)", main = "Cube of girth, lognormal noise")
points(1:10, tree_data$volume, pch = 16, col = "red", cex = 1.3)
boxplot(t(sim_null), names = tree_data$girth, col = "grey90",
        xlab = "Girth", ylab = "Volume", main = "No relationship, wide noise")
points(1:10, tree_data$volume, pch = 16, col = "red", cex = 1.3)
par(mfrow = c(1, 1))
#  [PAUSE] the no-relationship model also covers every point. is it as good as the geometric one?

# 6. A null hypothesis is a data-generating model ------------------------------

t8 <- yr[1:8]
y8 <- lek[1:8]
slope_obs <- unname(coef(lm(y8 ~ t8))[2])         # males per year
slope_obs

#  no trend: constant mean (the 8-year average), same NegBin noise (k = 8)
set.seed(5)
null_slopes <- replicate(5000, {
  y_null <- rnbinom(8, mu = mean(y8), size = 8)
  coef(lm(y_null ~ t8))[2]
})
p_8 <- mean(null_slopes <= slope_obs)             # fraction at least as steep
p_8

hist(null_slopes, breaks = 50, col = "grey85", border = "white", main = "",
     xlab = "Slope under the no-trend model (males per year)")
abline(v = slope_obs, col = "red", lwd = 2)

#  same simulation with ten years
t10 <- yr[1:10]
y10 <- lek[1:10]
slope_10 <- unname(coef(lm(y10 ~ t10))[2])
set.seed(5)
null_10 <- replicate(5000, coef(lm(rnbinom(10, mu = mean(y10), size = 8) ~ t10))[2])
mean(null_10 <= slope_10)

# TRENDS ON THE LOG SCALE =======================================================

# 1. The lognormal distribution ------------------------------------------------
#  median = exp(m), mean = exp(m + s^2 / 2), SD = sqrt(exp(s^2) - 1) * mean

s <- 0.15
median_I <- c(40, 24)                             # year 0 and year 20
mean_I <- median_I * exp(s^2 / 2)
sd_I <- mean_I * sqrt(exp(s^2) - 1)
data.frame(year = c(0, 20), median = median_I,
           mean = round(mean_I, 2), sd = round(sd_I, 2))

# 2. Fitting an exponential trend with a straight line --------------------------
#  I_t = I0 lambda^t exp(e_t)  ->  log I_t = log I0 + b t + e_t, b = log(lambda)

I0 <- 40
sigma <- 0.15
b_true <- log(0.6) / 20                           # 40% decline over 20 years
t_yr <- 1:20

set.seed(2201)
idx <- I0 * exp(b_true * t_yr + rnorm(20, 0, sigma))
fit_idx <- lm(log(idx) ~ t_yr)
coef(fit_idx)

b_hat <- unname(coef(fit_idx)["t_yr"])
par(mfrow = c(1, 2), mar = c(4, 4.5, 0.5, 0.5))
plot(t_yr, idx, pch = 16, ylim = c(0, 60), xlab = "Year", ylab = "Index")
curve(I0 * exp(b_true * x), add = TRUE, lty = 2, lwd = 2)
curve(exp(coef(fit_idx)[1] + b_hat * x), add = TRUE, lwd = 2, col = "steelblue4")
plot(t_yr, log(idx), pch = 16, xlab = "Year", ylab = "log(index)")
abline(log(I0), b_true, lty = 2, lwd = 2)
abline(fit_idx, lwd = 2, col = "steelblue4")
par(mfrow = c(1, 1))

exp(b_hat)                                        # fraction kept each year
1 - exp(20 * b_hat)                               # estimated decline over 20 years

# 3. What a confidence interval promises ----------------------------------------

confint(fit_idx, "t_yr", level = 0.90)
#  [PAUSE] widget: how close to 90% is the coverage? what changes under process noise?
