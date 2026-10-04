
#  NRES 746, Week 7, Lecture 1: How sure are we? -----------------------------
#  Uncertainty from the shape of the likelihood.
#  Run section by section in RStudio (use the outline to jump around).
#  Needs the 'emdbook' package for the last section.

# Pika data and the NLL function ---------------------------------------------

haypiles <- c(2, 5, 3, 0, 4, 4)   # active haypiles in six talus plots
n_plots <- length(haypiles)
lambda_hat <- mean(haypiles)      # MLE for a Poisson mean

nll_pois <- function(lambda, y) {
  -sum(dpois(y, lambda = lambda, log = TRUE))
}

lambda_grid <- seq(0.01, 8, by = 0.01)
lik_vals <- sapply(lambda_grid, function(l) prod(dpois(haypiles, l)))

# Bayesian preview: prior x likelihood ---------------------------------------

# posterior = likelihood * prior / (area under likelihood * prior)
step <- diff(lambda_grid)[1]
post_from_prior <- function(prior_vals) {
  numer <- lik_vals * prior_vals
  numer / sum(numer * step)
}

prior_flat <- rep(1, length(lambda_grid))
prior_info <- dgamma(lambda_grid, shape = 10, rate = 5)   # prior mean = 2
post_flat <- post_from_prior(prior_flat)
post_info <- post_from_prior(prior_info)
lik_scaled <- lik_vals / max(lik_vals) * max(post_flat)   # rescaled to share the axes

par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.5))
plot(lambda_grid, post_flat, type = "l", lwd = 2, ylim = c(0, 0.95),
     xlab = "lambda", ylab = "Density", main = "Flat prior", cex.main = 0.9)
lines(lambda_grid, lik_scaled, lwd = 2, lty = 2, col = "darkorange")
plot(lambda_grid, post_info, type = "l", lwd = 2, ylim = c(0, 0.95),
     xlab = "lambda", ylab = "Density", main = "Prior centered on 2", cex.main = 0.9)
lines(lambda_grid, prior_info, lwd = 2, col = "grey60")
lines(lambda_grid, lik_scaled, lwd = 2, lty = 2, col = "darkorange")
legend("topright", c("posterior", "prior", "likelihood (rescaled)"),
       lwd = 2, lty = c(1, 1, 2), col = c("black", "grey60", "darkorange"),
       bty = "n", cex = 0.7)
par(mfrow = c(1, 1))

lambda_grid[which.max(post_info)]   # posterior mode with the informative prior

# flat prior: posterior is Gamma(shape = sum(y) + 1, rate = n)
cred_int <- qgamma(c(0.025, 0.975), shape = sum(haypiles) + 1, rate = n_plots)
round(cred_int, 2)

# Curvature is information: 6 plots vs. 24 plots ----------------------------

haypiles_24 <- rep(haypiles, 4)   # same mean, four times the data
dnll_6 <- sapply(lambda_grid, nll_pois, y = haypiles) - nll_pois(3, haypiles)
dnll_24 <- sapply(lambda_grid, nll_pois, y = haypiles_24) - nll_pois(3, haypiles_24)
lr_cut <- qchisq(0.95, df = 1) / 2   # 3.84 / 2 = 1.92

par(mar = c(4, 4.5, 0.5, 0.5))
plot(lambda_grid, dnll_6, type = "l", lwd = 2, xlim = c(0.5, 7), ylim = c(0, 6),
     xlab = "lambda (haypiles per plot)", ylab = "NLL minus its minimum")
lines(lambda_grid, dnll_24, lwd = 2, lty = 2)
abline(h = lr_cut, col = "red")
legend("topright", c("6 plots", "24 plots"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)

# Likelihood ratio interval ---------------------------------------------------

lr_cut

# the interval ends are where NLL(lambda) - NLL(MLE) = 1.92
f_cut <- function(l, y) nll_pois(l, y) - nll_pois(mean(y), y) - lr_cut
lr_int <- c(uniroot(f_cut, c(0.5, 3), y = haypiles)$root,
            uniroot(f_cut, c(3, 8), y = haypiles)$root)
lr_int_24 <- c(uniroot(f_cut, c(1, 3), y = haypiles_24)$root,
               uniroot(f_cut, c(3, 6), y = haypiles_24)$root)
round(rbind(plots_6 = lr_int, plots_24 = lr_int_24), 2)

round(rbind(likelihood_ratio = lr_int, bayes_flat_prior = cred_int), 2)

# Quadratic approximation and the Hessian -------------------------------------

# NLL(lambda) ~ NLL(MLE) + 0.5 * NLL''(MLE) * (lambda - MLE)^2
# SE = 1 / sqrt(NLL''(MLE))
fit <- optim(par = 2, fn = nll_pois, y = haypiles, method = "Brent",
             lower = 0.01, upper = 20, hessian = TRUE)
fit$par
fit$hessian   # curvature at the MLE

se_lambda <- 1 / sqrt(fit$hessian[1, 1])
wald_int <- fit$par + c(-1, 1) * qnorm(0.975) * se_lambda
round(c(se = se_lambda, lower = wald_int[1], upper = wald_int[2]), 2)

quad_approx <- 0.5 * fit$hessian[1, 1] * (lambda_grid - lambda_hat)^2
plot(lambda_grid, dnll_6, type = "l", lwd = 2, xlim = c(0.5, 7), ylim = c(0, 5),
     xlab = "lambda (haypiles per plot)", ylab = "NLL minus its minimum")
lines(lambda_grid, quad_approx, lwd = 2, lty = 2, col = "blue")
abline(h = lr_cut, col = "red")
abline(v = lr_int, lty = 3)
abline(v = wald_int, lty = 3, col = "blue")
legend(x = 4.5, y = 1.3, c("true NLL", "quadratic"), lwd = 2, lty = c(1, 2),
       cex = 0.8, col = c("black", "blue"), bty = "n")

# Two parameters: slice vs. profile (Bolker's tadpole data) -------------------

library(emdbook)
data(ReedfrogFuncresp)
initial_n <- ReedfrogFuncresp$Initial
killed <- ReedfrogFuncresp$Killed

# Holling type II: p = a / (1 + a h N); number killed is binomial
nll_holling <- function(a, h) {
  p <- a / (1 + a * h * initial_n)
  -sum(dbinom(killed, size = initial_n, prob = p, log = TRUE))
}

tad_fit <- optim(c(0.5, 0.0125), function(par) nll_holling(par[1], par[2]),
                 method = "L-BFGS-B", lower = c(0.01, 0.0001),
                 upper = c(2, 0.1), hessian = TRUE)
a_hat <- tad_fit$par[1]
h_hat <- tad_fit$par[2]
min_nll <- tad_fit$value
round(tad_fit$par, 4)

a_seq <- seq(0.33, 0.8, length.out = 121)
h_seq <- seq(0.001, 0.035, length.out = 121)
nll_surf <- outer(a_seq, h_seq, Vectorize(nll_holling))

# profile: for each a, re-optimize h; slice: hold h at its MLE
profile_h <- sapply(a_seq, function(a) {
  optimize(function(h) nll_holling(a, h), c(0.0001, 0.1))$minimum
})
profile_a <- sapply(seq_along(a_seq), function(i) nll_holling(a_seq[i], profile_h[i]))
slice_a <- sapply(a_seq, nll_holling, h = h_hat)

par(mfrow = c(1, 2), mar = c(4, 4.5, 0.5, 0.5))
contour(a_seq, h_seq, nll_surf - min_nll, levels = c(1, 2, 3, 5, 8, 12),
        xlab = "attack rate a", ylab = "handling time h", col = "grey50")
contour(a_seq, h_seq, nll_surf - min_nll, levels = qchisq(0.95, 2) / 2,
        add = TRUE, lwd = 2, drawlabels = FALSE)   # 95% joint region
abline(h = h_hat, lty = 2, col = "purple")
lines(a_seq, profile_h, lwd = 2, col = "darkgreen")
points(a_hat, h_hat, pch = 16)
plot(a_seq, slice_a - min_nll, type = "l", lwd = 2, lty = 2, col = "purple",
     ylim = c(0, 5), xlab = "attack rate a", ylab = "NLL minus its minimum")
lines(a_seq, profile_a - min_nll, lwd = 2, col = "darkgreen")
abline(h = lr_cut, col = "red")
legend(x = 0.63, y = 1.1, c("slice", "profile"), lwd = 2, lty = c(2, 1),
       cex = 0.8, seg.len = 1.2, col = c("purple", "darkgreen"), bty = "n")
par(mfrow = c(1, 1))

# Hessian matrix -> variance-covariance matrix ---------------------------------

vcov_tad <- solve(tad_fit$hessian)
round(sqrt(diag(vcov_tad)), 4)   # standard errors of a and h
round(cov2cor(vcov_tad), 2)      # correlation between the two estimates
