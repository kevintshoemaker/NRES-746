
#  NRES 746, Week 8: Priors, posteriors, and derivative-based optimization ----
#  Code from Week8_PosteriorsAndOptimization.html. Package: emdbook.

# PRIORS AND POSTERIORS ========================================================

# 2. WOD: same data, three priors ----------------------------------------------
#  posterior = Beta(alpha + k, beta + n - k); 21 of 30 juvenile tortoises survived
#  [PAUSE] is the predator-proof pen a reasonable prior? what would you report?

# 3. Worked example: germination with two priors ------------------------------

k <- 30
n <- 40
prior_flat <- c(1, 1)
prior_info <- c(121, 81)
post_flat <- prior_flat + c(k, n - k)   # Beta(31, 11)
post_info <- prior_info + c(k, n - k)   # Beta(151, 91)

# mean = alpha / (alpha + beta); mode = (alpha - 1) / (alpha + beta - 2)
beta_mean <- function(ab) ab[1] / sum(ab)
beta_mode <- function(ab) (ab[1] - 1) / (sum(ab) - 2)

round(c(flat_mean = beta_mean(post_flat), flat_mode = beta_mode(post_flat),
        supplier_mean = beta_mean(post_info)), 3)

p_seq <- seq(0.001, 0.999, length.out = 500)
lik_scaled <- dbeta(p_seq, k + 1, n - k + 1)   # likelihood rescaled to area 1

# draws prior, rescaled likelihood, and posterior on one panel
plot_panel <- function(prior_ab, post_ab, title) {
  plot(p_seq, dbeta(p_seq, post_ab[1], post_ab[2]), type = "l", lwd = 2,
       ylim = c(0, 14), xlab = "Germination probability, p", ylab = "Density",
       main = title)
  lines(p_seq, dbeta(p_seq, prior_ab[1], prior_ab[2]), lwd = 2, col = "grey60")
  lines(p_seq, lik_scaled, lwd = 2, lty = 2, col = "darkorange")
}

par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.5))
plot_panel(prior_flat, post_flat, "Flat prior, Beta(1, 1)")
legend("topleft", c("prior", "likelihood (rescaled)", "posterior"),
       lwd = 2, lty = c(1, 2, 1), col = c("grey60", "darkorange", "black"),
       bty = "n", cex = 0.8)
plot_panel(prior_info, post_info, "Supplier prior, Beta(121, 81)")
par(mfrow = c(1, 1))

beta_mean(prior_info + c(300, 100))   # 300 of 400 germinate, same prior

post_two_steps <- (prior_info + c(k, n - k)) + c(k, n - k)   # 40 seeds, then 40 more
post_one_step <- prior_info + c(2 * k, 2 * (n - k))         # 80 seeds at once
rbind(post_two_steps, post_one_step)
#  [PAUSE] how many seeds (at 75%) to land halfway between 0.60 and 0.75? (widget)

# 4. Summarizing a posterior --------------------------------------------------

summ <- rbind(
  flat = c(mean = beta_mean(post_flat),
           median = qbeta(0.5, post_flat[1], post_flat[2]),
           mode = beta_mode(post_flat),
           lower95 = qbeta(0.025, post_flat[1], post_flat[2]),
           upper95 = qbeta(0.975, post_flat[1], post_flat[2])),
  supplier = c(mean = beta_mean(post_info),
               median = qbeta(0.5, post_info[1], post_info[2]),
               mode = beta_mode(post_info),
               lower95 = qbeta(0.025, post_info[1], post_info[2]),
               upper95 = qbeta(0.975, post_info[1], post_info[2]))
)
round(summ, 3)

# shortest interval holding a share 'mass' of the posterior draws
hpd <- function(draws, mass = 0.95) {
  draws <- sort(draws)
  n_in <- floor(mass * length(draws))
  widths <- draws[(n_in + 1):length(draws)] - draws[1:(length(draws) - n_in)]
  i <- which.min(widths)
  c(lower = draws[i], upper = draws[i + n_in])
}

set.seed(746)
draws_flat <- rbeta(1e5, post_flat[1], post_flat[2])
round(rbind(central = summ["flat", c("lower95", "upper95")],
            hpd = hpd(draws_flat)), 3)

# 5. Where new Bayesians get stuck --------------------------------------------

integrate(function(p) dbinom(k, n, p), 0, 1)$value   # 1 / (n + 1)

# a flat prior on p is peaked on the log-odds scale
set.seed(746)
p_draws <- runif(1e5)
par(mfrow = c(1, 2))
hist(p_draws, breaks = 40, freq = FALSE, col = "grey80", border = "white",
     main = "Flat prior on p", xlab = "p")
hist(qlogis(p_draws), breaks = 60, freq = FALSE, col = "grey80", border = "white",
     main = "Same prior, log-odds scale", xlab = "log(p / (1 - p))", xlim = c(-8, 8))
curve(dlogis(x), add = TRUE, lwd = 2)
par(mfrow = c(1, 1))

# pikas: likelihood ratio interval vs. flat-prior credible interval
haypiles <- c(2, 5, 3, 0, 4, 4)
nll_pois <- function(lambda) -sum(dpois(haypiles, lambda, log = TRUE))
f_cut <- function(l) nll_pois(l) - nll_pois(mean(haypiles)) - qchisq(0.95, 1) / 2
lr_int <- c(uniroot(f_cut, c(0.5, 3))$root, uniroot(f_cut, c(3, 8))$root)

# flat prior on lambda: posterior is Gamma(shape = sum(y) + 1, rate = n)
cred_int <- qgamma(c(0.025, 0.975), shape = sum(haypiles) + 1, rate = length(haypiles))
round(rbind(likelihood_ratio = lr_int, bayes_flat_prior = cred_int), 2)
qgamma(c(0.025, 0.975), 18, 6)   # prior flat on log(lambda) instead

# 6. Why Bayesian methods are so widely used ----------------------------------

odds_draws <- draws_flat / (1 - draws_flat)
round(quantile(odds_draws, c(0.025, 0.5, 0.975)), 2)
mean(draws_flat > 0.6)   # P(germination rate > 0.6), flat prior

# Discussion ------------------------------------------------------------------
#  [DISCUSS] supplier's prior; Bayesian interval for -a/b; "95% chance" for the pikas

# Bonus: when there's no conjugate prior --------------------------------------

library(emdbook)
titer <- subset(MyxoTiter_sum, grade == 1)$titer

shape_g <- seq(10, 160, length.out = 300)
scale_g <- seq(0.03, 0.45, length.out = 300)
log_lik <- outer(shape_g, scale_g, Vectorize(function(a, s)
  sum(dgamma(titer, shape = a, scale = s, log = TRUE))))
post_grid <- exp(log_lik - max(log_lik))   # flat prior over the grid
post_grid <- post_grid / sum(post_grid)

# draw parameter pairs from the grid in proportion to their posterior probability
set.seed(746)
idx <- sample(length(post_grid), 10000, replace = TRUE, prob = post_grid)
draws_grid <- cbind(shape = shape_g[row(post_grid)[idx]],
                    scale = scale_g[col(post_grid)[idx]])
apply(draws_grid, 2, quantile, c(0.025, 0.5, 0.975))

plot(draws_grid[1:2000, ], pch = 16, cex = 0.4, col = adjustcolor("steelblue", 0.3),
     xlim = range(shape_g), ylim = range(scale_g), xlab = "shape", ylab = "scale")
contour(shape_g, scale_g, post_grid, levels = max(post_grid) * c(0.05, 0.3, 0.7),
        drawlabels = FALSE, add = TRUE)

# Bonus: posteriors with two peaks --------------------------------------------

# equal mixture of Normal(-2, 0.5) and Normal(2, 0.5)
dens_mix <- function(x) 0.5 * dnorm(x, -2, 0.5) + 0.5 * dnorm(x, 2, 0.5)
x_mix <- seq(-5, 5, length.out = 20001)
dx <- diff(x_mix)[1]
d_mix <- dens_mix(x_mix)
o <- order(d_mix, decreasing = TRUE)
keep <- sort(x_mix[o[cumsum(d_mix[o] * dx) <= 0.95]])   # highest-density points, 95% of the area
gaps <- which(diff(keep) > 2 * dx)
round(cbind(lower = keep[c(1, gaps + 1)], upper = keep[c(gaps, length(keep))]), 2)

set.seed(1)
draws_mix <- c(rnorm(5e4, -2, 0.5), rnorm(5e4, 2, 0.5))
round(rbind(shortest_single = hpd(draws_mix),
            central = quantile(draws_mix, c(0.025, 0.975))), 2)
mean(draws_mix)

# DERIVATIVE-BASED OPTIMIZATION ================================================

# 1. Big picture --------------------------------------------------------------

# gamma NLL; Inf outside the allowed range so optimizers can't use it
nll_gamma <- function(p) {
  if (any(p <= 0)) return(Inf)
  -sum(dgamma(titer, shape = p[1], scale = p[2], log = TRUE))
}

shape_vec <- seq(10, 90, length.out = 120)
scale_vec <- seq(0.03, 0.32, length.out = 120)
surf <- outer(shape_vec, scale_vec, Vectorize(function(a, s) nll_gamma(c(a, s))))

# contour plot of NLL minus its minimum
plot_surface <- function(...) {
  contour(shape_vec, scale_vec, surf - min(surf),
          levels = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
          xlab = "shape", ylab = "scale", col = "grey60", labcex = 0.6, ...)
}

fit_nm <- optim(c(20, 0.05), nll_gamma)
plot_surface()
points(fit_nm$par[1], fit_nm$par[2], pch = 4, cex = 1.5, lwd = 2)

# 2. Newton's method in one dimension -----------------------------------------
#  pikas: NLL = 6 lambda - 18 log(lambda); lambda* = lambda - NLL' / NLL''

d1 <- function(l) 6 - 18 / l
d2 <- function(l) 18 / l^2

guess <- 2
for (i in seq_len(4)) {
  guess <- guess - d1(guess) / d2(guess)
  print(guess)
}

7 - d1(7) / d2(7)   # starting at 7 overshoots

# on the log scale, theta = log(lambda): NLL' = 6 e^theta - 18, NLL'' = 6 e^theta
theta <- log(7)
for (i in seq_len(4)) {
  theta <- theta - (6 * exp(theta) - 18) / (6 * exp(theta))
  print(exp(theta))
}
#  [PAUSE] why can't the log-scale version overshoot past zero? (widget)

# 3. Slopes without calculus --------------------------------------------------

nll_pika <- function(l) 6 * l - 18 * log(l)
delta <- 0.5
vals <- nll_pika(c(1.5, 2, 2.5))
fwd <- (vals[3] - vals[2]) / delta
c(forward_at_2 = fwd, exact_at_2 = d1(2), exact_at_2.25 = d1(2.25))

# f''(x) ~ [f(x + delta) - 2 f(x) + f(x - delta)] / delta^2
c(second_difference = (vals[3] - 2 * vals[2] + vals[1]) / delta^2,
  exact = d2(2))

# 4. Newton's method in two dimensions, and BFGS ------------------------------

# central differences for each parameter: f'(x) ~ [f(x + d) - f(x - d)] / (2 d)
grad_fd <- function(f, p, delta = 1e-4 * abs(p)) {
  sapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, delta[i])
    (f(p + e) - f(p - e)) / (2 * delta[i])
  })
}

# matrix of second derivatives (the Hessian), by finite differences
hess_fd <- function(f, p, delta = 1e-3 * abs(p)) {
  k <- length(p)
  H <- matrix(0, k, k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      ei <- replace(numeric(k), i, delta[i])
      ej <- replace(numeric(k), j, delta[j])
      H[i, j] <- (f(p + ei + ej) - f(p + ei - ej) - f(p - ei + ej) +
                    f(p - ei - ej)) / (4 * delta[i] * delta[j])
    }
  }
  H
}

# Newton's method; returns every point visited
newton <- function(f, start, n_steps = 8) {
  path <- matrix(start, nrow = 1)
  p <- start
  for (i in seq_len(n_steps)) {
    p <- p - solve(hess_fd(f, p), grad_fd(f, p))
    path <- rbind(path, p, deparse.level = 0)
  }
  path
}

path_good <- newton(nll_gamma, c(40, 0.2))
path_bad <- newton(nll_gamma, c(20, 0.05))
round(cbind(path_bad, nll = apply(path_bad, 1, nll_gamma)), 3)

plot_surface()
lines(path_good, type = "o", pch = 19, cex = 0.7, col = "steelblue", lwd = 2)
lines(path_bad, type = "o", pch = 19, cex = 0.7, col = "darkorange", lwd = 2)
legend("topright", c("start (40, 0.2)", "start (20, 0.05)"), lwd = 2,
       col = c("steelblue", "darkorange"), bty = "n", cex = 0.8)

fit_bfgs <- optim(c(20, 0.05), nll_gamma, method = "BFGS")
fit_bfgs_log <- optim(log(c(20, 0.05)), function(q) nll_gamma(exp(q)),
                      method = "BFGS")
cmp <- rbind(
  nelder_mead = c(fit_nm$par, fit_nm$value, fit_nm$counts[1], fit_nm$convergence),
  bfgs = c(fit_bfgs$par, fit_bfgs$value, fit_bfgs$counts[1], fit_bfgs$convergence),
  bfgs_log_scale = c(exp(fit_bfgs_log$par), fit_bfgs_log$value,
                     fit_bfgs_log$counts[1], fit_bfgs_log$convergence)
)
colnames(cmp) <- c("shape", "scale", "nll", "n_evals", "convergence")
round(cmp, 4)
#  [PAUSE] is the BFGS fit done? what would you do?

# Discussion ------------------------------------------------------------------
#  [DISCUSS] two starts, two NLLs; the "replace negatives with 0.001" fix; why not a tiny step

# Bonus: the shape of the valley ----------------------------------------------

fit_h <- optim(fit_nm$par, nll_gamma, hessian = TRUE)
vcov_myxo <- solve(fit_h$hessian)
sqrt(diag(vcov_myxo))         # approximate standard errors
cov2cor(vcov_myxo)[1, 2]      # correlation of the two estimates
prod(fit_nm$par)              # shape x scale = the Gamma mean
mean(titer)

nll_mean_shape <- function(p) nll_gamma(c(p[2], p[1] / p[2]))
fit_ms <- optim(c(7, 20), nll_mean_shape, hessian = TRUE)
fit_ms$par
cov2cor(solve(fit_ms$hessian))[1, 2]
