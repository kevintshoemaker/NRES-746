#  NRES 746, Week 7: Likelihood intervals and Bayes' rule for parameters ------
#  Code from Week7_IntervalsAndBayes.html. Package: emdbook.

# UNCERTAINTY FROM THE SHAPE OF THE LIKELIHOOD =================================

# 2. WOD: reading the interval ---------------------------------------------------
#  2 of 20 frogs chytrid-positive; binomial NLL minus its minimum

nll_chytrid <- function(p) -dbinom(2, size = 20, prob = p, log = TRUE)
p_grid <- seq(0.001, 0.4, by = 0.001)
plot(p_grid, nll_chytrid(p_grid) - nll_chytrid(0.1), type = "l", lwd = 2,
     ylim = c(0, 5), xlab = "Prevalence, p", ylab = "NLL minus its minimum")
abline(h = 2, lty = 2)
#  [PAUSE] is the range of plausible values symmetric around 0.1? why might it be lopsided?

# 3. Curvature is information --------------------------------------------------

haypiles <- c(2, 5, 3, 0, 4, 4)
lambda_hat <- mean(haypiles)

nll_pois <- function(lambda, y) {
  -sum(dpois(y, lambda = lambda, log = TRUE))
}

lambda_grid <- seq(0.01, 8, by = 0.01)
haypiles_24 <- rep(haypiles, 4)
dnll_6 <- sapply(lambda_grid, nll_pois, y = haypiles) - nll_pois(3, haypiles)
dnll_24 <- sapply(lambda_grid, nll_pois, y = haypiles_24) - nll_pois(3, haypiles_24)
lr_cut <- qchisq(0.95, df = 1) / 2   # 3.84 / 2 = 1.92

par(mar = c(4, 4.5, 0.5, 0.5))
plot(lambda_grid, dnll_6, type = "l", lwd = 2, xlim = c(0.5, 7), ylim = c(0, 6),
     xlab = "lambda (haypiles per plot)", ylab = "NLL minus its minimum")
lines(lambda_grid, dnll_24, lwd = 2, lty = 2)
abline(h = lr_cut, col = "red")
legend("topright", c("6 plots", "24 plots", "cutoff (1.92)"), lwd = c(2, 2, 1),
       lty = c(1, 2, 1), col = c("black", "black", "red"), bty = "n")

c(dbinom(1, 2, 0.5), dbinom(500, 1000, 0.5))   # more data, smaller maximum likelihood

# 4. The likelihood ratio interval ---------------------------------------------

# zero where the NLL is 1.92 above its minimum
f_cut <- function(l, y) nll_pois(l, y) - nll_pois(mean(y), y) - lr_cut
lr_int <- c(uniroot(f_cut, c(0.5, 3), y = haypiles)$root,
            uniroot(f_cut, c(3, 8), y = haypiles)$root)
lr_int_24 <- c(uniroot(f_cut, c(1, 3), y = haypiles_24)$root,
               uniroot(f_cut, c(3, 6), y = haypiles_24)$root)
round(rbind(plots_6 = lr_int, plots_24 = lr_int_24), 2)

# same cutoff for the chytrid data
f_cut_chytrid <- function(p) nll_chytrid(p) - nll_chytrid(0.1) - lr_cut
round(c(uniroot(f_cut_chytrid, c(1e-6, 0.1))$root,
        uniroot(f_cut_chytrid, c(0.1, 0.9))$root), 3)

# Discussion ---------------------------------------------------------------------
#  [DISCUSS] is lambda = 1.5 disproven? is lambda = 3 proven?
#  [DISCUSS] six plots: trust this interval? how would you check?

# CURVATURE, THE HESSIAN, AND TWO-PARAMETER LIKELIHOODS =========================

# 1. The quadratic approximation -------------------------------------------------
#  NLL(lambda) ~ NLL(MLE) + 0.5 * NLL''(MLE) * (lambda - MLE)^2;  SE = 1 / sqrt(NLL'')

fit <- optim(par = 2, fn = nll_pois, y = haypiles, method = "Brent",
             lower = 0.01, upper = 20, hessian = TRUE)
fit$hessian
se_lambda <- 1 / sqrt(fit$hessian[1, 1])
wald_int <- fit$par + c(-1, 1) * qnorm(0.975) * se_lambda
round(c(se = se_lambda, lower = wald_int[1], upper = wald_int[2]), 2)

quad_approx <- 0.5 * fit$hessian[1, 1] * (lambda_grid - lambda_hat)^2
par(mar = c(4, 4.5, 0.5, 0.5))
plot(lambda_grid, dnll_6, type = "l", lwd = 2, xlim = c(0.5, 7), ylim = c(0, 5),
     xlab = "lambda (haypiles per plot)", ylab = "NLL minus its minimum")
lines(lambda_grid, quad_approx, lwd = 2, lty = 2, col = "blue")
abline(h = lr_cut, col = "red")
abline(v = lr_int, lty = 3)
abline(v = wald_int, lty = 3, col = "blue")
legend("topright", c("NLL", "quadratic approximation"), lwd = 2, lty = c(1, 2),
       col = c("black", "blue"), bty = "n")
#  [PAUSE] widget, 2 of 20 positive: what goes wrong with the Hessian interval?

# 2. More than one parameter: slices and profiles --------------------------------

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

# profile: re-optimize h for each a; slice: hold h at its MLE
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
legend("topright", c("slice", "profile"), lwd = 2, lty = c(2, 1),
       col = c("purple", "darkgreen"), bty = "n")
par(mfrow = c(1, 1))
#  [PAUSE] which curve would you use for an interval on the attack rate?

vcov_tad <- solve(tad_fit$hessian)
round(sqrt(diag(vcov_tad)), 4)   # standard errors of a and h
round(cov2cor(vcov_tad), 2)      # correlation between the two estimates

# Bonus: why the MLE is approximately normal -----------------------------------

set.seed(746)
lambda0 <- 3
lam <- seq(1.6, 5.5, length.out = 300)
score <- function(l, y) sum(y) / l - length(y)

# score curves for n_sims simulated surveys of n plots, with the predicted bells
score_panel <- function(n, n_sims = 15, ylim = c(-12, 12)) {
  info <- n / lambda0
  plot(NA, xlim = range(lam), ylim = ylim, xlab = "lambda", ylab = "Score",
       main = paste0(n, " plots (I = ", round(info, 1), ")"))
  abline(h = 0, lwd = 1.5)
  abline(v = lambda0, lty = 3, col = "grey40")
  for (i in seq_len(n_sims)) {
    y <- rpois(n, lambda0)
    lines(lam, sapply(lam, score, y = y), col = adjustcolor("steelblue", 0.6))
    points(mean(y), 0, pch = 16, col = "red", cex = 0.8)
    points(lambda0, score(lambda0, y), pch = 16, col = "steelblue4", cex = 0.6)
  }
  yy <- seq(-2.5, 2.5, length.out = 200) * sqrt(info)
  lines(lambda0 + 0.55 * dnorm(yy, 0, sqrt(info)) / dnorm(0, 0, sqrt(info)), yy,
        col = "steelblue4", lwd = 2)
  dx <- dnorm(lam, lambda0, sqrt(1 / info))
  lines(lam, ylim[1] + 5 * dx / max(dx), col = "red", lwd = 2)
}

par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.5))
score_panel(6)
score_panel(24)
par(mfrow = c(1, 1))

# LOGISTIC REGRESSION AND BAYES' RULE FOR PARAMETERS ============================

# 1. Logistic regression and the log-odds -----------------------------------------
#  p = 1 / (1 + exp(-(a + b x)))   <->   log(p / (1 - p)) = a + b x

a_ex <- -2    # nest success vs. distance from a forest edge (100s of m)
b_ex <- 0.8
x_seq <- seq(0, 6, length.out = 200)
par(mfrow = c(1, 2), mar = c(4, 4.5, 2, 0.5))
plot(x_seq, plogis(a_ex + b_ex * x_seq), type = "l", lwd = 2, ylim = c(0, 1),
     xlab = "Distance from edge (100 m)", ylab = "P(nest succeeds)",
     main = "Probability scale")
abline(h = 0.5, lty = 3)
plot(x_seq, a_ex + b_ex * x_seq, type = "l", lwd = 2,
     xlab = "Distance from edge (100 m)", ylab = "log(p / (1 - p))",
     main = "Log-odds scale")
abline(h = 0, lty = 3)
par(mfrow = c(1, 1))

qlogis(0.75)           # log-odds of p = 0.75 is log(3)
plogis(qlogis(0.75))   # plogis() undoes qlogis()
#  [PAUSE] with b = 0.8, how do the odds change per 100 m? per 200 m? the probability?

# 3. Bayes' rule, from a few hypotheses to a whole distribution --------------------
#  nests: 1, 1, 0, 1, so L(p) = p^3 (1 - p)

p3 <- c(0.25, 0.5, 0.75)
lik3 <- p3^3 * (1 - p3)
prior3 <- rep(1/3, 3)
post3 <- lik3 * prior3 / sum(lik3 * prior3)
round(rbind(likelihood = lik3, posterior = post3), 3)
sum(lik3)   # likelihoods don't add to 1

# posterior probability at each of n_grid equally spaced values of p (flat prior)
grid_post <- function(n_grid) {
  p_g <- seq(0, 1, length.out = n_grid + 2)[-c(1, n_grid + 2)]
  w <- p_g^3 * (1 - p_g)
  list(p = p_g, prob = w / sum(w), spacing = diff(p_g)[1])
}

par(mfrow = c(1, 3), mar = c(4, 4, 2, 3.2))
for (n_grid in c(3, 11, 101)) {
  g <- grid_post(n_grid)
  plot(g$p, g$prob, type = "h", lwd = if (n_grid < 20) 4 else 1,
       col = if (n_grid < 20) "black" else "grey60", xlim = c(0, 1),
       xlab = "p", ylab = "Posterior probability", main = paste(n_grid, "hypotheses"))
  if (n_grid == 101) {
    lines(g$p, dbeta(g$p, 4, 2) * g$spacing, lwd = 2, col = "darkorange")
    axis(4, at = (0:2) * g$spacing, labels = 0:2, col.axis = "darkorange")
    mtext("density", side = 4, line = 2.1, cex = 0.7, col = "darkorange")
  }
}
par(mfrow = c(1, 1))

grid_tab <- t(sapply(c(3, 11, 101), function(n_grid) {
  g <- grid_post(n_grid)
  c(n = n_grid, tallest_bar = max(g$prob), bar_over_spacing = max(g$prob) / g$spacing)
}))
round(grid_tab, 3)
dbeta(0.75, 4, 2)   # peak of the orange curve
#  [PAUSE] widget, 3 to 201 hypotheses: which readout shrinks, which settles down?

# 4. The Beta prior and conjugacy -----------------------------------------------
#  Beta(alpha, beta) prior + k of n -> Beta(alpha + k, beta + n - k); nests: Beta(4, 2)

# Discussion ---------------------------------------------------------------------
#  [DISCUSS] 201 hypotheses: P(p = 0.75) is about 0.01. no support for 0.75?
#  [DISCUSS] odds x 2.23 per 100 m: "2.2 times as likely to succeed"?

max(grid_post(201)$prob)
exp(0.8)

# Bonus: brute-force Bayes with an informative prior -----------------------------

lik_frog <- function(p) dbinom(3, 10, p)   # frog heard on 3 of 10 visits

# posterior density = likelihood x prior / (area under likelihood x prior)
posterior_fn <- function(lik, prior) {
  numer <- function(p) lik(p) * prior(p)
  marg_lik <- integrate(numer, 0, 1)$value   # P(data): one number
  function(p) numer(p) / marg_lik
}

prior_flat <- function(p) dbeta(p, 1, 1)
prior_info <- function(p) dbeta(p, 15, 5)   # like 14 detections in 18 earlier visits
post_flat <- posterior_fn(lik_frog, prior_flat)
post_info <- posterior_fn(lik_frog, prior_info)
post_mean <- function(post) integrate(function(p) p * post(p), 0, 1)$value
round(c(flat = post_mean(post_flat), informative = post_mean(post_info)), 3)

detections <- c(3, 1, 6, 2, 3, 2, 6, 1, 3, 3)   # 10 visits at each of 10 ponds
lik_ponds <- function(p) sapply(p, function(t) prod(dbinom(detections, 10, t)))
prior_strong <- function(p) dbeta(p, 150, 50)   # like 149 detections in 198 visits

post_ponds <- posterior_fn(lik_ponds, prior_info)
post_strong <- posterior_fn(lik_ponds, prior_strong)
round(c(informative = post_mean(post_ponds), strong = post_mean(post_strong)), 3)

p_seq <- seq(0.001, 0.999, length.out = 500)

# prior (grey) and posterior (black) on one panel
frog_panel <- function(prior, post, title) {
  plot(p_seq, post(p_seq), type = "l", lwd = 2, ylim = c(0, 16),
       xlab = "Detection probability, p", ylab = "Density", main = title)
  lines(p_seq, prior(p_seq), lwd = 2, col = "grey60")
}

par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.5))
frog_panel(prior_info, post_info, "Beta(15, 5) prior, 3 of 10")
legend("topleft", c("prior", "posterior"), lwd = 2, col = c("grey60", "black"),
       bty = "n")
frog_panel(prior_info, post_ponds, "Beta(15, 5) prior, 30 of 100")
frog_panel(prior_strong, post_strong, "Beta(150, 50) prior, 30 of 100")
par(mfrow = c(1, 1))

# conjugate shortcut: Beta(15 + 3, 5 + 7), Beta(15 + 30, 5 + 70), Beta(150 + 30, 50 + 70)
round(c(one_pond = 18 / (18 + 12), ten_ponds = 45 / (45 + 75),
        strong_prior = 180 / (180 + 120)), 3)
