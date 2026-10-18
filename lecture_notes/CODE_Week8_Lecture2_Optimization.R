
#  NRES 746, Week 8, Lecture 2: Optimization ----------------------------------
#  Section numbers match the lecture notes. Needs the emdbook package.

# 0. Take-home debrief (~8 min) ----------------------------------------------
#  [BOARD] students' Newton steps; collect fixes -> [CODE]

nll_pika <- function(lambda) 6 * lambda - 18 * log(lambda)
d1_pika <- function(lambda) 6 - 18 / lambda
d2_pika <- function(lambda) 18 / lambda^2

lambda <- 2
for (i in seq_len(4)) {
  lambda <- lambda - d1_pika(lambda) / d2_pika(lambda)
  print(lambda)
}
7 - d1_pika(7) / d2_pika(7)   # overshoot

# Newton on theta = log(lambda): NLL = 6 exp(theta) - 18 theta
theta <- log(7)
for (i in seq_len(4)) {
  theta <- theta - (6 * exp(theta) - 18) / (6 * exp(theta))
  print(exp(theta))
}

# 1. Big picture: the workflow (~5 min) ----------------------------------------
#  [BOARD] full workflow, first time -> [CODE] myxo surface

library(emdbook)
titer <- subset(MyxoTiter_sum, grade == 1)$titer
hist(titer, freq = FALSE, main = "", xlab = "Titer")

# gamma NLL; Inf outside the allowed range so optimizers can't use it
nll_gamma <- function(p) {
  if (any(p <= 0)) return(Inf)
  -sum(dgamma(titer, shape = p[1], scale = p[2], log = TRUE))
}

shape_vec <- seq(10, 90, length.out = 120)
scale_vec <- seq(0.03, 0.32, length.out = 120)
surf <- outer(shape_vec, scale_vec, Vectorize(function(a, s) nll_gamma(c(a, s))))
length(surf)   # evaluations for this picture

# contour plot of NLL minus its minimum over the shape-scale grid
plot_surface <- function(...) {
  contour(shape_vec, scale_vec, surf - min(surf),
          levels = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
          xlab = "shape", ylab = "scale", col = "grey60", labcex = 0.6, ...)
}

fit_nm <- optim(c(20, 0.05), nll_gamma)
fit_nm$par
plot_surface()
points(fit_nm$par[1], fit_nm$par[2], pch = 4, cex = 1.5, lwd = 2)

30^10 * 0.001 / (60 * 60 * 24 * 365)   # years at 1 ms per evaluation

hist(titer, freq = FALSE, main = "", xlab = "Titer")
curve(dgamma(x, shape = fit_nm$par[1], scale = fit_nm$par[2]), add = TRUE, lwd = 2)

# 2. Newton in 2D (~10 min) ----------------------------------------------------
#  [BOARD] gradient, Hessian, Newton update -> [CODE]

# central differences: f'(x) ~ [f(x + delta) - f(x - delta)] / (2 delta)
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
round(cbind(path_good, nll = apply(path_good, 1, nll_gamma)), 3)
round(cbind(path_bad, nll = apply(path_bad, 1, nll_gamma)), 3)

plot_surface()
lines(path_good, type = "o", pch = 19, cex = 0.7, col = "steelblue", lwd = 2)
lines(path_bad, type = "o", pch = 19, cex = 0.7, col = "darkorange", lwd = 2)

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
#  [PAUSE] is the BFGS fit done?

# 5. Wrap (~1 min) -------------------------------------------------------------
#  [BOARD] point at the FIT box; hand out optional Curvature by Hand WOD


# Not covered Oct. 14: moved to CODE_Week9_MCMCAndModelSelection.R ----------
#  Nelder-Mead and simulated annealing (Week 9 Monday); reparameterizing (Lab 7)
