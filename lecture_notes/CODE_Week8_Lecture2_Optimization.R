
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

# 3. Nelder-Mead (~10 min) -----------------------------------------------------
#  [BOARD] sketch triangle moves -> [CODE]

# hand-coded Nelder-Mead for two parameters; stores every simplex
nelder_mead <- function(f, start, step, n_iter = 80) {
  simplex <- rbind(start, start + c(step[1], 0), start + c(0, step[2]))
  vals <- apply(simplex, 1, f)
  history <- vector("list", n_iter)
  moves <- character(n_iter)
  for (it in seq_len(n_iter)) {
    o <- order(vals)
    simplex <- simplex[o, ]
    vals <- vals[o]
    history[[it]] <- simplex
    centroid <- colMeans(simplex[1:2, ])
    refl <- centroid + (centroid - simplex[3, ])
    f_r <- f(refl)
    if (f_r < vals[1]) {
      expd <- centroid + 2 * (centroid - simplex[3, ])
      f_e <- f(expd)
      if (f_e < f_r) {
        simplex[3, ] <- expd; vals[3] <- f_e; moves[it] <- "expand"
      } else {
        simplex[3, ] <- refl; vals[3] <- f_r; moves[it] <- "reflect"
      }
    } else if (f_r < vals[2]) {
      simplex[3, ] <- refl; vals[3] <- f_r; moves[it] <- "reflect"
    } else {
      contr <- centroid + 0.5 * (simplex[3, ] - centroid)
      f_c <- f(contr)
      if (f_c < vals[3]) {
        simplex[3, ] <- contr; vals[3] <- f_c; moves[it] <- "contract"
      } else {
        simplex[2:3, ] <- (simplex[2:3, ] + rbind(simplex[1, ], simplex[1, ])) / 2
        vals[2:3] <- apply(simplex[2:3, ], 1, f)
        moves[it] <- "shrink"
      }
    }
  }
  list(best = simplex[which.min(vals), ], value = min(vals),
       history = history, moves = moves)
}

nm <- nelder_mead(nll_gamma, c(20, 0.05), step = c(5, 0.02))
round(c(shape = nm$best[1], scale = nm$best[2], nll = nm$value), 4)
table(nm$moves)

move_cols <- c(expand = "darkorange", reflect = "steelblue",
               contract = "firebrick", shrink = "purple")
plot_surface()
for (i in 1:30) polygon(nm$history[[i]], border = move_cols[nm$moves[i]], lwd = 1.5)
legend("topright", names(move_cols), col = move_cols, lwd = 2, bty = "n")

# correlation between shape and scale; Gamma mean = shape x scale
vc <- solve(optimHess(fit_nm$par, nll_gamma))
cov2cor(vc)
prod(fit_nm$par)
mean(titer)
#  [PAUSE] what can these data tell us?

# 4. Simulated annealing (~8 min; rest to Week 9) -----------------------------
#  [BOARD] two valleys; acceptance rule -> [CODE]

accept_tab <- outer(c(10, 1, 0.1), c(0.5, 2, 10),
                    function(temp, delta) exp(-delta / temp))
dimnames(accept_tab) <- list(paste("T =", c(10, 1, 0.1)),
                             paste("NLL worse by", c(0.5, 2, 10)))
signif(accept_tab, 2)

# simulated annealing: accept worse steps w.p. exp(-delta / temp); cool temp over time
anneal <- function(f, start, step_sd, n_steps = 3000, temp0 = 10,
                   cool = 0.9, every = 100) {
  current <- start
  f_current <- f(current)
  best <- current
  f_best <- f_current
  temp <- temp0
  track <- matrix(NA, n_steps, 4,
                  dimnames = list(NULL, c("shape", "scale", "nll", "temp")))
  for (i in seq_len(n_steps)) {
    proposal <- current + rnorm(length(current), 0, step_sd)
    f_prop <- f(proposal)
    if (is.finite(f_prop) && runif(1) < exp(-(f_prop - f_current) / temp)) {
      current <- proposal
      f_current <- f_prop
    }
    if (f_current < f_best) {
      best <- current
      f_best <- f_current
    }
    track[i, ] <- c(current, f_current, temp)
    if (i %% every == 0) temp <- temp * cool
  }
  list(best = best, value = f_best, track = track)
}

set.seed(746)
sa <- anneal(nll_gamma, c(20, 0.05), step_sd = c(2, 0.005))
round(c(shape = sa$best[1], scale = sa$best[2], nll = sa$value), 4)

par(mfrow = c(1, 2))
plot_surface()
lines(sa$track[, "shape"], sa$track[, "scale"], col = adjustcolor("darkorange", 0.5))
plot(sa$track[, "nll"], type = "l", log = "y", xlab = "Step", ylab = "Current NLL")
par(mfrow = c(1, 1))

# optim's built-in annealer, default settings
set.seed(746)
fit_sann <- optim(c(20, 0.05), nll_gamma, method = "SANN")
fit_sann[c("par", "value", "counts", "convergence")]
#  [PAUSE] is SANN broken?

# 5. Wrap (~1 min) -------------------------------------------------------------
#  [BOARD] point at the FIT box

# If time: reparameterize ------------------------------------------------------

nll_mean_shape <- function(p) nll_gamma(c(p[2], p[1] / p[2]))
fit_ms <- optim(c(7, 20), nll_mean_shape)
fit_ms$par
cov2cor(solve(optimHess(fit_ms$par, nll_mean_shape)))
