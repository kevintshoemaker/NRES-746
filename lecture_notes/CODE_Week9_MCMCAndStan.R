
#  NRES 746, Week 9: Optimization to MCMC (Mon.); Hamiltonian MC and Stan (Wed.) =
#  R demos only; section numbers match the NOTES files. Needs emdbook, coda, cmdstanr.

# MONDAY, Oct. 19: Nelder-Mead, simulated annealing, Metropolis ===============

rm(list=ls())

# 2. Nelder-Mead (~7 min) ------------------------------------------------------

library(emdbook)
titer <- subset(MyxoTiter_sum, grade == 1)$titer   # myxo titers, 27 rabbits

# gamma NLL; Inf outside the allowed range so optimizers can't use it
nll_gamma <- function(p) {
  if (any(p <= 0)) return(Inf)
  -sum(dgamma(titer, shape = p[1], scale = p[2], log = TRUE))
}

# NLL surface over a shape x scale grid, for the background contours
shape_vec <- seq(10, 90, length.out = 120)
scale_vec <- seq(0.03, 0.32, length.out = 120)
surf <- outer(shape_vec, scale_vec, Vectorize(function(a, s) nll_gamma(c(a, s))))

# contour plot of NLL minus its minimum
plot_surface <- function(...) {
  contour(shape_vec, scale_vec, surf - min(surf),
          levels = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
          xlab = "shape", ylab = "scale", col = "grey60", labcex = 0.6, ...)
}

# hand-coded Nelder-Mead for two parameters; stores every simplex
# The algorithm, in words:
#   1. Start with a triangle (k + 1 = 3 points for k = 2 parameters) around the starting guess.
#   2. Rank the points by NLL: best, good, worst.
#   3. Reflect the worst point through the midpoint of the other two.
#   4. If the reflection beats the best point, try going twice as far (expand);
#      keep whichever of the two is lower.
#   5. Else, if the reflection beats the good point, keep it.
#   6. Else, try a point halfway between the worst point and the midpoint (contract);
#      keep it if it beats the worst point.
#   7. If even that fails, shrink the whole triangle halfway toward the best point.
#   8. Repeat from step 2 until the NLL values at the three points are nearly equal
#      (optim()'s stopping rule; here we just run a fixed number of iterations).
nelder_mead <- function(f, start, step, n_iter = 80) {
  simplex <- rbind(start, start + c(step[1], 0), start + c(0, step[2]))
  vals <- apply(simplex, 1, f)                # compute nll
  history <- vector("list", n_iter)
  moves <- character(n_iter)
  for (it in seq_len(n_iter)) {
    o <- order(vals)                                 # best first, worst last
    simplex <- simplex[o, ]
    vals <- vals[o]
    history[[it]] <- simplex
    centroid <- colMeans(simplex[1:2, ])
    refl <- centroid + (centroid - simplex[3, ])     # reflect the worst point
    f_r <- f(refl)
    if (f_r < vals[1]) {
      expd <- centroid + 2 * (centroid - simplex[3, ])   # best yet: try expanding
      f_e <- f(expd)
      if (f_e < f_r) {
        simplex[3, ] <- expd; vals[3] <- f_e; moves[it] <- "expand"
      } else {
        simplex[3, ] <- refl; vals[3] <- f_r; moves[it] <- "reflect"
      }
    } else if (f_r < vals[2]) {
      simplex[3, ] <- refl; vals[3] <- f_r; moves[it] <- "reflect"
    } else {
      contr <- centroid + 0.5 * (simplex[3, ] - centroid)   # still bad: contract
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

nm <- nelder_mead(f=nll_gamma, start=c(20, 0.05), step = c(5, 0.02))
round(c(shape = nm$best[1], scale = nm$best[2], nll = nm$value), 4)
table(nm$moves)

# first 30 triangles, colored by move
move_cols <- c(expand = "darkorange", reflect = "steelblue",
               contract = "firebrick", shrink = "purple")
plot_surface()
for (i in 1:30) polygon(nm$history[[i]], border = move_cols[nm$moves[i]], lwd = 1.5)
legend("topright", names(move_cols), col = move_cols, lwd = 2, bty = "n")

# the target: the MLE from optim(), marked with an X in the animation
mle <- optim(c(20, 0.05), nll_gamma)$par

# animate the simplex: current triangle bold (colored by its next move), earlier ones
# faded; n_show up to 80; pass xlim/ylim through ... to zoom in on later iterations
animate_nm <- function(nm, n_show = 10, pause = 0.6, target = mle, ...) {
  for (i in seq_len(n_show)) {
    dev.hold()
    plot_surface(main = paste0("Iteration ", i, " (next move: ", nm$moves[i],
                               "),  best NLL = ", round(nll_gamma(nm$history[[i]][1, ]), 2)),
                 ...)
    for (j in seq_len(i - 1)) polygon(nm$history[[j]], border = adjustcolor(move_cols[nm$moves[j]], 0.25))
    polygon(nm$history[[i]], border = move_cols[nm$moves[i]], lwd = 3)
    points(nm$history[[i]][1, , drop = FALSE], pch = 19)   # best vertex
    points(target[1], target[2], pch = 4, cex = 2, lwd = 3, col = "darkgreen")
    legend("topright", c(names(move_cols), "MLE"), col = c(move_cols, "darkgreen"),
           lwd = c(2, 2, 2, 2, 3), lty = c(1, 1, 1, 1, NA), pch = c(NA, NA, NA, NA, 4), bty = "n")
    dev.flush()
    Sys.sleep(pause)
  }
}
animate_nm(nm, n_show = 30)
#  animate_nm(nm, n_show = 80, pause = 0.25)   # the whole search
#  animate_nm(nm, n_show = 80, pause = 0.25, xlim = c(35, 65), ylim = c(0.1, 0.2))   # zoomed

# optim()'s default is Nelder-Mead
fit_nm <- optim(c(20, 0.05), nll_gamma)
fit_nm[c("par", "value", "counts", "convergence")]

# 3. Simulated annealing (~8 min) ----------------------------------------------

# acceptance probability exp(-delta / T) for a few temperatures and NLL increases
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

# left: the path on the surface. right: NLL of the chain's current point at each step (log
# y-axis). It plunges from ~1,800 at the start to under 100 within ~50 steps, then bounces
# above the minimum (37.67) because worse moves are still accepted. As T falls (x 0.9 every
# 100 steps, from 10 to ~0.5) the bounces shrink and the trace settles toward the minimum.
par(mfrow = c(1, 2))
plot_surface()
lines(sa$track[, "shape"], sa$track[, "scale"], col = adjustcolor("darkorange", 0.5))
plot(sa$track[, "nll"], type = "l", log = "y", xlab = "Step", ylab = "Current NLL")
par(mfrow = c(1, 1))

# optim()'s built-in annealer, default settings
set.seed(746)
fit_sann <- optim(c(20, 0.05), nll_gamma, method = "SANN")
fit_sann[c("par", "value", "counts", "convergence")]
#  [PAUSE] is SANN broken? should you trust convergence = 0?

# 4. Metropolis MCMC: annealing at T = 1, keep every step (~15 min) -----------------

# germination trial: 30 of 40 seeds; flat prior, so log posterior = log likelihood + const
k <- 30
n <- 40
log_post_germ <- function(p) {
  if (p <= 0 || p >= 1) return(-Inf)
  dbinom(k, n, p, log = TRUE) + dbeta(p, 1, 1, log = TRUE)
}

# Metropolis sampler: accept w.p. min(1, Post(new) / Post(old)); returns the whole chain
metropolis <- function(log_post, start, step_sd, n_steps) {
  chain <- matrix(NA, n_steps, length(start), dimnames = list(NULL, names(start)))
  current <- start
  lp_current <- log_post(current)
  n_accept <- 0
  for (i in seq_len(n_steps)) {
    proposal <- current + rnorm(length(start), 0, step_sd)
    lp_prop <- log_post(proposal)
    if (log(runif(1)) < lp_prop - lp_current) {   # log of the posterior ratio
      current <- proposal
      lp_current <- lp_prop
      n_accept <- n_accept + 1
    }
    chain[i, ] <- current   # a rejected step records the current value again
  }
  attr(chain, "accept_rate") <- n_accept / n_steps
  chain
}
#  [PAUSE] anneal() used runif(1) < exp(-delta / temp); this uses log(runif(1)) < ...; different rules?

set.seed(746)
chain_germ <- metropolis(log_post_germ, start = 0.5, step_sd = 0.1, n_steps = 5000)
attr(chain_germ, "accept_rate")

# trace plot, then histogram of the chain vs. the exact posterior, Beta(31, 11)
plot(chain_germ, type = "l", xlab = "Step", ylab = "p")
hist(chain_germ, breaks = 40, freq = FALSE, main = "", xlab = "Germination probability, p")
curve(dbeta(x, 31, 11), add = TRUE, lwd = 2, col = "darkorange")

# posterior summaries from the chain vs. exact
rbind(chain = c(mean = mean(chain_germ), quantile(chain_germ, c(0.025, 0.975))),
      exact = c(31 / 42, qbeta(c(0.025, 0.975), 31, 11)))

# two parameters: myxo Gamma with flat priors on shape and scale
set.seed(746)
chain_myxo <- metropolis(function(p) -nll_gamma(p),
                         start = c(shape = 20, scale = 0.05),
                         step_sd = c(4, 0.01), n_steps = 20000)
attr(chain_myxo, "accept_rate")

plot_surface()
points(chain_myxo[seq(1, 20000, by = 10), ], pch = 16, cex = 0.3,
       col = adjustcolor("steelblue", 0.4))

# 95% credible intervals, dropping the first 2000 steps
apply(chain_myxo[-(1:2000), ], 2, quantile, c(0.025, 0.5, 0.975))

# Not covered Oct. 19: the three-chains discussion moved to Wednesday (section 0) ----





# WEDNESDAY, Oct. 21: Hamiltonian Monte Carlo and Stan =========================
#  Uses titer, nll_gamma, plot_surface, log_post_germ, and metropolis() from Monday's section.

# 0. Lab 7 debrief: which chains to trust? (~3 min) ------------------------------

# germination data, same start (p = 0.05), 1000 steps, three proposal widths
step_sds <- c(0.005, 0.1, 3)
chains <- lapply(step_sds, function(s) {
  set.seed(1)
  metropolis(log_post_germ, start = 0.05, step_sd = s, n_steps = 1000)
})
round(rbind(step_sd = step_sds,
            accept_rate = sapply(chains, attr, "accept_rate"),
            post_mean = sapply(chains, mean)), 3)
#  [PAUSE] first is best, since nearly all were accepted? which Lab 7 chain is each like?

# reveal: trace plots
par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 0.5))
for (i in seq_along(chains)) {
  plot(chains[[i]], type = "l", ylim = c(0, 1), xlab = "Step", ylab = "p",
       main = paste("step sd =", step_sds[i]))
  abline(h = 31 / 42, lty = 2, col = "darkorange")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# 2. Why random walks struggle (~7 min) ----------------------------------------

# Monday's myxo chain: 20,000 steps, but how many independent draws' worth?
set.seed(746)
chain_myxo <- metropolis(function(p) -nll_gamma(p),
                         start = c(shape = 20, scale = 0.05),
                         step_sd = c(4, 0.01), n_steps = 20000)
keep <- chain_myxo[-(1:2000), ]
coda::effectiveSize(keep)

# a narrow, tilted target: bivariate normal with correlation 0.95
rho <- 0.95
sigma_inv <- solve(matrix(c(1, rho, rho, 1), 2))
log_post_bvn <- function(p) -0.5 * drop(t(p) %*% sigma_inv %*% p)
grad_bvn <- function(p) -drop(sigma_inv %*% p)   # gradient of the log posterior

# contours of the target, for the background of both plots
plot_bvn <- function(main) {
  g <- seq(-3, 3, length.out = 100)
  z <- outer(g, g, Vectorize(function(a, b) log_post_bvn(c(a, b))))
  contour(g, g, z, levels = c(-0.5, -2, -4.5), drawlabels = FALSE, col = "grey60",
          xlab = expression(theta[1]), ylab = expression(theta[2]), main = main)
}

set.seed(3)
rw <- metropolis(log_post_bvn, start = c(-2, -2), step_sd = 0.25, n_steps = 200)
plot_bvn("Random-walk Metropolis, 200 steps")
lines(rw, col = "steelblue")
coda::effectiveSize(rw[, 1])

# 3. Hamiltonian Monte Carlo (~10 min) -----------------------------------------

# one HMC step: random momentum kick, L leapfrog steps along the gradient, then accept/reject
hmc_step <- function(p, log_post, grad, eps, L) {
  m <- rnorm(length(p))
  p_new <- p
  m_new <- m
  path <- matrix(p, nrow = 1)
  for (i in seq_len(L)) {
    m_new <- m_new + 0.5 * eps * grad(p_new)   # half step for momentum
    p_new <- p_new + eps * m_new                # full step for position
    m_new <- m_new + 0.5 * eps * grad(p_new)   # half step for momentum
    path <- rbind(path, p_new)
  }
  # total energy = -log posterior + kinetic energy m^2 / 2; accept w.p. exp(-change)
  log_accept <- log_post(p_new) - 0.5 * sum(m_new^2) - log_post(p) + 0.5 * sum(m^2)
  list(p = if (log(runif(1)) < log_accept) p_new else p, path = path)
}

# run n HMC steps; keep the first few trajectories for plotting
hmc <- function(start, log_post, grad, eps, L, n_steps, n_paths = 5) {
  chain <- matrix(NA, n_steps, length(start))
  paths <- list()
  p <- start
  for (i in seq_len(n_steps)) {
    s <- hmc_step(p, log_post, grad, eps, L)
    p <- s$p
    chain[i, ] <- p
    if (i <= n_paths) paths[[i]] <- s$path
  }
  list(chain = chain, paths = paths)
}

set.seed(3)
hm <- hmc(c(-2, -2), log_post_bvn, grad_bvn, eps = 0.1, L = 10, n_steps = 200)

# side by side: random walk vs. HMC, same number of steps
par(mfrow = c(1, 2))
plot_bvn("Random walk")
lines(rw, col = "steelblue")
plot_bvn("HMC")
points(hm$chain, pch = 16, cex = 0.5)
for (pa in hm$paths) lines(pa, col = "darkorange", lwd = 2)
par(mfrow = c(1, 1))

c(random_walk = coda::effectiveSize(rw[, 1]), hmc = coda::effectiveSize(hm$chain[, 1]))
#  [PAUSE] each HMC step costs 10 gradient evaluations; is it still worth it?

# 4. Stan: the myxo model (~12 min) --------------------------------------------

library(cmdstanr)

# same model as Monday: Gamma(shape, scale), flat priors on the positive values
myxo_code <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] titer;
}
parameters {
  real<lower=0> shape;
  real<lower=0> scale;
}
model {
  titer ~ gamma(shape, 1 / scale);  // Stan's gamma uses shape and rate
}
"
myxo_model <- cmdstan_model(write_stan_file(myxo_code))   # compiles to C++ (slow, once)

fit_stan <- myxo_model$sample(data = list(N = length(titer), titer = titer),
                              chains = 4, parallel_chains = 4, seed = 746)

# R-hat (want about 1.00) and effective sample sizes
fit_stan$summary()
fit_stan$diagnostic_summary()

# Stan's 4,000 draws on Monday's surface
draws <- fit_stan$draws(c("shape", "scale"), format = "draws_matrix")
plot_surface()
points(draws[, "shape"], draws[, "scale"], pch = 16, cex = 0.3,
       col = adjustcolor("darkorange", 0.4))

# effective draws: Monday's 18,000 kept steps vs. Stan's 4,000 draws
rbind(metropolis = coda::effectiveSize(keep),
      stan = fit_stan$summary(c("shape", "scale"))$ess_bulk)

# 95% credible intervals for shape
rbind(metropolis = quantile(keep[, 1], c(0.025, 0.975)),
      stan = quantile(draws[, "shape"], c(0.025, 0.975)))

# Discussion (~5 min): is Stan worth it? ----------------------------------------
#  [DISCUSS] numbers are in the NOTES prompt
