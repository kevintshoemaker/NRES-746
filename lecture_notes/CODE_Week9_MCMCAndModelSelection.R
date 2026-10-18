
#  NRES 746, Week 9: Optimization to MCMC (Mon.); model selection (Wed.) ======
#  R demos only; section numbers match the NOTES files. Needs emdbook.

# MONDAY, Oct. 19: Nelder-Mead, simulated annealing, Metropolis ===============

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
nelder_mead <- function(f, start, step, n_iter = 80) {
  simplex <- rbind(start, start + c(step[1], 0), start + c(0, step[2]))
  vals <- apply(simplex, 1, f)
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

nm <- nelder_mead(nll_gamma, c(20, 0.05), step = c(5, 0.02))
round(c(shape = nm$best[1], scale = nm$best[2], nll = nm$value), 4)
table(nm$moves)

# first 30 triangles, colored by move
move_cols <- c(expand = "darkorange", reflect = "steelblue",
               contract = "firebrick", shrink = "purple")
plot_surface()
for (i in 1:30) polygon(nm$history[[i]], border = move_cols[nm$moves[i]], lwd = 1.5)
legend("topright", names(move_cols), col = move_cols, lwd = 2, bty = "n")

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

# path on the surface (left) and current NLL as the temperature falls (right)
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

# 4. Metropolis: annealing at T = 1, keep every step (~15 min) -----------------

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

# Discussion (~6 min): three germination chains --------------------------------

# same data, same start (p = 0.05), 1000 steps, three proposal widths
step_sds <- c(0.005, 0.1, 3)
chains <- lapply(step_sds, function(s) {
  set.seed(1)
  metropolis(log_post_germ, start = 0.05, step_sd = s, n_steps = 1000)
})
round(rbind(step_sd = step_sds,
            accept_rate = sapply(chains, attr, "accept_rate"),
            post_mean = sapply(chains, mean)), 3)
#  [PAUSE] which chain would you trust? exact posterior mean is 0.738

# reveal: trace plots
par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 0.5))
for (i in seq_along(chains)) {
  plot(chains[[i]], type = "l", ylim = c(0, 1), xlab = "Step", ylab = "p",
       main = paste("step sd =", step_sds[i]))
  abline(h = 31 / 42, lty = 2, col = "darkorange")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# WEDNESDAY, Oct. 21: Model selection ==========================================
#  (to be added)
