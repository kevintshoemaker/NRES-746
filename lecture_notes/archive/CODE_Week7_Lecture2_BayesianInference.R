
#  NRES 746, Week 7, Lecture 2: Bayesian inference ---------------------------
#  Section numbers match the lecture notes. No extra packages needed.

# Opening: logistic regression and the log-odds (~10 min) ---------------------
#  [BOARD] logistic curve <-> log-odds line -> [CODE]

a_ex <- -2   # nest success vs. distance from edge (100s of m)
b_ex <- 0.8
x_seq <- seq(0, 6, length.out = 200)

par(mfrow = c(1, 2), mar = c(4, 4.5, 2, 0.5))
plot(x_seq, plogis(a_ex + b_ex * x_seq), type = "l", lwd = 2, ylim = c(0, 1),
     xlab = "Distance from edge (100 m)", ylab = "P(nest succeeds)",
     main = "Probability scale")
abline(h = 0.5, v = -a_ex / b_ex, lty = 3)
plot(x_seq, a_ex + b_ex * x_seq, type = "l", lwd = 2,
     xlab = "Distance from edge (100 m)", ylab = "log(p / (1 - p))",
     main = "Log-odds scale")
abline(h = 0, v = -a_ex / b_ex, lty = 3)
par(mfrow = c(1, 1))

qlogis(0.75)            # log-odds of p = 0.75 is log(3)
plogis(qlogis(0.75))    # plogis() undoes qlogis()

# 1. Big picture (~3 min) ------------------------------------------------------
#  [BOARD] workflow loop; frequentist vs. Bayesian uncertainty

# 2a. Bayes' rule with three hypotheses ---------------------------------------
#  [BOARD] Bayes' rule, three hypotheses -> [CODE]

p3 <- c(0.25, 0.5, 0.75)
lik3 <- p3^3 * (1 - p3)
prior3 <- rep(1/3, 3)
post3 <- lik3 * prior3 / sum(lik3 * prior3)
round(rbind(likelihood = lik3, posterior = post3), 3)

sum(lik3)    # 46/256: likelihoods don't sum to 1
lik3 * 256   # 3, 16, 27
post3 * 46   # posterior = 3/46, 16/46, 27/46

# 2b. From a few hypotheses to a whole distribution ---------------------------
#  [BOARD] bars -> density; denominator -> [CODE]

par(mfrow = c(1, 3), mar = c(4, 4, 2, 2.5))
for (n_grid in c(3, 11, 101)) {
  p_g <- seq(0, 1, length.out = n_grid + 2)[-c(1, n_grid + 2)]
  w <- p_g^3 * (1 - p_g)
  prob <- w / sum(w)
  spacing <- diff(p_g)[1]
  plot(p_g, prob, type = "h", lwd = if (n_grid < 20) 4 else 1,
       col = if (n_grid < 20) "black" else "grey60",
       xlim = c(0, 1), xlab = "p", ylab = "Posterior probability",
       main = paste(n_grid, "hypotheses"))
  if (n_grid == 101) {
    lines(p_g, dbeta(p_g, 4, 2) * spacing, lwd = 2, col = "darkorange")
    axis(4, at = (0:2) * spacing, labels = 0:2, col.axis = "darkorange")
    mtext("density", side = 4, line = 1.6, cex = 0.7, col = "darkorange")
  }
  cat(n_grid, "values: tallest bar", round(max(prob), 3),
      "; bar / spacing", round(max(prob) / spacing, 2), "\n")
}
par(mfrow = c(1, 1))
dbeta(0.75, 4, 2)   # peak of the limiting density

# 3. The Beta prior and conjugacy (~8 min) -------------------------------------
#  [BOARD] conjugacy algebra: add counts -> [CODE] check numerically

k_nest <- 3
n_nest <- 4
p_seq <- seq(0.001, 0.999, length.out = 500)
numer <- p_seq^k_nest * (1 - p_seq)^(n_nest - k_nest) * dbeta(p_seq, 1, 1)
step <- diff(p_seq)[1]
plot(p_seq, numer / sum(numer * step), type = "l", lwd = 4, col = "grey70",
     xlab = "p", ylab = "Posterior density")   # prior x likelihood, rescaled
lines(p_seq, dbeta(p_seq, 1 + k_nest, 1 + n_nest - k_nest), lwd = 2, lty = 2)
legend("topleft", c("prior x likelihood, rescaled", "Beta(4, 2)"),
       lwd = c(4, 2), lty = c(1, 2), col = c("grey70", "black"), bty = "n")

# 4. Wrap (~1 min) -------------------------------------------------------------
#  [BOARD] workflow loop
