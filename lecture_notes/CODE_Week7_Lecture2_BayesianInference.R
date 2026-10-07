
#  NRES 746, Week 7, Lecture 2: Bayesian inference ---------------------------
#  Logistic regression; prior x likelihood = posterior, from Bayes' rule.
#  Run section by section in RStudio (use the outline to jump around).
#  Section numbers match the lecture notes. No extra packages needed.

# Opening: logistic regression and the log-odds ------------------------------
#  p = 1 / (1 + exp(-(a + b x)))   <->   log(p / (1 - p)) = a + b x
#  Odds = p / (1 - p): p = 0.75 is odds of 3 ("3 to 1")
#  b = change in log-odds per unit x; exp(b) multiplies the odds
#  p = 0.5 where a + b x = 0, i.e. x = -a / b (the LD50 in Lab 6)

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

# 1. Big picture --------------------------------------------------------------
#  THEORY -> MODEL(S) -> DATA (fit) -> COMPARE -> repeat
#  Same model, same likelihood as Monday; new way to describe uncertainty.
#  Frequentist: parameter fixed, data random, CI = a procedure.
#  Bayesian: uncertainty about the parameter = a probability distribution.
#  prior -> (likelihood) -> posterior, via Bayes' rule (Week 3, Lab 3).

# 2a. Bayes' rule with three hypotheses ---------------------------------------
#  P(H_i | D) = P(D | H_i) P(H_i) / sum_j P(D | H_j) P(H_j)
#  Week 5 WOD nests: outcomes 1, 1, 0, 1, so L(p) = p^3 (1 - p)

p3 <- c(0.25, 0.5, 0.75)
lik3 <- p3^3 * (1 - p3)
prior3 <- rep(1/3, 3)
post3 <- lik3 * prior3 / sum(lik3 * prior3)
round(rbind(likelihood = lik3, posterior = post3), 3)

sum(lik3)    # 46/256: likelihoods don't sum to 1
lik3 * 256   # 3, 16, 27
post3 * 46   # posterior = 3/46, 16/46, 27/46

# 2b. From a few hypotheses to a whole distribution ---------------------------
#  Same calculation with 3, 11, 101 values of p; flat prior on every grid
#  Why do the bars shrink? One value of p has probability 0
#  Bar / spacing -> density (Beta(4, 2)); areas give probabilities
#  P(p | y) = L(p) P(p) / integral of L(p) P(p) dp
#  Denominator = one number (marginal likelihood); hard in many dimensions -> MCMC
#  Sticking points: posterior = our knowledge of p, not p varying;
#  a density is not a probability (it can exceed 1)

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

# 3. The Beta prior and conjugacy --------------------------------------------
#  Binomial likelihood: L(p) ~ p^k (1 - p)^(n - k)
#  Beta(alpha, beta) prior: P(p) ~ p^(alpha - 1) (1 - p)^(beta - 1)
#  Numerator: p^(k + alpha - 1) (1 - p)^(n - k + beta - 1) -> Beta(alpha + k, beta + n - k)
#  Conjugate = just add counts; alpha - 1, beta - 1 = prior successes, failures
#  Nests (3 of 4) with a flat Beta(1, 1) prior -> Beta(4, 2)

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

# 4. Return to the big picture ------------------------------------------------
#  Same model, same likelihood; multiply by a prior, read off the posterior
#  Conjugate = add counts
#  Next time: germination example with two priors (Week 8, Lecture 1)
