#  NRES 746, Week 4: Counting, expectation, and probability distributions -----
#  Code from Week4_ProbabilityDistributions.html. Base R only.

# COUNTING, RANDOM VARIABLES, AND EXPECTATION ===================================

# 2. WOD: a probability puzzle -------------------------------------------------
#  two dice: 36 ordered pairs; sums of 7 vs. doubles

sums <- outer(1:6, 1:6, "+")
sum(sums == 7) / 36
sum(diag(6)) / 36
#  [PAUSE] two eight-sided dice: which sum is most likely, and does it tie, beat, or lose to doubles?

# 3. Counting arrangements ------------------------------------------------------
#  permutations P(N, k) = N! / (N - k)!; combinations choose(N, k) = P(N, k) / k!

12 * 11 * 10 * 9          # ordered visits to 4 of 12 ponds
choose(12, 4)             # sets of 4 ponds
choose(30, 5)             # lottery tickets from the WOD

#  binomial: P(X = k) = choose(N, k) p^k (1 - p)^(N - k); 6 nests, p = 0.25
combn(6, 2)                                  # the 15 ways 2 of 6 nests can fail
choose(6, 2) * 0.25^2 * 0.75^4
dbinom(2, size = 6, prob = 0.25)

# 4. The birthday problem -------------------------------------------------------
#  P(all different) = prod_{i = 0}^{n - 1} (365 - i) / 365

n_people <- 25
p_all_diff <- prod((365 - 0:(n_people - 1)) / 365)
1 - p_all_diff
choose(n_people, 2)
1 - exp(-choose(n_people, 2) / 365)          # pairs approximation
#  [PAUSE] smallest group with better than even odds, for 365 and for 10,000 values? pairs?

# 5. Random variables and their distributions ----------------------------------

x_nest <- 0:6
p_nest <- dbinom(x_nest, size = 6, prob = 0.25)
round(p_nest, 3)
sum(p_nest)
barplot(p_nest, names.arg = x_nest, xlab = "Failed nests (of 6)", ylab = "Probability")

#  densities can exceed 1: uniform pH on [7.0, 7.1] has density 1 / 0.1
dunif(7.05, min = 7.0, max = 7.1)

#  trout lengths ~ Normal(22, 4): P(20 <= X <= 24) = F(24) - F(20)
pnorm(24, mean = 22, sd = 4) - pnorm(20, mean = 22, sd = 4)
dnorm(22, mean = 22, sd = 4)                  # a density, per cm
#  [PAUSE] widget: shrink the SD to 0.1. is the model broken? narrow the interval: area vs. height x width?

# 6. Expectation and variance ---------------------------------------------------
#  mean of tabulated data = sum_x x c(x) / N; replace c(x) / N with P(X = x) to get E[X]

set.seed(4)
routes <- rbinom(40, size = 6, prob = 0.25)
counts <- table(factor(routes, levels = 0:6))
counts
mean(routes)
sum(0:6 * counts / 40)

sum(x_nest * p_nest)                                      # binomial: N p = 1.5
integrate(function(x) x * dexp(x, rate = 1.5), 0, Inf)    # exponential: 1 / rate

#  Var(X) = E[(X - mu)^2] = E[X^2] - mu^2
mu_nest <- sum(x_nest * p_nest)
sum((x_nest - mu_nest)^2 * p_nest)
sum(x_nest^2 * p_nest) - mu_nest^2
6 * 0.25 * 0.75                                           # N p (1 - p)
#  [PAUSE] no survey can record 1.5 failures. in what sense is 1.5 "expected"?

# Discussion ---------------------------------------------------------------------

dnorm(7.9, mean = 7.9, sd = 0.02)            # pH density
#  [PAUSE] is a density of 19.9 a problem?
#  [PAUSE] 150 bears, chance match 1 in 10,000 per pair: how worried should the study be?

# Bonus: computing on the log scale ---------------------------------------------
#  P(all different) = 365! / ((365 - n)! 365^n)

factorial(365)                               # overflows to Inf
365^25
log_p <- lfactorial(365) - lfactorial(365 - 25) - 25 * log(365)
exp(log_p)
p_all_diff

# PROBABILITY DISTRIBUTIONS =====================================================

# 1. WOD: expectation and the sample mean ---------------------------------------
#  E[xbar] = mu, Var(xbar) = sigma^2 / n (independent draws), SE = sigma / sqrt(n)

set.seed(746)
xbar_5 <- replicate(10000, mean(rnorm(5, mean = 22, sd = 4)))
var(xbar_5)                                  # theory: 16 / 5 = 3.2
sd(xbar_5)                                   # theory: 4 / sqrt(5)
hist(xbar_5, breaks = 50, xlab = "Mean length of 5 trout (cm)", main = "")
#  [PAUSE] five trout from one pool, not independent: is Var(xbar) larger or smaller than sigma^2 / 5?

# 2. Four functions for every distribution --------------------------------------
#  d = pmf/pdf, p = P(X <= x), q = inverse of p, r = random draws

dnorm(22, mean = 22, sd = 4)                 # density at 22 cm
1 - pnorm(30, mean = 22, sd = 4)             # P(longer than 30 cm)
qnorm(0.9, mean = 22, sd = 4)                # 90% of trout are shorter than this
set.seed(1)
rnorm(5, mean = 22, sd = 4)                  # five simulated trout

pbinom(2, size = 6, prob = 0.25)             # P(X <= 2) includes 2
1 - pbinom(2, size = 6, prob = 0.25)         # P(X >= 3)

# 3. Discrete distributions -----------------------------------------------------
#  Poisson: mean = variance = lambda; NegBin(mu, k): variance = mu + mu^2 / k

x_seed <- 0:15
par(mfrow = c(1, 2))
barplot(dpois(x_seed, 4), names.arg = x_seed, main = "Poisson(4)",
        xlab = "Seeds per quadrat", ylab = "Probability", ylim = c(0, 0.22))
barplot(dnbinom(x_seed, mu = 4, size = 1.5), names.arg = x_seed, main = "NegBin(4, k = 1.5)",
        xlab = "Seeds per quadrat", ylab = "Probability", ylim = c(0, 0.22))
par(mfrow = c(1, 1))

c(poisson = 4, negbin = 4 + 4^2 / 1.5)                             # variances
c(poisson = dpois(0, 4), negbin = dnbinom(0, mu = 4, size = 1.5))  # P(no seeds)

# 4. Continuous distributions ---------------------------------------------------
#  gamma (shape a, scale s): mean a s, variance a s^2; beta(a, b): mean a / (a + b)

par(mfrow = c(1, 4))
curve(dgamma(x, shape = 2, scale = 2), 0, 15, lwd = 2, xlab = "Biomass (g)",
      ylab = "Density", main = "Gamma(a = 2, s = 2)")
curve(dexp(x, rate = 1.5), 0, 4, lwd = 2, xlab = "Hours between visits",
      ylab = "Density", main = "Exponential(1.5)")
curve(dlnorm(x, meanlog = 2, sdlog = 0.5), 0, 25, lwd = 2, xlab = "Body mass (g)",
      ylab = "Density", main = "Lognormal(2, 0.5)")
curve(dbeta(x, 2, 5), 0, 1, lwd = 2, xlab = "Proportion cover",
      ylab = "Density", main = "Beta(2, 5)")
par(mfrow = c(1, 1))
#  [PAUSE] widget: which distribution fits counts with variance 3.5 times the mean? what k?

# 5. Choosing a distribution ----------------------------------------------------
#  method of moments for the negative binomial: k = xbar / (s^2 / xbar - 1)

set.seed(14)
seeds <- rnbinom(20, mu = 4, size = 1.5)
seeds
c(mean = mean(seeds), variance = var(seeds), ratio = var(seeds) / mean(seeds))

m <- mean(seeds)
k_hat <- m / (var(seeds) / m - 1)
k_hat

x_q <- 0:16
mids <- barplot(table(factor(seeds, levels = x_q)), ylim = c(0, 5),
                xlab = "Seeds per quadrat", ylab = "Number of quadrats")
points(mids, 20 * dpois(x_q, m), pch = 16, type = "b", lty = 2)
points(mids, 20 * dnbinom(x_q, mu = m, size = k_hat), pch = 17, type = "b", col = "firebrick")

# Discussion ---------------------------------------------------------------------
#  [PAUSE] clutch size always 2 to 5 eggs, modeled as Poisson(3.4): what does it get wrong?
#  [PAUSE] seed counts with variance 3.5 times the mean: which processes could cause it?
#  [PAUSE] fraction of leaf area eaten: beta or binomial? what must each student believe?

# Bonus: where the negative binomial comes from ---------------------------------
#  lambda_i ~ Gamma(k, mu / k), Y_i ~ Pois(lambda_i)  ->  Y_i ~ NegBin(mu, k)

set.seed(5)
rates <- rgamma(10000, shape = 1.5, scale = 4 / 1.5)
y_mix <- rpois(10000, rates)
c(mean = mean(y_mix), variance = var(y_mix))              # NegBin: 4 and 4 + 4^2 / 1.5

x_c <- 0:20
mids <- barplot(table(factor(y_mix, levels = x_c)) / 10000, ylim = c(0, 0.17),
                xlab = "Seeds per quadrat", ylab = "Proportion")
points(mids, dnbinom(x_c, mu = 4, size = 1.5), pch = 17, col = "firebrick")
