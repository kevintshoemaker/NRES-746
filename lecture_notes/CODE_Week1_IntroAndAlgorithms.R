#  NRES 746, Week 1: Course introduction and algorithms for inference -----------
#  Code from Week1_IntroAndAlgorithms.html. Base R only.

# WHAT THIS COURSE IS ABOUT =====================================================

# 3. WOD: a self-diagnostic -----------------------------------------------------

#  N_t = N_0 exp(r t)  ->  r = ln(N_t / N_0) / t
log(550 / 100) / 8

#  chain rule: d/dx (4x^2 - 3x)^5 = 5(4x^2 - 3x)^4 (8x - 3)
fprime <- D(quote((4 * x^2 - 3 * x)^5), "x")
fprime
eval(fprime, list(x = 1))

#  R(t) = 2t^2, so the total from 0 to 5 is R(5) - R(0) = 50
integrate(function(t) 4 * t, lower = 0, upper = 5)

#  P(both girls | at least one girl) = 1/3
set.seed(1)
kids <- matrix(sample(c("G", "B"), 2e5, replace = TRUE), ncol = 2)
at_least_one <- kids[, 1] == "G" | kids[, 2] == "G"
both <- kids[, 1] == "G" & kids[, 2] == "G"
mean(both[at_least_one])

#  the three R bugs
site_data <- matrix(c(12, 7, 19, 4, 45, 23, 61, 15), nrow = 4, ncol = 2)
site_data[2]          # the 2nd element, counting down the columns
site_data[, 2]        # the whole 2nd column

values <- c(10, 20, 30, 40)
doubled <- numeric(length(values))
for (i in seq_along(values)) {     # loop over positions 1 to 4, not values 10 to 40
  doubled[i] <- values[i] * 2
}
doubled

get_range <- function(x) {         # braces keep all three lines inside the function
  max_val <- max(x)
  min_val <- min(x)
  max_val - min_val
}
get_range(values)

#  [PAUSE] A 95% CI for nest survival is 0.42 to 0.61. Is "a 95% chance the true value is in there" wrong?


# MATH TOOLS FOR LAB 1 ==========================================================

# 2. Exponents and logarithms ---------------------------------------------------

log10(2) + log10(50)
log10(500) - log10(5)
c(log10(3^4), 4 * log10(3))
log10(7e6)                         # order of magnitude: in the millions

#  power law B = a M^b  ->  log10(B) = log10(a) + b log10(M)
mass <- c(1, 10, 100, 1000)
rate <- 3 * mass^0.75
data.frame(mass, rate = round(rate, 1), log10_rate = round(log10(rate), 3))

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
curve(3 * x^0.75, from = 1, to = 1000, lwd = 2, xlab = "Body mass, M",
      ylab = "Metabolic rate, B", main = "Original scale")
curve(3 * x^0.75, from = 1, to = 1000, lwd = 2, log = "xy",
      xlab = "Body mass, M (log scale)", ylab = "B (log scale)",
      main = "Log-log scale")
par(mfrow = c(1, 1))
#  [PAUSE] How much does log10(B) rise per tenfold increase in mass? What would a slope of 1 mean?

# 3. Derivatives: the power, product, and chain rules ---------------------------

#  f = x^2 e^-x; product + chain rules give f' = x(2 - x)e^-x, zero at x = 0 and 2
D(quote(x^2 * exp(-x)), "x")
f <- function(x) x^2 * exp(-x)
optimize(f, interval = c(0, 10), maximum = TRUE)$maximum

curve(f(x), from = 0, to = 8, lwd = 2, ylim = c(-0.2, 0.6),
      xlab = "x", ylab = "")
curve(x * (2 - x) * exp(-x), add = TRUE, lty = 2, lwd = 2, col = "steelblue")
abline(h = 0, col = "grey60")
abline(v = 2, lty = 3)
legend("topright", c("f(x)", "f'(x)"), lty = 1:2, lwd = 2,
       col = c("black", "steelblue"), bty = "n")

# 4. Adding up a rate: a first look at integration ------------------------------

#  ramp from 0 to 4 per day over 5 days: area = (1/2)(5)(4) = 10; F(t) = 0.4 t^2
integrate(function(t) 0.8 * t, lower = 0, upper = 5)
F <- function(t) 0.4 * t^2
F(5) - F(0)


# ALGORITHMS FOR INFERENCE ======================================================

# 2. WOD: pseudocoding a permutation test ---------------------------------------

lizards <- data.frame(
  mass = c(175, 168, 168, 190, 156, 181, 182, 175, 174, 179,
           185, 169, 173, 173, 188, 186, 175, 174, 179, 180),
  group = rep(c("trtA", "Control"), each = 10)
)
group_means <- tapply(lizards$mass, lizards$group, mean)
group_means
observed_dif <- group_means[["trtA"]] - group_means[["Control"]]
observed_dif

# 3. Two routes to a z-test -----------------------------------------------------

pop_mean <- 4.5
pop_sd <- 0.9
salmon <- c(3.14, 3.27, 2.56, 3.77, 3.34, 4.32, 3.84, 2.19, 5.24, 3.09)
n <- length(salmon)
sample_mean <- mean(salmon)
sample_mean

hist(salmon, freq = FALSE, breaks = seq(1.5, 7, by = 0.5), col = "grey85",
     xlim = c(1.5, 7), ylim = c(0, 0.6), main = "",
     xlab = "Body mass (kg)")
curve(dnorm(x, pop_mean, pop_sd), add = TRUE, lwd = 2, col = "firebrick")
abline(v = sample_mean, lwd = 3, col = "steelblue")

#  route 1, analytical: SE = sigma / sqrt(n)
pop_se <- pop_sd / sqrt(n)
z <- (sample_mean - pop_mean) / pop_se
p_analytical <- pnorm(sample_mean, mean = pop_mean, sd = pop_se)
c(z = z, p = p_analytical)

#  route 2, computational: simulate samples of 10 fish with no diet effect
set.seed(746)
n_sims <- 100000
null_means <- numeric(n_sims)
for (i in seq_len(n_sims)) {
  fake_sample <- rnorm(n, mean = pop_mean, sd = pop_sd)
  null_means[i] <- mean(fake_sample)
}
p_simulated <- mean(null_means <= sample_mean)
p_simulated
sd(null_means)                     # compare with pop_se

hist(null_means, breaks = 60, col = "grey85", main = "",
     xlab = "Mean of 10 fish, simulated with no diet effect (kg)",
     xlim = range(c(null_means, sample_mean)))
abline(v = sample_mean, lwd = 3, col = "steelblue")

# 4. Shuffling labels: a permutation test ---------------------------------------

set.seed(746)
reps <- 10000
null_difs <- numeric(reps)
for (i in seq_len(reps)) {
  shuffled <- sample(lizards$group)
  null_difs[i] <- mean(lizards$mass[shuffled == "trtA"]) -
    mean(lizards$mass[shuffled == "Control"])
}
p_perm <- mean(abs(null_difs) >= abs(observed_dif))
p_perm

hist(null_difs, breaks = 40, col = "grey85", main = "",
     xlab = "trtA mean minus Control mean, labels shuffled (g)")
abline(v = c(observed_dif, -observed_dif), lwd = 3, col = "steelblue",
       lty = c(1, 2))
#  [PAUSE] Many runs of 1,000 salmon replicates contain no mean as low as 3.476. What p-value would you report?

#  the same loop as a function
permutation_test <- function(value, group, reps = 10000) {
  labels <- unique(group)
  dif <- function(g) mean(value[g == labels[1]]) - mean(value[g == labels[2]])
  observed <- dif(group)
  null_difs <- numeric(reps)
  for (i in seq_len(reps)) {
    null_difs[i] <- dif(sample(group))
  }
  list(observed = observed, p_value = mean(abs(null_difs) >= abs(observed)))
}

set.seed(746)
permutation_test(lizards$mass, lizards$group)
t.test(mass ~ group, data = lizards, var.equal = TRUE)$p.value

# Bonus: bootstrapping a confidence interval ------------------------------------

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(Volume ~ Girth, data = trees, pch = 16, col = "firebrick")
plot(Volume ~ Height, data = trees, pch = 16, col = "steelblue")
par(mfrow = c(1, 1))

#  R^2 of volume on each predictor, one at a time
r_squared <- function(data, response = "Volume") {
  predictors <- setdiff(names(data), response)
  sapply(predictors, function(p) {
    summary(lm(reformulate(p, response), data = data))$r.squared
  })
}
r_squared(trees)

#  resample trees with replacement, recompute R^2, take quantiles
set.seed(746)
boot_rsq <- t(replicate(2000, {
  rows <- sample.int(nrow(trees), replace = TRUE)
  r_squared(trees[rows, ])
}))
apply(boot_rsq, 2, quantile, probs = c(0.025, 0.975))
