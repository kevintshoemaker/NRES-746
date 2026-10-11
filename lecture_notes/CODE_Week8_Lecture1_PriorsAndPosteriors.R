
#  NRES 746, Week 8, Lecture 1: Priors and posteriors -------------------------
#  Section numbers match the lecture notes. No extra packages needed.

# WOD debrief: tortoise analysts (~8 min) --------------------------------------
#  [BOARD] collect likelihood sketches -> [CODE] reveal

survived <- 21
collared <- 30
wod_priors <- list(A = c(2, 2), B = c(46, 4), C = c(10, 10))
s_seq <- seq(0.001, 0.999, length.out = 500)

par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 0.5))
for (an in names(wod_priors)) {
  pr <- wod_priors[[an]]
  po <- pr + c(survived, collared - survived)
  plot(s_seq, dbeta(s_seq, po[1], po[2]), type = "l", lwd = 2, ylim = c(0, 12),
       xlab = "Annual survival probability", ylab = "Density",
       main = paste0("Analyst ", an, ": posterior mean ", round(po[1] / sum(po), 2)))
  lines(s_seq, dbeta(s_seq, pr[1], pr[2]), lwd = 2, col = "grey60")
  lines(s_seq, dbeta(s_seq, survived + 1, collared - survived + 1),
        lwd = 2, lty = 2, col = "darkorange")   # the reveal: shared likelihood
}
par(mfrow = c(1, 1))

# 1. Big picture and recap (~3 min) -------------------------------------------
#  [BOARD] workflow loop; Beta + binomial -> Beta

# 2. Germination, two priors (~10 min) -----------------------------------------

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
        info_mean = beta_mean(post_info)), 3)

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

beta_mean(prior_info + c(300, 100))   # more data: 300 of 400, same prior

post_two_steps <- (prior_info + c(k, n - k)) + c(k, n - k)   # 40 seeds, then 40 more
post_one_step <- prior_info + c(2 * k, 2 * (n - k))         # 80 seeds at once
rbind(post_two_steps, post_one_step)

# 3. Summarizing a posterior (~4 min) ------------------------------------------
#  [CODE]

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

# 4. Where new Bayesians get stuck (~8 min) ------------------------------------
#  [BOARD] list (a)-(f); [CODE] for (a), (b), (e)

# (a) area under the nest likelihood is not 1
integrate(function(p) p^3 * (1 - p), 0, 1)

# (b) uniform on p -> peaked prior on the log-odds, log(p / (1 - p))
set.seed(746)
p_draws <- runif(1e5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.5))
hist(p_draws, breaks = 40, freq = FALSE, col = "grey80", border = "white",
     main = "Flat prior on p", xlab = "p")
hist(qlogis(p_draws), breaks = 60, freq = FALSE, col = "grey80",
     border = "white", main = "Same prior, log-odds scale",
     xlab = "log(p / (1 - p))", xlim = c(-8, 8))
curve(dlogis(x), add = TRUE, lwd = 2)
par(mfrow = c(1, 1))

# (e) Week 7 pikas: LR interval vs. flat-prior credible interval
haypiles <- c(2, 5, 3, 0, 4, 4)
nll_pois <- function(lambda) -sum(dpois(haypiles, lambda, log = TRUE))
f_cut <- function(l) nll_pois(l) - nll_pois(mean(haypiles)) - qchisq(0.95, 1) / 2
lr_int <- c(uniroot(f_cut, c(0.5, 3))$root, uniroot(f_cut, c(3, 8))$root)
cred_int <- qgamma(c(0.025, 0.975), shape = sum(haypiles) + 1,      
                   rate = length(haypiles))   # posterior Gamma(19, 6)  NOTE: flat, improper prior
round(rbind(likelihood_ratio = lr_int, bayes_flat_prior = cred_int), 2)

# 5. Why Bayesian methods are so widely used (~5 min) --------------------------
#  [BOARD] reasons and costs -> [CODE] derived quantities

# derived quantity: posterior for the odds of germination, p / (1 - p)
odds_draws <- draws_flat / (1 - draws_flat)
round(quantile(odds_draws, c(0.025, 0.5, 0.975)), 2)
hist(odds_draws, breaks = 60, freq = FALSE, col = "grey80", border = "white",
     main = "Posterior for the odds", xlab = "p / (1 - p)")

# probability statements come straight from the draws
mean(draws_flat > 0.6)   # P(germination rate > 0.6 | data, flat prior)

# 6. Return to the big picture (~1 min) ----------------------------------------
#  [BOARD] workflow loop; hand out Newton take-home

# Discussion (~8 min) ----------------------------------------------------------
#  [DISCUSS]
