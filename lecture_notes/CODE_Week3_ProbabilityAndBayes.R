#  NRES 746, Week 3: Probability rules and Bayes' rule -------------------------
#  Code from Week3_ProbabilityAndBayes.html. Base R only.

# ADDING AND MULTIPLYING PROBABILITIES ==========================================

# 2. Sample spaces and events ---------------------------------------------------
#  four mutually exclusive seed fates; P(A or B) = P(A) + P(B); P(not A) = 1 - P(A)

fates <- c(rodent = 0.35, bird = 0.15, germinate = 0.10, soil = 0.40)
sum(fates)                              # rule 2
fates["rodent"] + fates["bird"]         # P(eaten), by rule 3
1 - fates["germinate"]                  # P(does not germinate)

# 3. Multiplying probabilities --------------------------------------------------
#  expected count = n p; independent events: P(A and B) = P(A) P(B)

200 * 0.35                              # expected rodent-eaten seeds of 200

p_find <- 0.6
outcomes <- c(both     = p_find * p_find,
              only_1   = p_find * (1 - p_find),
              only_2   = (1 - p_find) * p_find,
              neither  = (1 - p_find) * (1 - p_find))
outcomes
sum(outcomes)
1 - outcomes["neither"]                 # P(at least one observer finds it)
p_find + p_find - p_find^2              # general addition rule, same answer
#  [PAUSE] if hidden nests are hard for both observers, is P(both miss) more or less than 0.16?

# CONDITIONAL PROBABILITY AND BAYES' RULE =======================================

# 2. Joint, marginal, and conditional probability -------------------------------

trees <- matrix(c(39, 76, 101, 25), nrow = 2,
                dimnames = list(status = c("infested", "healthy"),
                                species = c("juniper", "pinyon")))
trees
p_trees <- trees / sum(trees)           # joint probabilities
round(p_trees, 3)

p_status  <- rowSums(trees) / sum(trees)
p_species <- colSums(trees) / sum(trees)
p_status
p_species

p_status["healthy"] * p_species["pinyon"]   # product of marginals
p_trees["healthy", "pinyon"]                # actual joint probability

#  P(A | B) = P(A and B) / P(B)
p_trees["healthy", "pinyon"] / p_species["pinyon"]    # = 25 / 126
p_trees["healthy", "juniper"] / p_species["juniper"]  # = 76 / 115

#  law of total probability: P(A) = P(A | B) P(B) + P(A | not B) P(not B)
(25 / 126) * p_species["pinyon"] + (76 / 115) * p_species["juniper"]
#  [PAUSE] widget: make the top of B level. what happens to P(A | B), compared with P(A)?

# 3. Bayes' rule ----------------------------------------------------------------
#  P(A | B) = P(B | A) P(A) / P(B)

p_healthy_given_pinyon <- 25 / 126
p_healthy_given_pinyon * p_species["pinyon"] / p_status["healthy"]
25 / 101                                    # straight from the table

#  medical test (Bolker 4.3.1): P(+) = P(+ | I) P(I) + P(+ | U) P(U)
p_inf <- 1e-6
p_pos_inf <- 1
p_pos_unf <- 0.01

p_pos <- p_pos_inf * p_inf + p_pos_unf * (1 - p_inf)
p_pos
p_inf_pos <- p_pos_inf * p_inf / p_pos
p_inf_pos
#  [PAUSE] widget: sensitivity or false positive rate, which matters more for P(present | +)? what if no false positives?

# 4. More than two hypotheses ---------------------------------------------------
#  P(H_i | D) = P(D | H_i) P(H_i) / sum_j P(D | H_j) P(H_j); liana classes (Bolker 4.3.2)

prior <- c(0.50, 0.25, 0.20, 0.05)
lik <- c(0.05, 0.10, 0.30, 0.60)       # P(3 infested neighbors | class)
post <- lik * prior / sum(lik * prior)
round(rbind(prior, likelihood = lik, posterior = post), 3)
sum(lik)                               # likelihoods need not add to 1
#  [PAUSE] with equal priors of 0.25, what is the posterior for heavy infestation?

barplot(rbind(prior, post), beside = TRUE, col = c("grey75", "steelblue"),
        border = NA, names.arg = c("liana-free", "light", "moderate", "heavy"),
        xlab = "Liana class", ylab = "Probability", ylim = c(0, 0.55))
legend("topright", c("prior", "posterior, 3 infested neighbors"),
       fill = c("grey75", "steelblue"), border = NA, bty = "n")

#  [PAUSE] sun detector lies with probability 1/36: what prior makes P(nova | "yes") > 0.5?
#  [PAUSE] 80% of failed nests had a predator on camera: P(failed | predator) = 0.8?
#  [PAUSE] two eDNA positives on the same day: how independent are they?

# Bonus: the Monty Hall problem -------------------------------------------------

prior_door <- rep(1/3, 3)
lik_door <- c(1/2, 0, 1)       # P(host opens door 2 | car behind door i)
lik_door * prior_door / sum(lik_door * prior_door)

# one game; returns TRUE if the final choice has the car
play_monty <- function(switch_door) {
  doors <- 1:3
  car <- sample(doors, 1)
  pick <- sample(doors, 1)
  goats <- setdiff(doors, c(pick, car))          # doors the host may open
  opened <- goats[sample(length(goats), 1)]
  if (switch_door) pick <- setdiff(doors, c(pick, opened))
  pick == car
}

set.seed(31)
mean(replicate(10000, play_monty(switch_door = FALSE)))
mean(replicate(10000, play_monty(switch_door = TRUE)))

# Bonus: many zeros from two processes ------------------------------------------
#  P(0) = (1 - v) + v (1 - p)^N; P(x > 0) = v choose(N, x) p^x (1 - p)^(N - x)

v <- 0.7
p_take <- 0.4
n_seeds <- 5
x_taken <- 0:n_seeds
p_zi <- v * dbinom(x_taken, n_seeds, p_take)
p_zi[1] <- p_zi[1] + (1 - v)
round(p_zi, 3)
sum(p_zi)

barplot(rbind(dbinom(x_taken, n_seeds, p_take), p_zi), beside = TRUE,
        names.arg = x_taken, col = c("grey75", "steelblue"), border = NA,
        xlab = "Seeds taken (of 5)", ylab = "Probability", ylim = c(0, 0.4))
