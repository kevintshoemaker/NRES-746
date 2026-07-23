# ---- required libraries ----

library(rstan) #interface to stan for bayesian modeling
library(stringr) #for easy string manipulation
library(bayesplot) #for mcmc diagnostic and visualization tools
library(magrittr) #provides the pipe (%>%) operator 
library(cmdstanr)
library(ggplot2)

rstan_options(auto_write = TRUE) #tell rstan to automatically save compiled models for faster reuse
#use all available cpu cores for faster sampling 
options(mc.cores = 4)

# ---- simulation settings ----
N=50; K=5
sim_values <- list(
  N = N,  # observations
  K = K,    # groups
  X = rnorm(N, 0, 5), #predictor variable
  gg = sample(1:K, N, replace = TRUE), #group IDs (1 to 5), randomly assigned
  tau = 0.1, #standard deviation of group-level effects
  sigma = 1, #residual standard deviation
  alpha0 = 2, #true intercept
  b1 = 0.5 #true slope 
)

# ---- generate data (Fixed_param) ----

simmod = cmdstan_model("data_simyes.stan")
sim_data1 = simmod$sample(data=sim_values,chains=1,iter_sampling=1,fixed_param=TRUE,seed=42)

#extract the simulated data and paramters from the stan output
sim_draws <- sim_data1$draws(format="draws_df")

#pull the simulated response variable (y)
sim_y     <- sapply(1:sim_values$N,function(t) sim_draws[[sprintf("y[%s]",t)]]  )

#pull the simulated group-level effects (gamma)
sim_alpha <- sapply(1:sim_values$K,function(t) sim_draws[[sprintf("alpha[%s]",t)]]  )

plot(sim_y~sim_values$X)

# ---- prepare data for model fit ----

     #create a list to pass to stan for model fitting 
data_for_fit <- list( 
  N = sim_values$N, #number of observations
  K = sim_values$K, #number of groups
  y = sim_y, #simulated response variable
  gg = sim_values$gg, #group IDs
  X = sim_values$X #predictor variable 
)


# ---- CENTERED MODEL ----

mod_ctr = cmdstan_model("centered_modelyes.stan")

fit_ctr = mod_ctr$sample(data=data_for_fit,chains=4,iter_warmup = 500, iter_sampling=500)

draws_ctr = fit_ctr$draws(format="draws_df")

fit_ctr$summary()

# ---- Trace plots for alpha ----

library(bayesplot)

K <- data_for_fit$K  # number of groups
#create strings like "gamma[1]", "gamma[2]",..., for plotting
alpha_string <- str_c("alpha[", 1:K, "]")

#plot mcmc trace plots for gamma, tau, and sigma (diagnostic check for convergence)
bayesplot::mcmc_trace(draws_ctr, pars = c("tau", "sigma"))
bayesplot::mcmc_trace(draws_ctr, pars = c(alpha_string,"tau", "sigma"))+
  ggtitle("Centered graph")

# ---- NON-CENTERED MODEL ----

mod_nctr = cmdstan_model("noncenteredyes.stan")

fit_nctr = mod_nctr$sample(data=data_for_fit,chains=4,iter_warmup = 500, iter_sampling=500)

draws_nctr = fit_nctr$draws(format="draws_df")

fit_nctr$summary()

# ---- Trace plots for alpha ----

library(bayesplot)
bayesplot::mcmc_trace(draws_nctr,"alpha0")

K <- data_for_fit$K  # number of groups

#plot mcmc trace plots for gamma, tau, and sigma (diagnostic check for convergence)
bayesplot::mcmc_trace(draws_nctr, pars = c(alpha_string,"tau", "sigma"))+
  ggtitle("Noncentered graph")

# ---- PARAMETER RECOVERY PLOTS ----
library(tidyverse) #data manipulation and plotting 
library(tidybayes) #for tidy extraction of bayesian model results 

# Create a tibble of the true Gamma values
alpha_values <- tibble(
  group = 1:data_for_fit$K, #group numbers (1 through K)
  .variable = str_c("alpha[", group, "]"), #label each group
  values = sim_alpha #the true Gamma values from the simulated data
)

# ---- Centered Model: Gamma Recovery ----
fit_ctr %>%
  gather_draws(alpha[k]) %>%                       # Extract posterior draws for each Gamma
  mutate(.variable = str_c("alpha[", k, "]")) %>%  # Label each group
  ggplot(aes(x = .value, y = .variable)) +        #plot the posterior draws for each gamma
  geom_halfeyeh(.width = 0.95,                    #draw 95% credible intervals
                fill = "skyblue", alpha = 0.6) +  # fill color for centered model
  geom_vline(aes(xintercept = values),            
             alpha_values, color = "red", linewidth = 1) +  # add red lines showing "true value"
  facet_wrap(~ .variable, ncol = 1, scales = "free_y") +  #make a separate panel for each group
  labs(
    title = "Alpha Recovery (Centered Model)",
    x = "Posterior Estimates",
    y = NULL
  ) +
  theme_bw()

# ---- Non-Centered Model: Gamma Recovery ----
fit_nctr %>%
  gather_draws(alpha[k]) %>%                       # Extract posterior draws for each Gamma
  mutate(.variable = str_c("alpha[", k, "]")) %>%  # Label each group
  ggplot(aes(x = .value, y = .variable)) + #plot posteriors by group
  geom_halfeyeh(.width = 0.95, fill = "lightgreen", alpha = 0.6) + #green fill for non-centered
  geom_vline(aes(xintercept = values), alpha_values, color = "red", linewidth = 1) + #true values
  facet_wrap(~ .variable, ncol = 1, scales = "free_y") + #separate panels for each group 
  labs(
    title = "alpha Recovery (Non-Centered Model)", #plot title 
    x = "Posterior Estimates",
    y = NULL
  ) +
  theme_bw()
