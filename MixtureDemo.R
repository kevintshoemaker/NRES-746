#install these if needed----
#install.packages("mlmRev")
#install.packages("loo")
#install.packages("bayesplot")
#install.packages("glmmTMB")
#install.packages("rstan")


#load in data/packages----
library(rstan) #tell people to add this and glmmTMB
library(cmdstanr)
library(loo)
library(bayesplot)
library(mlmRev)
library(glmmTMB)

#look at owl data
data(Owls)
str(Owls)
summary(Owls)

rstan_options(auto_write = TRUE) #tell rstan to automatically save compiled models for faster reuse

#Response variable: Sibling Negotiation (call count)
#Potential Predictors: FoodTreatment, Nest, SexParent, ArrivalTime, BroodSize, etc.
#define variables for the model
stan_data <- list(
  N = nrow(Owls),
  y = Owls$SiblingNegotiation,
  food = Owls$FoodTreatment,
  J = length(unique(Owls$Nest)),
  nest_id=Owls$Nest
)

#hurdle model-------

#compile stan model
hurdle_model = cmdstan_model("MixtureHurdleDemo.stan")

#fit the stan model
fit_hurdle <- hurdle_model$sample(
  data=stan_data,
  chains = 4,           # number of MCMC chains that are run at the same time
  parallel_chains = 4,  # number of computer cores used
  
  iter_warmup = 1000,   #number of warmup iterations per chain 
  iter_sampling = 1000, #number of iterations per chain (after the warm up)
)

#summarize what we just did
fit_hurdle$summary()

#visualize with trace plots--------
posterior_draws <- fit_hurdle$draws(c("alpha0", "alpha_food", "beta0", "beta_food", "sigma_a", "sigma_b"))
mcmc_trace(posterior_draws)

# model selection criteria ------

fit_hurdle$loo()


# zero inflated model -------
zim_model <- cmdstan_model("ZeroInflatedDemo.stan")

fit_zim <- zim_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

#summarize what we just did
fit_zim$summary()

# Compare fits
loo_hurdle <- loo(fit_hurdle$draws("log_lik", format="matrix"))
loo_zim <- loo(fit_zim$draws("log_lik", format="matrix"))
loo_compare(loo_hurdle, loo_zim)


###Exercise 1-------
#a. Create a zero inflated model in a new stan with the same variable that was used for the hurdle
#b. Bring the model to R and see how the model did 
#c. Compare LOOIC between the two models 

###Exercise 2-------
#a. Using a hurdle or zero-inflated model (whichever one did better), change or 
#add predictor variables in stan
#b. Compare your new model with the previous ones

###Questions---------
#1. Out of the three models we made today, which one was the best?
#2. Which variables did you add/change to make your model and how did that affect it's performance?



