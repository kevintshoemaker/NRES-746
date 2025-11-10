
## [[NOTE: claude.ai provided the initial draft of this code]]

rm(list=ls())

# Load required libraries -----

library(rstan)
library(cmdstanr)
library(titanic)  # For titanic dataset
library(dplyr)
library(faux)
library(ISLR)

options(mc.cores=4)

# Load and prepare data ----------

# data("titanic_train")
# df <- titanic_train %>%
# select(Survived, Pclass, Sex, Age, SibSp, Parch, Fare) %>%
# na.omit()

# data(Hitters)            # alternative: hitters dataset?
# df <- na.omit(Hitters)
# dim(df)

# Prepare predictor matrix ------------

# X <- model.matrix(~ scale(Fare), data = df)[, -1, drop=F]  #  Sex + scale(Age) + scale(SibSp) + scale(Parch) +
# X <- cbind(X, rnorm_multi(n = nrow(X), mu = 0, sd = 1, r = 0.99, vars = 10) )         #matrix(runif(10*nrow(X),-2,2),ncol=10))

# X = scale(model.matrix(Salary ~ ., df)[, -1, drop=F])
# rownames(X) = NULL
# colnames(X) = NULL

#  made up data...

set.seed(1234)
n = 500
p = 750
X = replicate(p, rnorm(n = n))
dim(X)
beta = c(1, 1, 1, rep(0, p-3))
z = X %*% beta
prob = exp(z) / (1 + exp(z))
y = rbinom(length(z), size = 1, prob = prob)


# prepare other variables -----

# y <- df$Survived
# y <- df$Salary

# class = df$Pclass

cc = rep(1, times=nrow(X))

N <- nrow(X)
K <- ncol(X)
nclass = max(cc)

set.seed(2025)

test_frac  = .5

ndx_test = sample(1:nrow(X),nrow(X)*test_frac)
ndx_train = setdiff(1:nrow(X),ndx_test)

# Prepare data for Stan
stan_data <- list(
  N = length(ndx_train),
  N_test = length(ndx_test),
  K = K,
  C = nclass,
  cc = cc[ndx_train], #df$Pclass[ndx_train],
  cc_test = cc[ndx_test], #df$Pclass[ndx_test],
  X = X[ndx_train,],
  X_test = X[ndx_test,],
  y = y[ndx_train],
  y_test = y[ndx_test],
  lambda = 2  # Regularization strength (higher is more regularized/less complex)
)

mod = cmdstan_model("stan_demo_titanic.stan")

# Fit the model
# fit <- stan(
#   model_code = "stan_demo_titanic.stan",
#   data = stan_data,
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   seed = 123
# )

fit = mod$sample(
  data=stan_data,
  chains=4,
  iter_warmup=500,
  iter_sampling=500
)

# Print results
# print(fit, pars = c("alpha", "beta"))
a= fit$summary()

# Extract posterior samples
draws = fit$draws(format="draws_df")

names(draws)

library(bayesplot)

mcmc_trace(draws,"alpha[1]")
# mcmc_trace(draws,"alpha[2]")
# mcmc_trace(draws,"alpha[3]")
mcmc_trace(draws,"beta[400]")

# Plot coefficient distributions

# coefnames= colnames(X)
# coef_ndx = sprintf("beta[%s]",1:length(coefnames))

# library(ggplot2)
# sapply(1:length(coef_ndx), function(t) mcmc_dens(draws,coef_ndx[t]) + labs(title=coefnames[t])  )


# Evaluate out-of-sample performance -------------

mcmc_dens(draws,"accuracy")
mean(draws$accuracy)   # 75% accuracy with lambda = 1



# optimize lambda 
l=1
get_accuracy = function(l=1.0){
  stan_data$lambda = l
  
  thisfit = mod$sample(
    data=stan_data,
    chains=4,
    iter_warmup=500,
    iter_sampling=500
  )
  
  # Extract posterior samples
  draws = thisfit$draws(format="draws_df")
  c(quantile(draws$accuracy,c(0.5)),quantile(draws$tp,c(0.5) ) ) 
  
}

lams= c(0.1,1,2,5,10,50,100)
tune_lambda = sapply(lams, get_accuracy)

cbind(lambda=lams, t(tune_lambda))

plot(log(lams),tune_lambda[1,],type="l")

lams[which.max(tune_lambda[2,])]
