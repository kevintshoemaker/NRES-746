
library(cmdstanr)
library(bayesplot)
library(posterior)
library(ggplot2)

options(mc.cores=4)

N=200
sig = 0.5
b0 = -1
b1 = 1.1
x = runif(N)
y = rnorm(N,b0+b1*x,sig)
plot(y~x)


stan_data <- list(
  N = N,
  x=x,
  y=y
)

stanmod = cmdstan_model("stan_demo1.stan") # Compile stan model

fit <- stanmod$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
) 

fit$summary()

samples <- fit$draws(format="draws_df")
bayesplot::mcmc_trace(samples,"b0")

bayesplot::mcmc_dens(samples,"b0") + geom_vline(xintercept=b0,lwd=2)

bayesplot::mcmc_dens(samples,"b1") + geom_vline(xintercept=b1,lwd=2)

bayesplot::mcmc_dens(samples,"sigma") + geom_vline(xintercept=sig,lwd=2)




