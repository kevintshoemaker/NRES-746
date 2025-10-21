
# In class stan demo examples for NRES 746

# load packages --------

library(cmdstanr)
library(bayesplot)
library(posterior)
library(ggplot2)
library(mvtnorm)

options(mc.cores=4)

# set global simulation parameters -----------

N=200                            # sample size
sig = 0.5                        # residual standard dev 
alpha0 = -1                      # global mean intercept
b1 = 1.1                         # global mean slope term for effect of covariate x1
G = 15                           # 
tau = 0.25                       # hyperparam for random intercept

# generate data -------------

x1 = runif(N)                     # generate x1 covariate
gg = sample(1:G,N,replace = T)    # group index for each observation
alpha = rnorm(G,alpha0,tau)       # generate random intercepts for each group
y = rnorm(N, alpha[gg]+b1*x1, sig)    # generate scalar response y 
plot(y~x1)                        # plot to make sure it looks right!

ggplot(data.frame(y=y,x=x1,G=as.factor(gg)),aes(x,y,colour = G)) + geom_point() + theme_classic()


# package data for stan ---------

stan_data <- list(
  N = N,
  G = G,
  gg = gg,
  x1=x1,
  y=y
)

# compile stan model -------
stanmod = cmdstan_model("stan_demo1.stan") # Compile stan model

fit <- stanmod$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
) 

# summary of posterior samples --------

fit$summary()

# package samples for further analysis ---------

samples <- fit$draws(format="draws_df")

# visualize posterior --------

bayesplot::mcmc_trace(samples,"alpha0")
bayesplot::mcmc_trace(samples,"alpha[2]")

bayesplot::mcmc_dens(samples,"alpha0") + geom_vline(xintercept=alpha0,lwd=2)

bayesplot::mcmc_dens(samples,"b1") + geom_vline(xintercept=b1,lwd=2)

bayesplot::mcmc_dens(samples,"sigma") + geom_vline(xintercept=sig,lwd=2)

bayesplot::mcmc_dens(samples,"tau") + geom_vline(xintercept=tau,lwd=2)


bayesplot::mcmc_pairs(samples,"alpha0","tau")

bayesplot::mcmc_pairs(samples,pars=c("alpha0","alpha[2]") )


# posterior predictive check --------

df = data.frame(
  y = y,
  x1 = x1,
  G = gg
)

# set range of x1
x1_seq=seq(0,1,length=10)
gg_seq = rep(1,length=10)
nMCMC = length(samples$sigma)

# make predictions from posterior
n_samp = 500
do_pred = function(){
  grab = samples[sample(1:nMCMC,1),]     # grab sample from joint posterior
  thisalpha = sapply(gg_seq, function(t) grab[[sprintf("alpha[%s]",t )]] ) 
  thismu = with(grab, thisalpha + b1 *x1_seq)
  with(grab,rnorm(length(x1_seq), thismu, sigma) )
}
post_pred = as.data.frame(t(replicate(n_samp,do_pred())))

df_sub= subset(df,G==1)
plot(y~x1,type="n",data=df_sub,ylim=c(-3,2),xlim=c(-0.1,1.1))

expected <- with(samples, mean(alpha[1]) + mean(b1) * x1_seq) 
points(x1_seq,expected,type="l",col="darkred",lwd=2)

boxplot(x=as.list(post_pred),at=x1_seq,add=T,boxwex=0.025,xaxt="n",range=0,border="darkred")
points(df_sub$x1,df_sub$y, cex=1.5,pch=20)


## Bayesian p-value ----------

do_pred2 = function(){
  grab = samples[sample(1:nMCMC,1),]
  thisalpha = sapply(df$G, function(t) grab[[sprintf("alpha[%s]",t )]] ) 
  thismu = with(grab, thisalpha + b1 *df$x1)
  thissim = with(grab, rnorm(nrow(df),thismu,sigma ))
  RMSE_obs = sqrt(mean((df$y-thismu)^2))
  RMSE_sim = sqrt(mean((thissim-thismu)^2))
  c(RMSE_obs = RMSE_obs,RMSE_sim = RMSE_sim)
}
post_pred2 = as.data.frame(t(replicate(n_samp,do_pred2())))

plot(post_pred2$RMSE_sim~post_pred2$RMSE_obs, main="posterior predictive check")
abline(0,1,col="red",lwd=2)


























