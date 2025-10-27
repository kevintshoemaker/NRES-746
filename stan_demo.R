
# In class stan demo examples for NRES 746

rm(list=ls())

rstan=FALSE

# load packages --------

library(cmdstanr)
library(rstudioapi)
library(rstan)
library(bayesplot)
library(posterior)
library(ggplot2)
library(mvtnorm)

options(mc.cores=4)

# set global simulation parameters -----------

N=100                            # sample size
sig = 0.5                        # residual standard dev 
alpha0 = -1                      # global mean intercept
b1 = 1.1                         # global mean slope term for effect of covariate x1
G = 20                           # 
tau = c(0.3,0.4)                 # hyperparam for random intercept
rho = 0.5                        # hyperparam for correlation among slope and intercept

# generate data -------------

library(mvtnorm)

x1 = runif(N)                         # generate x1 covariate
gg = sample(1:G,N,replace = T)        # group index for each observation
vcv = diag(tau)
vcv[1,2] = rho; vcv[2,1] = rho
alpha = rmvnorm(G,c(alpha0,b1), )           # generate random slopes and intercepts for each group
y = rnorm(N, alpha[gg,1]+alpha[gg,2]*x1, sig)    # generate scalar response y 
plot(y~x1)                            # plot to make sure it looks right!

df = data.frame(   # package data into data frame for later
  y = y,
  x1 = x1,
  G = factor(gg,levels=1:G)
)

ggplot(df,aes(x1,y,colour = G)) + geom_point() + theme_classic()


# do prior predictive checks --------

library(ggdist)
K=2
eta =4

     # explore lkj correlation prior
{
  temp = rlkjcorr_marginal(n=sum(2:K),K=K,eta=eta)
  mat = diag(rep(1,K))
  ctr = 0
  for(r in 1:(K-1)){ mat[r,(r+1):K] = temp[(ctr+1):(ctr+(K-r))] ; ctr= (ctr+(K-r)) }
  mat
}

# package data for stan ---------

stan_data <- list(
  N = N,
  G = G,
  gg = gg,
  x=x1,
  y=y
)

# compile and fit stan model -------

if(!rstan){
  stanmod = cmdstan_model("stan_demo2.stan") # Compile stan model (cmdstanr)
  fit <- stanmod$sample(
    data = stan_data,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 500
  ) 
}else{
  stanmod <- stan_model("stan_demo2.stan")
  fit <- sampling(stanmod, data = stan_data, chains=4, iter=1000)
}

# summary of posterior samples --------

if(!rstan){
  fit$summary()
}else{
  a=summary(fit)
  a$summary
}

# package samples for further analysis ---------

if(!rstan){
  samples <- fit$draws(format="draws_df")
}else{
  samples <- as_draws_df(fit)   # convert to 'draws' object for visualization with bayesplot
}

# visualize posterior --------

bayesplot::mcmc_trace(samples,"alpha0")
bayesplot::mcmc_trace(samples,"alpha[2]")

bayesplot::mcmc_dens(samples,"alpha0") + geom_vline(xintercept=alpha0,lwd=2)

bayesplot::mcmc_dens(samples,"beta0") + geom_vline(xintercept=b1,lwd=2)

bayesplot::mcmc_dens(samples,"sigma") + geom_vline(xintercept=sig,lwd=2)

bayesplot::mcmc_dens(samples,"tau[1]") + geom_vline(xintercept=tau[1],lwd=2)
bayesplot::mcmc_dens(samples,"tau[2]") + geom_vline(xintercept=tau[2],lwd=2)

bayesplot::mcmc_pairs(samples,pars=c("alpha0","alpha[2]") )


# visualize random effect and shrinkage -----------

# true group-level intercept
alpha

# estimated group-level intercept

names(samples)
alpha_mc = sapply(1:G, function(t) samples[[sprintf("alpha[%s]",t)]] )

library(tidyr)
alpha_mc2 = pivot_longer(as.data.frame(alpha_mc),everything(), names_to = "G", values_to = "alpha")
alpha_mc2$G = gsub("V","",alpha_mc2$G)
alpha_mc2$G = factor(alpha_mc2$G,levels=c(1:G))
# View(alpha_mc2)

library(dplyr)
alpha_mc3 = alpha_mc2 |> 
  group_by(G) |> 
  summarise(alpha=mean(alpha))

sampsize = data.frame(G=1:G,N=as.numeric(table(df$G)),y=min(alpha_mc2$alpha) )

mod1=lm(y~0+x1+G,df)   # estimate of intercept with no shrinkage

df2 = data.frame(G=1:G,alpha=coef(mod1)[-1]) 

## visualize shrinkage... 
ggplot(alpha_mc2,aes(G,alpha)) + geom_violin(fill=gray(0.7),colour=NA) +
  geom_point(data=alpha_mc3,aes(G,alpha),pch="X",col="darkblue",size=3) +
  geom_point(data=df2,aes(G,alpha),size=3) +
  geom_point(data=data.frame(G=(1:G),alpha=alpha),aes(G,alpha),size=3,pch="-",alpha=.5,) +
  geom_hline(yintercept = mean(samples$alpha0),col="darkgreen",lwd=2) +
  geom_text(data=sampsize,aes(x=G,y=y,label=N)) +
  ylim(min(alpha_mc2$alpha)-0.1,max(alpha_mc2$alpha)+0.1) +
  theme_classic()



# posterior predictive check --------


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

plot(post_pred2$RMSE_sim~post_pred2$RMSE_obs, main="posterior predictive check",
         ylab="RMSE, simulated", xlab="RMSE, observed")
abline(0,1,col="red",lwd=2)

p_val = with(post_pred2, mean(RMSE_sim > RMSE_obs)  )
p_val



## try running model in lme4



mod

























