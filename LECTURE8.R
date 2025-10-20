
#  NRES 746, Lecture 8                     
#   University of Nevada, Reno                       
#   Model selection and multi-model inference    -------------------      



# Load the balsam fir dataset
library(ggplot2)
library(emdbook)
data(FirDBHFec)
fir <- na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
fir$TOTCONES <- round(fir$TOTCONES)
head(fir)


plot(fir$TOTCONES ~ fir$DBH)   # fecundity as a function of tree size (diameter at breast height)


# tree fecundity by size, categorized into two site-level categories: "wave" and "non-wave" 

ggplot(fir,aes(DBH,TOTCONES)) + geom_point(aes(colour = WAVE_NON)) + theme_classic()


# build likelihood function for the full model: CONES ~ negBINOM( a(wave)*DBH^b(wave), dispersion(wave))
   
NLL_full <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos
  a <- c(params[1],params[2])[wave.code]     # a parameters (two for wave and one for non-wave)
  b <- c(params[3],params[4])[wave.code]      # b parameter (two for wave and one for non-wave)
  k <- c(params[5],params[6])[wave.code]       # over-dispersion parameters (two for wave and one for non-wave)
  expcones <- a*fir$DBH^b   # expected number of cones (deterministic component)
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))     # add stochastic component: full data likelihood
}

params <- c(a.n=1,a.w=1,b.n=1,b.w=1,k.n=1,k.w=1)

NLL_full(params)


## Find the MLE -----------------------

pars = c(a.n=1,a.w=1,b.n=1,b.w=1,k.n=1,k.w=1)
opt_full <- optim(pars,NLL_full,method="BFGS")

MLE_full <- opt_full$par
MinNLL_full <- opt_full$value
MLE_full


# build likelihood function for a reduced model: CONES ~ negBINOM( a(wave)*DBH^b, dispersion(wave))

NLL_constb <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos
  a <- c(params[1],params[2])[wave.code]      # a parameters
  b <- params[3]                              # b parameter (not a function of wave/nonwave)
  k <- c(params[4],params[5])[wave.code]      # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a.n=1,a.w=1,b=1,k.n=1,k.w=1)

NLL_constb(params)


### Find the MLE

opt_constb <- optim(fn=NLL_constb,par=c(a.n=1,a.w=1,b=1,k.n=1,k.w=1),method="L-BFGS-B")

MLE_constb = opt_constb$par

MinNLL_constb = opt_constb$value
MLE_constb


# compute -2*loglik for each model at the MLE

ms_full <- 2*MinNLL_full     # this is 2 * min.nll = -2*logLik_at_MLE
ms_constb <- 2*MinNLL_constb

ms_full
ms_constb


# Likelihood-Ratio test -----------------------

Deviance <- ms_constb - ms_full 
Deviance

Chisq.crit <- qchisq(0.95,1)
Chisq.crit

Deviance>=Chisq.crit   # perform the LRT

1-pchisq(Deviance,1)   # p-value


# Visualize the likelihood ratio test- compare the observed deviance with the distribution of deviances expected under the null hypothesis

curve(dchisq(x,df=1),0,5)
abline(v=Deviance,col="red",lwd=4)


# Try a different reduced model: CONES ~ negBINOM( a*DBH^b, dispersion)

NLL_nowave <- function(params){
  a <- params[1]      # a parameters
  b <- params[2]      # b parameter (not a function of wave/nonwave)
  k <- params[3]      # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a=1,b=1,k=1)

NLL_nowave(params)


# Find the MLE

opt_nowave <- optim(fn=NLL_nowave,par=params,method="L-BFGS-B")

MLE_nowave = opt_nowave$par

MinNLL_nowave = opt_nowave$value

MLE_nowave


# Perform LRT -- this time with three fewer free parameters in the reduced model

ms_full <- 2*MinNLL_full

ms_nowave <- 2*MinNLL_nowave

Deviance <- ms_nowave - ms_full 
Deviance

Chisq.crit <- qchisq(0.95,df=3)   # now three additional params in the more complex model!
Chisq.crit

Deviance>=Chisq.crit

1-pchisq(Deviance,df=3)   # p-value


# Visualize the likelihood ratio test (test statistic and sampling distribution under the null)
curve(dchisq(x,df=3),0,15)
abline(v=Deviance,col="red",lwd=4)


# Information-theoretic metrics for model-selection ------------------------------

# Akaike's Information Criterion (AIC)

## First, let's build another likelihood function: whereby only the "b" parameter differs by "wave" sites

NLL_constak <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos
  a <- params[1]                             # a parameters
  b <- c(params[2],params[3])[wave.code]                              # b parameter (not a function of wave/nonwave)
  k <- params[4]                               # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a=1,b.n=1,b.w=1,k=1)  

NLL_constak(params)


### Fit the new model

opt_constak <- optim(fn=NLL_constak,par=params)

MLE_constak= opt_constak$par

MinNLL_constak = opt_constak$value

ms_constak <- 2*MinNLL_constak
MLE_constak

### Now, let's build and fit one more final model- this time with no wave effect and a Poisson error distribution

PoisLik_nowave <- function(params){
  a <- params[1]      # a parameters
  b <- params[2]      # b parameter (not a function of wave/nonwave)
  expcones <- a*fir$DBH^b
  -sum(dpois(fir$TOTCONES,lambda=expcones,log=TRUE))
}

params <- c(a=1,b=1)

PoisLik_nowave(params)

opt_pois <- optim(fn=PoisLik_nowave,par=params)

MLE_pois= opt_pois$par

MinNLL_pois= opt_pois$value

ms_pois <- 2*MinNLL_pois

MLE_pois


# Compare all five models using AIC!

AIC_constak <- ms_constak + 2*4
AIC_full <- ms_full + 2*6
AIC_constb <- ms_constb + 2*5
AIC_nowave <- ms_nowave + 2*3
AIC_pois <- ms_pois + 2*2

AICtable <- data.frame(
  Model = c("Full","Constant b","Constant a and k","All constant","Poisson"),
  AIC = c(AIC_full,AIC_constb,AIC_constak,AIC_nowave,AIC_pois),
  LogLik = c(ms_full/-2,ms_constb/-2,ms_constak/-2,ms_nowave/-2,ms_pois/-2),
  params = c(6,5,4,3,2),
  stringsAsFactors = F
)

AICtable$DeltaAIC <- AICtable$AIC-AICtable$AIC[which.min(AICtable$AIC)]

AICtable$Weights <- round(exp(-0.5*AICtable$DeltaAIC) / sum(exp(-0.5*AICtable$DeltaAIC)),3)

AICtable$AICc <- AICtable$AIC + ((2*AICtable$params)*(AICtable$params+1))/(nrow(fir)-AICtable$params-1)

AICtable[order(AICtable$AIC),c(1,7,2,5,6,4,3)]

# Bayes factor example  ---------------------

# take a basic binomial distribution with parameter p fixed at 0.5:

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability")

## Q: What is the *marginal likelihood* under this simple model for an observation of 2 mortalities out of 10? 

## A:

dbinom(2,10,0.5)



## Now we can consider a model whereby "p" is a free parameter
curve(dbeta(x,1,1))  # uniform prior on "p"


# Compute the marginal likelihood of observing 2 mortalities

# ?integrate
binom2 <- function(x) dbinom(x=2,size=10,prob=x)
marginal_likelihood <- integrate(f=binom2,0,1)$value    # use "integrate" function in R
marginal_likelihood  # equal to 0.0909 = 1/11


# Compute the marginal likelihood of observing 3 mortalities

binom3 <- function(x) dbinom(x=3,size=10,prob=x)
marginal_likelihood <- integrate(f=binom3,0,1)$value    # use "integrate" function
marginal_likelihood   # equal to 0.0909 = 1/11


# simulate data from the model across all possible values of the parameter "p"

lots=1000000
a_priori_data <- rbinom(lots,10,prob=rbeta(lots,1,1))   # no particular observation is favored
for_hist <- table(a_priori_data)/lots
barplot(for_hist,xlab="Potential Observation",ylab="Marginal likelihood")


# Visualize the marginal likelihood of all possible observations

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))


# Overlay the marginal likelihood of the simpler model, with p fixed at 0.5

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability",add=T,col="red",density=20)


# Finally, compute the bayes factor given that we observed 2 mortalities. Which model is better?

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability",add=T,col="red",density=20)

abline(v=3,col="green",lwd=4 )


BayesFactor = (1/11)/dbinom(2,10,0.5)   
BayesFactor


# Compute the bayes factor given that we observed 3 mortalities. Which model is better now?

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability",add=T,col="red",density=20)

abline(v=4.3,col="green",lwd=4 )


BayesFactor = dbinom(3,10,0.5)/(1/11)
BayesFactor


# Visualize the likelihood ratio

# probs2 <- rep(1/11,times=11)          
# names(probs2) = 0:10
# barplot(probs2,ylab="probability",ylim=c(0,1))

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability",col="red",density=20,ylim=c(0,1))

probs3 <- dbinom(0:10,10,0.3)          
names(probs3) = 0:10
barplot(probs3,ylab="probability",add=T,col="green",density=10,angle = -25)

abline(v=4.3,col="green",lwd=4 )


# LRT: simple model (p fixed at 0.5) vs complex model (p is free parameter)

Likelihood_simple <- dbinom(3,10,0.5)
Likelihood_complex <- dbinom(3,10,0.3)
Likelihood_simple
Likelihood_complex
-2*log(Likelihood_simple)--2*log(Likelihood_complex)

qchisq(0.95,1)

pchisq(1.64,1)    # very high p value, simpler model is preferred


# AIC: simple model (p fixed at 0.5) vs complex model (p is free parameter)

AIC_simple <- -2*log(Likelihood_simple) + 2*0
AIC_complex <-  -2*log(Likelihood_complex) + 2*1

AIC_simple
AIC_complex    


### Alternatively, use AICc

AICc_simple <- -2*log(Likelihood_simple) + 0 + 0
AICc_complex <-  -2*log(Likelihood_complex) + 1 + ((2*2)/(3-1-1))

AICc_simple
AICc_complex    


# Alternatively, try BIC

BIC_simple <- -2*log(Likelihood_simple) + log(10)*0
BIC_complex <-  -2*log(Likelihood_complex) + log(10)*1

BIC_simple
BIC_complex    



data {
  int<lower = 1> N;
  array [N] int<lower=0> obs_cones;
  vector<lower=0> [N] DBH;
  array [N] int<lower=1,upper=2> wave_ndx; 
}

parameters {
  vector[2] loga, b, logbeta;
}

transformed parameters {
  vector[2] a = exp(loga);
  vector[2] beta = exp(logbeta);
}

model  {
  vector[N] mean_cones = exp(loga[wave_ndx] + b[wave_ndx] .* DBH);   // power function: a*DBH^b
  vector[N] alpha = mean_cones .* beta[wave_ndx];
  obs_cones ~ neg_binomial(alpha,beta[wave_ndx]);
}

generated quantities {   // need log_lik of each data point for model selection
  vector[N] log_lik; // N is the number of data points
  {
    real m2, a2, b2;
    for (n in 1:N) {
       m2 = exp(loga[wave_ndx[n]] + b[wave_ndx[n]] * DBH[n]);   // power function: a*DBH^b
       a2 = m2 * beta[wave_ndx[n]];
       log_lik[n] = neg_binomial_lpmf(obs_cones[n] | a2, beta[wave_ndx[n]]);
    }
  }
}




# Package the data for stan

stan_data <- list(
  N = nrow(fir),
  obs_cones = fir$TOTCONES,
  wave_ndx = as.numeric(fir$WAVE_NON),
  DBH = fir$DBH
)
#data.package


# Run the model in stan

library(cmdstanr)    # load packages
library(bayesplot)
library(posterior)
options(mc.cores = 4) 

# compile the model using:
firmodel_full <- cmdstan_model("firmodel_full.stan") # Compile stan model

fit_full <- suppressMessages( firmodel_full$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)  )

# fit_full$summary()
samples <- fit_full$draws(format="draws_df")
bayesplot::mcmc_trace(samples,"a[1]")
bayesplot::mcmc_trace(samples,"b[2]")
bayesplot::mcmc_trace(samples,"beta[1]")


meanrep_wave = exp(samples$`loga[2]` + samples$`b[2]`*mean(stan_data$DBH))
meanrep_nonwave = exp(samples$`loga[1]` + samples$`b[1]`*mean(stan_data$DBH))

hist(meanrep_nonwave,ylab="Prob Density",xlab="number of cones",freq = F,xlim=c(25,60),ylim=c(0,0.15),main="")
hist(meanrep_wave,density=20,col="darkgreen",add=T,freq=F)
legend("topleft",col=c("darkgreen","white"),density=c(20,0),legend=c("wave","nonwave"),bty="n")


# Extract the WAIC for the full model!

library(loo)

loglik_names <- sapply (1:stan_data$N,
                          function(t) sprintf("log_lik[%s]",t) )
lls = as.matrix(samples)[,loglik_names]

waic_full = loo::waic(lls)
waic_full$estimates["waic",]


loo_full <- fit_full$loo()
loo_full


data {
  int<lower = 1> N;
  array [N] int<lower=0> obs_cones;
  vector<lower=0> [N] DBH;
}

parameters {
  real loga, b, logbeta;
}

transformed parameters {
  real a = exp(loga);
  real beta = exp(logbeta);
}

model  {
  vector[N] mean_cones = exp(loga + b .* DBH);   // power function: a*DBH^b
  vector[N] alpha = mean_cones .* beta;
  obs_cones ~ neg_binomial(alpha,beta);
}

generated quantities {   // need log_lik of each data point for model selection
  vector[N] log_lik; // N is the number of data points
  {
    real m2, a2, b2;
    for (n in 1:N) {
       m2 = exp(loga + b * DBH[n]);   // power function: a*DBH^b
       a2 = m2 * beta;
       log_lik[n] = neg_binomial_lpmf(obs_cones[n] | a2, beta);
    }
  }
}

# compile the model using:
firmodel_reduced <- cmdstan_model("firmodel_reduced.stan") # Compile stan model

fit_reduced <- suppressMessages( firmodel_reduced$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 0
)  )

# fit_full$summary()
samples_red <- fit_reduced$draws(format="draws_df")
bayesplot::mcmc_trace(samples_red,"a")
# bayesplot::mcmc_trace(samples_red,"b")
# bayesplot::mcmc_trace(samples_red,"beta")

# Compute LOOIC

loo_reduced <- fit_reduced$loo()

loo_reduced


# Goodness of fit

N <- stan_data$N

plot(fir$TOTCONES~fir$DBH,ylim=c(0,900),cex=2)

# which.max(stan_data$DBH)
# which.min(stan_data$DBH)

a_param <- samples_red$a
b_param <- samples_red$b
beta_param <- samples_red$beta
mu_rep <- sapply(1:N,function(t) a_param * stan_data$DBH[t]^b_param)

sim_dat <- function(s){
  thismean = mu_rep[s,]
  thisbeta <- samples_red$beta[s]
  thisalpha <- thismean * thisbeta 
  thismu = thisalpha / thisbeta
  rnbinom(N,size=thisalpha,mu=thismu)
}

nMCMC = length(samples_red$b)

for(d in 1:N){
  dat= sim_dat(sample(1:nMCMC,1))
  points(stan_data$DBH,dat,pch=20,col="gray",cex=0.4)
}
points(fir$DBH,fir$TOTCONES,cex=2)


# Posterior Predictive Checks!

nreps = 500 
SSE_obs= numeric(nreps)
SSE_sim = numeric(nreps)
r=1
for(r in 1:nreps){
  this= sample(1:nMCMC,1,replace = T)
  SSE_obs[r] = sum((stan_data$obs_cones - mu_rep[this,])^2)
  simdat = sim_dat(this)
  SSE_sim[r] = sum((simdat - mu_rep[this,])^2)
}

plot(SSE_sim~SSE_obs,xlab="SSE, real data",ylab="SSE, simulated data",main="Posterior Predictive Check")
abline(0,1,col="red")
p.value=mean(SSE_sim>SSE_sim)
p.value 

