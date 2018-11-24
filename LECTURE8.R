
############################################################
####                                                    ####  
####  NRES 746, Lecture 8                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Model selection and multi-model inference         ####
############################################################



#######
# Load the balsam fir dataset (finally, no more rabbits and virus titers!)

library(emdbook)
data(FirDBHFec)
fir <- na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
fir$TOTCONES <- round(fir$TOTCONES)
head(fir)


plot(fir$TOTCONES ~ fir$DBH)   # fecundity as a function of tree size (diameter at breast height)


#########
# tree fecundity by size, categorized into two site-level categories: "wave" and "non-wave" 

ndx <- fir$WAVE_NON=="w"   # logical vector indicating which observations were from "wave" sites
plot(fir$TOTCONES[ndx] ~ fir$DBH[ndx],xlab="DBH",ylab="Tot Cones")
points(fir$DBH[!ndx],fir$TOTCONES[!ndx],pch=4,col="red")
legend("topleft",pch=c(1,4),col=c("black","red"),legend=c("Wave","Non-wave"),bty="n")


########
# build likelihood function for the full model: CONES ~ negBINOM( a(wave)*DBH^b(wave), dispersion(wave))

   
NegBinomLik_full <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos    # note: we are hard-coding the data into our likelhood function here!
  a <- c(params[1],params[2])[wave.code]     # a parameters (one for wave and one for non-wave)
  b <- c(params[3],params[4])[wave.code]      # b parameter (one for wave and one for non-wave)
  k <- c(params[5],params[6])[wave.code]       # over-dispersion parameters (one for wave and one for non-wave)
  expcones <- a*fir$DBH^b   # expected number of cones (deterministic component)
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))     # add stochastic component: full data likelihood
}

params <- c(a.n=1,a.w=1,b.n=1,b.w=1,k.n=1,k.w=1)

NegBinomLik_full(params)


#### Find the MLE

MLE_full <- optim(fn=NegBinomLik_full,par=c(a.n=1,a.w=1,b.n=1,b.w=1,k.n=1,k.w=1),method="L-BFGS-B")

MLE_full$par

MLE_full$value


########
# build likelihood function for a reduced model: CONES ~ negBINOM( a(wave)*DBH^b, dispersion(wave))


NegBinomLik_constb <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos
  a <- c(params[1],params[2])[wave.code]      # a parameters
  b <- params[3]                              # b parameter (not a function of wave/nonwave)
  k <- c(params[4],params[5])[wave.code]      # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a.n=1,a.w=1,b=1,k.n=1,k.w=1)

NegBinomLik_constb(params)


### Find the MLE

MLE_constb <- optim(fn=NegBinomLik_constb,par=c(a.n=1,a.w=1,b=1,k.n=1,k.w=1),method="L-BFGS-B")

MLE_constb$par

MLE_constb$value


#######
# compute -2*loglik for each model at the MLE

ms_full <- 2*MLE_full$value     # this is 2 * min.nll = -2*logLik_at_MLE

ms_constb <- 2*MLE_constb$value

ms_full
ms_constb


#############
# Likelihood-Ratio test (frequentist)

Deviance <- ms_constb - ms_full 
Deviance

Chisq.crit <- qchisq(0.95,1)
Chisq.crit

Deviance>=Chisq.crit   # perform the LRT

1-pchisq(Deviance,1)   # p-value


####### Visualize the likelihood ratio test- compare the observed deviance with the distribution of deviances expected under the null hypothesis

curve(dchisq(x,df=1),0,5)
abline(v=Deviance,col="red",lwd=4)


#############
# Try a different reduced model: CONES ~ negBINOM( a*DBH^b, dispersion)

NegBinomLik_nowave <- function(params){
  a <- params[1]      # a parameters
  b <- params[2]      # b parameter (not a function of wave/nonwave)
  k <- params[3]      # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a=1,b=1,k=1)

NegBinomLik_nowave(params)


### Find the MLE

MLE_nowave <- optim(fn=NegBinomLik_nowave,par=params,method="L-BFGS-B")

MLE_nowave$par

MLE_nowave$value


#########
# Perform LRT -- this time with three fewer free parameters in the reduced model

ms_full <- 2*MLE_full$value

ms_nowave <- 2*MLE_nowave$value

Deviance <- ms_nowave - ms_full 
Deviance

Chisq.crit <- qchisq(0.95,df=3)   # now three additional params in the more complex model!
Chisq.crit

Deviance>=Chisq.crit

1-pchisq(Deviance,df=3)   # p-value


###### Visualize the likelihood rato test
curve(dchisq(x,df=3),0,15)
abline(v=Deviance,col="red",lwd=4)


##############
# Information-theoretic metrics for model-selection
##############


#########
# Akaike's Information Criterion (AIC)

#### First, let's build another likelihood function: whereby only the "b" parameter differs by "wave" sites

NegBinomLik_constak <- function(params){
  wave.code <- as.numeric(fir$WAVE_NON)      # convert to ones and twos
  a <- params[1]                             # a parameters
  b <- c(params[2],params[3])[wave.code]                              # b parameter (not a function of wave/nonwave)
  k <- params[4]                               # dispersion parameters
  expcones <- a*fir$DBH^b
  -sum(dnbinom(fir$TOTCONES,mu=expcones,size=k,log=TRUE))
}

params <- c(a=1,b.n=1,b.w=1,k=1)  

NegBinomLik_constak(params)


### Fit the new model

MLE_constak <- optim(fn=NegBinomLik_constak,par=params)

MLE_constak$par

MLE_constak$value

ms_constak <- 2*MLE_constak$value


###########
### Now, let's build and fit one more final model- this time with no wave effect and a Poisson error distribution

PoisLik_nowave <- function(params){
  a <- params[1]      # a parameters
  b <- params[2]      # b parameter (not a function of wave/nonwave)
  expcones <- a*fir$DBH^b
  -sum(dpois(fir$TOTCONES,lambda=expcones,log=TRUE))
}

params <- c(a=1,b=1)

PoisLik_nowave(params)

MLE_pois <- optim(fn=PoisLik_nowave,par=params)

MLE_pois$par

MLE_pois$value

ms_pois <- 2*MLE_pois$value


###########
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

###########
# Bayes factor example
###########

##### take a basic binomial distribution with parameter p fixed at 0.5:

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability")

## Q: What is the *marginal likelihood* under this simple model for an observation of 2 mortalities out of 10? 

## A:

dbinom(2,10,0.5)



## Now we can consider a model whereby "p" is a free parameter
curve(dbeta(x,1,1))  # uniform prior on "p"


###########
# Compute the marginal likelihood of observing 2 mortalities

# ?integrate
binom2 <- function(x) dbinom(x=2,size=10,prob=x)
marginal_likelihood <- integrate(f=binom2,0,1)$value    # use "integrate" function
marginal_likelihood  # equal to 0.0909 = 1/11


###########
# Compute the marginal likelihood of observing 3 mortalities

binom3 <- function(x) dbinom(x=3,size=10,prob=x)
marginal_likelihood <- integrate(f=binom3,0,1)$value    # use "integrate" function
marginal_likelihood   # equal to 0.0909 = 1/11


#########
# simulate data from the model across all possible values of the parameter "p"

lots=100000
hist(rbinom(lots,10,prob=rbeta(lots,1,1)),freq = F,xlab="potential observations")   # no particular observation is favored


#########
# Visualize the marginal likelihood of all possible observations

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))


###########
# Overlay the marginal likelihood of the simpler model, with p fixed at 0.5

probs2 <- rep(1/11,times=11)          
names(probs2) = 0:10
barplot(probs2,ylab="probability",ylim=c(0,1))

probs1 <- dbinom(0:10,10,0.5)          
names(probs1) = 0:10
barplot(probs1,ylab="probability",add=T,col="red",density=20)


############
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


############
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


#############
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


#########
# LRT: simple model (p fixed at 0.5) vs complex model (p is free parameter)

Likelihood_simple <- dbinom(3,10,0.5)
Likelihood_complex <- dbinom(3,10,0.3)
Likelihood_simple
Likelihood_complex
-2*log(Likelihood_simple)--2*log(Likelihood_complex)

qchisq(0.95,1)

pchisq(1.64,1)    # very high p value


#########
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


######
# Alternatively, try BIC

BIC_simple <- -2*log(Likelihood_simple) + 0
BIC_complex <-  -2*log(Likelihood_complex) + log(3)*1

BIC_simple
BIC_complex    



##############
# Bayesian model selection: Bolker's fir dataset

cat("

model  {

### Likelihood

  for(i in 1:n.obs){
    expected.cones[i] <- a[wave[i]]*pow(DBH[i],b[wave[i]])   # power function: a*DBH^b
    p[i] <- r[wave[i]] / (r[wave[i]] + expected.cones[i])
    observed.cones[i] ~ dnegbin(p[i],r[wave[i]])
  }
  
  
  ### Priors
  for(j in 1:2){   # estimate separately for wave and non-wave
    a[j] ~ dunif(0.001,2)
    b[j] ~ dunif(0.5,4)
    r[j] ~ dunif(0.5,5)
  }
  
}
    
",file="BUGS_fir.txt")


#######
# Package the data for JAGS

data.package1 <- list(
  observed.cones = fir$TOTCONES,
  n.obs = nrow(fir),
  wave = as.numeric(fir$WAVE_NON),
  DBH = fir$DBH
)
#data.package


##########
# Make a function for generating initial guesses

init.generator1 <- function(){ list(
  a = runif(2, 0.2,0.5),
  b = runif(2, 2,3),
  r = runif(2, 1,2)
  
  )
}
init.generator1()



###########
# Run the model in JAGS

library(R2jags)    # load packages
library(coda)
library(lattice)

params.to.monitor <- c("a","b","r")

jags.fit1 <- jags(data=data.package1,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=10000,model.file="BUGS_fir.txt",n.chains = 2,n.burnin = 2000,n.thin=5 )

jagsfit1.mcmc <- as.mcmc(jags.fit1)   # convert to "MCMC" object (coda package)

summary(jagsfit1.mcmc)

plot(jagsfit1.mcmc)



########
# Visualize the model fit

densityplot(jagsfit1.mcmc)


hist(jags.fit1$BUGSoutput$sims.list$r[,1],main="dispersion param",ylab="Prob Density",xlab="dispersion param",freq = F,ylim=c(0,2),xlim=c(0.5,2.5))
hist(jags.fit1$BUGSoutput$sims.list$r[,2],density=20,col="green",add=T,freq=F)
legend("topright",col=c("green","white"),density=c(20,0),legend=c("wave","nonwave"),bty="n")


#######
# Extract the DIC for the full model!

DIC_full <- jags.fit1$BUGSoutput$DIC
DIC_full


#################
# Build JAGS code for the reduced model

cat("

model  {

### Likelihood

  for(i in 1:n.obs){
    expected.cones[i] <- a*pow(DBH[i],b)   # a*DBH^b
    p[i] <- r / (r + expected.cones[i])
    observed.cones[i] ~ dnegbin(p[i],r)
  }
  
  
  ### Priors
  
  a ~ dunif(0.001,2)
  b ~ dunif(0.5,4)
  r ~ dunif(0.5,5)

  
}
    
",file="BUGS_fir_reduced.txt")


##########
# Package data for JAGS

data.package2 <- list(
  observed.cones = fir$TOTCONES,
  n.obs = nrow(fir),
  #wave = as.numeric(fir$WAVE_NON),
  DBH = fir$DBH
)


############
# Function for generating initial guesses for all params

init.generator2 <- function(){ list(
  a = runif(1, 0.2,0.5),
  b = runif(1, 2,3),
  r = runif(1, 1,2)
  
  )
}
init.generator2()


###########
# Run the reduced model and visualize the JAGS fit

params.to.monitor <- c("a","b","r")

jags.fit2 <- jags(data=data.package2,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=10000,model.file="BUGS_fir_reduced.txt",n.chains = 2,n.burnin = 2000,n.thin=5 )

jagsfit2.mcmc <- as.mcmc(jags.fit2)   # convert to "MCMC" object (coda package)

summary(jagsfit2.mcmc)

plot(jagsfit2.mcmc[,"a"])
plot(jagsfit2.mcmc[,"b"])
plot(jagsfit2.mcmc[,"r"])



densityplot(jagsfit2.mcmc)


########
# Compute DIC

DIC_reduced <- jags.fit2$BUGSoutput$DIC

DIC_reduced
DIC_full


#############
# Use WAIC for bayesian model selection!

library(loo)    # load the "loo" package, which allows us to compute WAIC from JAGS output'


####
# First, re-make the JAGS code, this time recording the likelihood as a derived parameter

cat("

model  {

### Likelihood

  for(i in 1:n.obs){
    expected.cones[i] <- a[wave[i]]*pow(DBH[i],b[wave[i]])   # power function: a*DBH^b
    p[i] <- r[wave[i]] / (r[wave[i]] + expected.cones[i])
    observed.cones[i] ~ dnegbin(p[i],r[wave[i]])
    LogLik[i] <- log(dnegbin(observed.cones[i],p[i],r[wave[i]]))   # add log likelihood computation for each observation!
  }
  
  
  ### Priors
  for(j in 1:2){   # estimate separately for wave and non-wave
    a[j] ~ dunif(0.001,2)
    b[j] ~ dunif(0.5,4)
    r[j] ~ dunif(0.5,5)
  }
  
}
    
",file="BUGS_fir.txt")


#################
# Build JAGS code for the reduced model

cat("

model  {

### Likelihood

  for(i in 1:n.obs){
    expected.cones[i] <- a*pow(DBH[i],b)   # a*DBH^b
    p[i] <- r / (r + expected.cones[i])
    observed.cones[i] ~ dnegbin(p[i],r)
    LogLik[i] <- log(dnegbin(observed.cones[i],p[i],r))   # add log likelihood computation for each observation!
  }
  
  
  ### Priors
  
  a ~ dunif(0.001,2)
  b ~ dunif(0.5,4)
  r ~ dunif(0.5,5)

  
}
    
",file="BUGS_fir_reduced.txt")


############
# re-fit the models

params.to.monitor <- c("a","b","r","LogLik")    # now monitor the log likelihood

jags.fit1 <- jags(data=data.package1,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=10000,model.file="BUGS_fir.txt",n.chains = 2,n.burnin = 2000,n.thin=5 )

jags.fit2 <- jags(data=data.package2,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=10000,model.file="BUGS_fir_reduced.txt",n.chains = 2,n.burnin = 2000,n.thin=5 )



############
# Compute WAIC!

loglik_full <- jags.fit1$BUGSoutput$sims.list$LogLik
loglik_red <- jags.fit2$BUGSoutput$sims.list$LogLik

waic_full <- waic(loglik_full)
waic_red <- waic(loglik_red)

waic_full$waic
waic_red$waic

compare(waic_full, waic_red)


#############
# Explicit Bayesian model selection

cat("

model  {

  ### Likelihood for model 1: full

  for(i in 1:n.obs){
    expected.cones[i,1] <- a1[wave[i]]*pow(DBH[i],b1[wave[i]])       # a*DBH^b
    spread.cones[i,1] <- r1[wave[i]]
    p[i,1] <- spread.cones[i,1] / (spread.cones[i,1] + expected.cones[i,1])
    observed.cones[i,1] ~ dnegbin(p[i,1],spread.cones[i,1])
    predicted.cones[i,1] ~ dnegbin(p[i,1],spread.cones[i,1])
    SE_obs[i,1] <- pow(observed.cones[i,1]-expected.cones[i,1],2)
    SE_pred[i,1] <- pow(predicted.cones[i,1]-expected.cones[i,1],2)
  }
  
  
  ### Priors, model 1
  for(j in 1:2){   # estimate separately for wave and non-wave
    a1[j] ~ dunif(0.001,2)
    b1[j] ~ dunif(0.5,4)
    r1[j] ~ dunif(0.5,5)
  }

  ### Likelihood for model 2: reduced

  for(i in 1:n.obs){
    expected.cones[i,2] <- a2*pow(DBH[i],b2)       # a*DBH^b
    spread.cones[i,2] <- r2
    p[i,2] <- spread.cones[i,2] / (spread.cones[i,2] + expected.cones[i,2])
    observed.cones[i,2] ~ dnegbin(p[i,2],spread.cones[i,2])
    predicted.cones[i,2] ~ dnegbin(p[i,2],spread.cones[i,2])
    SE_obs[i,2] <- pow(observed.cones[i,2]-expected.cones[i,2],2)
    SE_pred[i,2] <- pow(predicted.cones[i,2]-expected.cones[i,2],2)
  }
  
  
  ### Priors, model 2
  a2 ~ dunif(0.001,2)
  b2 ~ dunif(0.5,4)
  r2 ~ dunif(0.5,5)

  ### Likelihood for model 3: constant a and b

  for(i in 1:n.obs){
    expected.cones[i,3] <- a3*pow(DBH[i],b3)       # a*DBH^b
    spread.cones[i,3] <- r3[wave[i]]
    p[i,3] <- spread.cones[i,3] / (spread.cones[i,3] + expected.cones[i,3])
    observed.cones[i,3] ~ dnegbin(p[i,3],spread.cones[i,3])
    predicted.cones[i,3] ~ dnegbin(p[i,3],spread.cones[i,3])
    SE_obs[i,3] <- pow(observed.cones[i,3]-expected.cones[i,3],2)
    SE_pred[i,3] <- pow(predicted.cones[i,3]-expected.cones[i,3],2)
  }
  
  SSE_obs[1] <- sum(SE_obs[,1]) 
  SSE_pred[1] <- sum(SE_pred[,1])
  SSE_obs[2] <- sum(SE_obs[,2]) 
  SSE_pred[2] <- sum(SE_pred[,2])
  SSE_obs[3] <- sum(SE_obs[,3]) 
  SSE_pred[3] <- sum(SE_pred[,3])

  ### Priors, model 3
  for(j in 1:2){   # estimate separately for wave and non-wave
    r3[j] ~ dunif(0.5,5)
  }
  a3 ~ dunif(0.001,2)
  b3 ~ dunif(0.5,4)

  #####################
  ### SELECT THE BEST MODEL!!! 
  #####################

  for(i in 1:n.obs){
    observed.cones2[i] ~ dnegbin(p[i,selected],spread.cones[i,selected])
    predicted.cones2[i] ~ dnegbin(p[i,selected],spread.cones[i,selected])     # for posterior predictive check!
    SE2_obs[i] <- pow(observed.cones2[i]-expected.cones[i,selected],2)
    SE2_pred[i] <- pow(predicted.cones2[i]-expected.cones[i,selected],2)
  }
  
  SSE2_obs <- sum(SE2_obs[])
  SSE2_pred <- sum(SE2_pred[])


  ### Priors
  
    # model selection...
  prior[1] <- 1/3
  prior[2] <- 1/3     # put substantially more weight because fewer parameters (there are more rigorous ways to do this!!)
  prior[3] <- 1/3
  selected ~ dcat(prior[])   
  
  
}
    
",file="BUGS_fir_modelselection.txt")



#########
# Package the data for JAGS

data.package3 <- list(
  observed.cones = matrix(rep(fir$TOTCONES,times=3),ncol=3,byrow=F),
  observed.cones2 = fir$TOTCONES,
  n.obs = nrow(fir),
  wave = as.numeric(fir$WAVE_NON),
  #n.models = 3,
  DBH = fir$DBH
)
#data.package


#########
# Run JAGS

params.to.monitor <- c("a1","b1","r1","a2","b2","r2","a3","b3","r3","selected","predicted.cones2","predicted.cones","SSE_obs","SSE_pred","SSE2_obs","SSE2_pred")

jags.fit3 <- jags(data=data.package3,parameters.to.save=params.to.monitor,n.iter=5000,model.file="BUGS_fir_modelselection.txt",n.chains = 2,n.burnin = 1000,n.thin=2 )

jagsfit3.mcmc <- as.mcmc(jags.fit3)   # convert to "MCMC" object (coda package)

BUGSlist <- as.data.frame(jags.fit3$BUGSoutput$sims.list)
#summary(jagsfit.mcmc)

#plot(jagsfit.mcmc)



##########
# Visualize the model fit

#plot(jagsfit.mcmc[,"selected"])

plot(jagsfit3.mcmc[,"a1[1]"])
plot(jagsfit3.mcmc[,"a1[2]"])
plot(jagsfit3.mcmc[,"a2"])
plot(jagsfit3.mcmc[,"a3"])

plot(jagsfit3.mcmc[,"r1[1]"])
plot(jagsfit3.mcmc[,"r1[2]"])
plot(jagsfit3.mcmc[,"r2"])
plot(jagsfit3.mcmc[,"r3[1]"])



##########
# Perform explicit model selection

n.iterations <- length(jags.fit3$BUGSoutput$sims.list$selected)
selected <- table(jags.fit3$BUGSoutput$sims.list$selected)
names(selected) <- c("Full model","No wave","Fixed a&b")
selected

barplot(selected/n.iterations,ylab="Degree of belief")


##########
# Goodness of fit

n.data <- length(fir$DBH)

plot(fir$TOTCONES~fir$DBH,ylim=c(0,900),cex=2)

for(d in 1:n.data){
  tofind <- sprintf("predicted.cones[%s,1]",d)
  model1 <- as.vector(jagsfit3.mcmc[,tofind])
  points(rep(fir$DBH[d],times=100),sample(model1[[1]],100),pch=20,col="gray",cex=0.4)
}


#########
# Perform posterior predictive check

plot(fir$TOTCONES~fir$DBH,ylim=c(0,900),cex=2)

for(d in 1:n.data){
  tofind <- sprintf("predicted.cones[%s,2]",d)
  model1 <- as.vector(jagsfit3.mcmc[,tofind])
  points(rep(fir$DBH[d],times=100),sample(model1[[1]],100),pch=20,col="gray",cex=0.4)
}


plot(fir$TOTCONES~fir$DBH,ylim=c(0,900),cex=2)

for(d in 1:n.data){
  tofind <- sprintf("predicted.cones[%s,3]",d)
  model1 <- as.vector(jagsfit3.mcmc[,tofind])
  points(rep(fir$DBH[d],times=100),sample(model1[[1]],100),pch=20,col="gray",cex=0.4)
}


###############
# Posterior Predictive Checks!

plot(as.vector(jagsfit3.mcmc[,"SSE_pred[1]"][[1]])~as.vector(jagsfit3.mcmc[,"SSE_obs[1]"][[1]]),xlab="SSE, real data",ylab="SSE, perfect data",main="Posterior Predictive Check")
abline(0,1,col="red")
p.value=length(which(as.vector(jagsfit3.mcmc[,"SSE_pred[1]"][[1]])>as.vector(jagsfit3.mcmc[,"SSE_obs[1]"][[1]])))/length(as.vector(jagsfit3.mcmc[,"SSE_pred[1]"][[1]]))
p.value 


plot(as.vector(jagsfit3.mcmc[,"SSE_pred[2]"][[1]])~as.vector(jagsfit3.mcmc[,"SSE_obs[2]"][[1]]),xlab="SSE, real data",ylab="SSE, perfect data",main="Posterior Predictive Check")
abline(0,1,col="red")
p.value=length(which(as.vector(jagsfit3.mcmc[,"SSE_pred[2]"][[1]])>as.vector(jagsfit3.mcmc[,"SSE_obs[2]"][[1]])))/length(as.vector(jagsfit3.mcmc[,"SSE_pred[2]"][[1]]))
p.value 


plot(as.vector(jagsfit3.mcmc[,"SSE_pred[3]"][[1]])~as.vector(jagsfit3.mcmc[,"SSE_obs[3]"][[1]]),xlab="SSE, real data",ylab="SSE, perfect data",main="Posterior Predictive Check")
abline(0,1,col="red")
p.value=length(which(as.vector(jagsfit3.mcmc[,"SSE_pred[3]"][[1]])>as.vector(jagsfit3.mcmc[,"SSE_obs[3]"][[1]])))/length(as.vector(jagsfit3.mcmc[,"SSE_pred[3]"][[1]]))
p.value   


plot(fir$TOTCONES~fir$DBH,ylim=c(0,900),cex=2)

for(d in 1:n.data){
  tofind <- sprintf("predicted.cones2[%s]",d)
  model1 <- as.vector(jagsfit3.mcmc[,tofind])
  points(rep(fir$DBH[d],times=100),sample(model1[[1]],100),pch=20,col="gray",cex=0.4)
}


###########
# Posterior predictive check with model-averaged model!

plot(as.vector(jagsfit3.mcmc[,"SSE2_pred"][[1]])~as.vector(jagsfit3.mcmc[,"SSE2_obs"][[1]]),xlab="SSE, real data",ylab="SSE, perfect data",main="Posterior Predictive Check")
abline(0,1,col="red")
p.value=length(which(as.vector(jagsfit3.mcmc[,"SSE2_pred"][[1]])>as.vector(jagsfit3.mcmc[,"SSE2_obs"][[1]])))/length(as.vector(jagsfit3.mcmc[,"SSE2_pred"][[1]]))
p.value 

