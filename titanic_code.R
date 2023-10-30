
# TITANIC example ---------------------

## clear workspace ---------------------

rm(list=ls())

## load packages ----------------------

library(titanic)
library(ggplot2)
library(tidyverse)
library(emdbook)
library(coda)

## load data --------------------

titanic <- titanic::titanic_train

## write jags code ------------------- 

fn <- "titanic_jags.txt"
cat("

model{

# likelihood

for(p in 1:npassengers){
  logit(psurv[p]) <- psurv0.l[class[p]] + b.fare * fare[p]  + b.age * age[p] + b.female * is.fem[p]
  survived[p] ~ dbern(psurv[p])
}

# priors

b.fare ~ dnorm(0,0.1)  # slightly regularized prior
b.age ~ dnorm(0,0.1)
b.female ~ dnorm(0,0.1)

for(i in 1:3){    
   psurv0[i] ~ dunif(0,1)                       # flat prior from 0 to 1 on p scale
   psurv0.l[i] <- log(psurv0[i]/(1-psurv0[i]))  # convert to logit scale
}


## interpolate missing data

meanage ~ dnorm(0,.1)
sdage ~ dunif(0,5)
precage <- pow(sdage,-2)

for(p in 1:npassengers){
  age[p] ~ dnorm(meanage,precage)   # DATA NODE
}


} ",file=fn)



## package the data for jags --------------------


# prep for standardization 

meanfare <- mean(titanic$Fare)
sdfare <- sd(titanic$Fare)
meanage <- mean(titanic$Age,na.rm=T)
sdage <- sd(titanic$Age,na.rm=T)

dat <- list(
  npassengers = nrow(titanic),
  class = titanic$Pclass,
  fare = (titanic$Fare-meanfare)/sdfare,
  age = (titanic$Age-meanage)/sdage,
  is.fem = ifelse(titanic$Sex=="female",1,0),
  survived = titanic$Survived
)
dat

## set inits -----------------  

inits <- function(){
  list(
    meanage = runif(1,-0.5,0.5),
    sdage = runif(1,1,1.5),
    b.fare = runif(1,-0.5,0.5),
    b.age = runif(1,-0.5,0.5),
    b.female = runif(1,-0.5,0.5),
    psurv0 = runif(3,0.4,0.6)   # note: three random numbers!
  )
}
inits()

# params to store

params <- c("psurv0","b.fare","b.age","b.female","sdage","meanage")

library(jagsUI)


## run jags -----------------------

# ?jags

mod <- jags(dat, parameters.to.save=params, model.file=fn,
     n.chains=3, n.adapt=1000, n.iter=10000, n.burnin=5000, n.thin=2,
     parallel=T)

mod   # check Rhat (same as psrf)

sims <- mod$sims.list
samps <- mod$samples

plot(mod,"psurv0")
plot(mod,"b.fare")
plot(mod, "b.female")
plot(mod, "b.age")

gelman.diag(samps)   # gelman-rubin diagnostic

## visualize the effect of age and sex  ------------------

maleprob <- plogis( qlogis(sims$psurv0[,2]) + sims$b.fare*0  + sims$b.age*0 + sims$b.female*0)
femprob <- plogis( qlogis(sims$psurv0[,2]) + sims$b.fare*0  + sims$b.age*0 + sims$b.female*1)

df <- data.frame(
  sex <- c(rep("M",length(maleprob)),rep("F",length(femprob)) ),
  probsurv <- c(maleprob,femprob)
)

ggplot(df,aes(sex,probsurv)) +
  geom_violin()

# and age...

age.st <- seq(min(dat$age,na.rm=T),max(dat$age,na.rm=T),length=20)
age.real <- age.st*sd(titanic$Age,na.rm=T) + mean(titanic$Age,na.rm = T)
df2 <- NULL

i=1
for(i in 1:length(age.st)){
  thisa <- age.st[i]
  thisa_real <- age.real[i]
  thisp <- plogis(qlogis(sims$psurv[,2]) + sims$b.fare * 0  + sims$b.age * thisa + 
                      sims$b.female * 0)
  temp <- data.frame(
    age = thisa_real,
    probsurv = thisp 
  )
  df2 <- rbind(df2,temp)
}

df3 <- df2 %>% 
  group_by(age) %>% 
  summarize(
    lwr = quantile(probsurv,0.025),
    med = quantile(probsurv,0.5),
    upr = quantile(probsurv,0.975)
  )

ggplot(df3,aes(age,med)) +
  geom_ribbon(aes(ymin=lwr,ymax=upr)) +
  geom_path() +
  labs(y="prob surv")


## posterior predictive check (bayes p-val)  ---------------------------

# goal: summarize the observed vs expected number of mortalities per class and sex

# observed

titanic$died <- 1-titanic$Survived

obsdat <- with(titanic, tapply(died,list(Sex,Pclass), sum   ) )
nandx <- which(is.na(dat$age))

expdat <- function(thissim){
  dat$age[nandx] <- rnorm(length(nandx),sims$meanage[thissim],sims$sdage[thissim])
  theseps <- sapply(1:nrow(titanic), function(t) plogis(qlogis(sims$psurv0[thissim,titanic$Pclass[t]]) + sims$b.fare[thissim] * dat$fare[t]  
              + sims$b.age[thissim] * dat$age[t] + sims$b.female[thissim] * dat$is.fem[t] ) )

  thisexp <- tapply(1-theseps,list(titanic$Sex,titanic$Pclass), sum   )
  return(thisexp)
}


chisq <- function(obs,exp){
  sum((obs-exp)^2/exp ) 
}

simdat <- function(thissim){
  dat$age[nandx] <- rnorm(length(nandx),sims$meanage[thissim],sims$sdage[thissim])
  theseps <- sapply(1:nrow(titanic), function(t) plogis(qlogis(sims$psurv0[thissim,titanic$Pclass[t]]) + sims$b.fare[thissim] * dat$fare[t]  
                                                        + sims$b.age[thissim] * dat$age[t] + sims$b.female[thissim] * dat$is.fem[t] ) )
  
  thismort <- rbinom(nrow(titanic),1,1-theseps)
  thissim <- tapply(thismort,list(titanic$Sex,titanic$Pclass), sum   )
  return(thissim)
}

obsdat
expdat(3)

chisq(obs=obsdat,exp=expdat(2))
chisq(obs=simdat(2),exp=expdat(2))

chisq(obs=obsdat,exp=expdat(3))
chisq(obs=simdat(3),exp=expdat(3))

nsims=300
errors_obs <- numeric(nsims)
errors_sim <- numeric(nsims)

for(i in 1:nsims){
  thissim <- sample(1:length(sims$b.fare),1)
  errors_obs[i] <- chisq(obs=obsdat,exp=expdat(thissim))
  errors_sim[i] <- chisq(obs=simdat(thissim),exp=expdat(thissim))
}

plot(errors_obs,errors_sim,xlim=c(0,30),ylim=c(0,30))
abline(1,1)


# MTCARS example -----------------------

## clear workspace --------------

rm(list=ls())

## load data ----------------

data(mtcars)

## write jags code --------------

fn <- "mtcars_jags.txt"
cat("
model{

b0~dunif(10,40)
b1~dunif(-0.05,0)
mpg.sd~dunif(0,10)
mpg.prec<-pow(mpg.sd,-2)
for(i in 1:nobs){
  mpg.exp[i] <- b0 * exp(b1*disp[i])
  mpg[i] ~ dnorm(mpg.exp[i],mpg.prec)  # data node
}

}",file=fn)


## package data for jags  ---------------------

dat <- list(
  nobs=nrow(mtcars),
  disp=mtcars$disp,
  mpg=mtcars$mpg
)
dat

inits <- function(){
  list(
    b0=runif(1,30,35),
    b1=runif(1,-0.005,-0.001),
    mpg.sd=runif(1,2,3)
  )
}
inits()


pars <- c("b0","b1","mpg.sd")

## run jags ---------------

mod <- jags(dat, inits, pars, model.file=fn,
            n.chains=3, n.adapt=1000, n.iter=10000, n.burnin=5000, n.thin=2,
            parallel=T)

plot(mod,"b0")
plot(mod,"b1")
plot(mod,"mpg.sd")


sims = mod$sims.list   # joint posterior distribution (easy to use!)
samps <- mod$samples   # MCMC chains
nmcmc <- length(sims$b0)

## posterior predictive check -----------------

mcmc_ndx <- 1

expvals <- function(mcmc_ndx){
  sapply(1:nrow(mtcars), function(t) sims$b0[mcmc_ndx] * exp(sims$b1[mcmc_ndx]*mtcars$disp[t]) )
}

simdat <- function(mcmc_ndx){
  exp_mpg <- sapply(1:nrow(mtcars), function(t) sims$b0[mcmc_ndx] * exp(sims$b1[mcmc_ndx]*mtcars$disp[t]) )
  rnorm(nrow(mtcars), exp_mpg, sims$mpg.sd[mcmc_ndx])
}

expvals(100)
simdat(100)

rmse <- function(obs,exp){
  sqrt(sum((obs-exp)^2))
}

rmse(mtcars$mpg,expvals(100))
rmse(simdat(100),expvals(100))

nsims = 1000


## write for loop in class!  ------------------














