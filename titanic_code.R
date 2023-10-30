


library(titanic)

## load data

titanic <- titanic::titanic_train



## write jags code
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

psurv0[1] ~ dunif(0,1)  # flat prior from 0 to 1 on p scale
psurv0[2] ~ dunif(0,1)
psurv0[3] ~ dunif(0,1)

for(i in 1:3){    # convert to logit scale
   psurv0.l[i] <- log(psurv0[i]/(1-psurv0[i]))
}


## interpolate missing data

meanage ~ dnorm(0,.1)
sdage ~ dunif(0,5)
precage <- pow(sdage,-2)

for(p in 1:npassengers){
  age[p] ~ dnorm(meanage,precage)   # DATA NODE
}


} ",file=fn)



# package the data for jags

dat <- list(
  npassengers = nrow(titanic),
  class = titanic$Pclass,
  fare = (titanic$Fare-mean(titanic$Fare))/sd(titanic$Fare),
  age = (titanic$Age-mean(titanic$Age,na.rm = T))/sd(titanic$Age,na.rm=T),
  is.fem = ifelse(titanic$Sex=="female",1,0),
  survived = titanic$Survived
)
dat


# set inits

# params to store

params <- c("psurv0","b.fare","b.age","b.female","sdage","meanage")

library(jagsUI)

# ?jags

mod <- jags(dat, parameters.to.save=params, model.file=fn,
     n.chains=3, n.adapt=1000, n.iter=10000, n.burnin=5000, n.thin=2,
     parallel=T)

mod

sims <- mod$sims.list
samps <- mod$samples

plot(mod,"psurv0")
plot(mod,"b.fare")
plot(mod, "b.female")
plot(mod, "b.age")


library(coda)
gelman.diag(samps)   # gelman-rubin diagnostic



# visualize the effect of age and sex

maleprob <- plogis( qlogis(sims$psurv0[,2]) + sims$b.fare*0  + sims$b.age*0 + sims$b.female*0)
femprob <- plogis( qlogis(sims$psurv0[,2]) + sims$b.fare*0  + sims$b.age*0 + sims$b.female*1)

df <- data.frame(
  sex <- c(rep("M",length(maleprob)),rep("F",length(femprob)) ),
  probsurv <- c(maleprob,femprob)
)

library(ggplot2)

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

library(tidyverse)

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



### posterior predictive check (bayes p-val)

# goal: summarize the observed vs expected number of mortalities per class and sex

# observed

titanic$died <- 1-titanic$Survived

obsdat <- with(titanic, tapply(died,list(Sex,Pclass), sum   ) )
nandx <- which(is.na(dat$age))

thissim <- 1

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

obsdat
expdat(2)


simdat <- function(thissim){
  dat$age[nandx] <- rnorm(length(nandx),sims$meanage[thissim],sims$sdage[thissim])
  theseps <- sapply(1:nrow(titanic), function(t) plogis(qlogis(sims$psurv0[thissim,titanic$Pclass[t]]) + sims$b.fare[thissim] * dat$fare[t]  
                                                        + sims$b.age[thissim] * dat$age[t] + sims$b.female[thissim] * dat$is.fem[t] ) )
  
  thismort <- rbinom(nrow(titanic),1,1-theseps)
  thissim <- tapply(thismort,list(titanic$Sex,titanic$Pclass), sum   )
  return(thissim)
}

chisq(obs=obsdat,exp=expdat(2))
chisq(obs=simdat(2),exp=expdat(2))

chisq(obs=obsdat,exp=expdat(3))
chisq(obs=simdat(3),exp=expdat(3))

nsims=100
errors_obs <- numeric(nsims)
errors_sim <- numeric(nsims)

for(i in 1:nsims){
  thissim <- sample(1:length(sims$b.fare),1)
  errors_obs[i] <- chisq(obs=obsdat,exp=expdat(thissim))
  errors_sim[i] <- chisq(obs=simdat(thissim),exp=expdat(thissim))
}

plot(errors_obs,errors_sim,xlim=c(0,30),ylim=c(0,30))
abline(1,1)




# next: try mtcars example with posterior predictive check?






