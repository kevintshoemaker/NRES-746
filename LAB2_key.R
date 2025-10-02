
#  NRES 746, Lab 2               
#   University of Nevada, Reno                      
#   Data generating models  ------------------------------      



# exercise 2 ---------

# exercise 2.1a

LMDataSim <- function(N,xlims,parms){
  xvals <- runif(N,xlims[1],xlims[2])
  data.frame(
    x = xvals,
    y = rnorm(N, parms[1] + parms[2]*xvals, parms[3]) 
  )
}


LMDataSim(10,xlims=c(40,120),parms=c(55,-0.09,6.8))


temp <- LMDataSim(10,xlims=c(40,120),parms=c(55,-0.09,6.8))

coef(lm(y~x,data=temp))[2]


## exercise 2.1b
# parms = c(55,-0.09,6.8)
# simdat=LMDataSim(5,xlims=c(40,120),parms=parms)
LMVerify <- function(simdat,parms,plot=T){
  mod <- lm(y~x,data=simdat)
  out = list()
  out$fitted_parms <- c(coef(mod),sigma=summary(mod)$sigma)
  ci = confint(mod,level = 0.9)
  out$slope_in <- dplyr::between(parms[2],ci[2,1],ci[2,2] )
  if(plot==T){
    a = ggplot(simdat,aes(x=x,y=y)) +
      geom_smooth(method="lm",level=0.9) +
      geom_point() +
      geom_abline(intercept=parms[1],slope=parms[2], col="darkgreen",lty=2) +
      theme_classic()
    print(a)  
  }
  return(out)
}


trueparms=c(55,-0.09,6.8)
df <- LMDataSim(15,xlims=c(40,120),parms=trueparms)

LMVerify(df, trueparms, plot=T)


simdat2 <- LMDataSim(N=100,xlims=c(0,15),parms=c(40,-2.2,1.8) )

LMVerify(simdat2, c(40,-2.2,1.8), plot=T)


## exercise 2.1c

lots=1000
trueparams <- c(55,-0.09,6.8)

inside = replicate(lots, {
  d=LMDataSim(10,xlims=c(40,120),parms=trueparams)
  LMVerify(d,trueparams,plot=FALSE)$slope_in 
})
sum(inside)/lots     #should be close to 90%



# Exercise 2 -------

## exercise 2.2a

# N=100;xlims=c(0,10);parms=c(0,5,2)
LMDataSim2 <- function(N,xlims,parms){
  xvals <-runif(N,xlims[1],xlims[2])
  data.frame(
    x=xvals,
    y=rnorm(N,parms[1]+parms[2]*xvals,parms[3]*abs(parms[1]+parms[2]*xvals))
  )
}


test <- LMDataSim2(500,xlims=c(1,10),parms=c(1,4,.5))
plot(test)


simdat <- LMDataSim2(100,xlims=c(1,10),parms=c(3,3,0.3))
LMVerify(simdat, c(3,3,0.3), plot=T)


##  Exercise 2.2b

### pt 1.

lots=1000
difs <- numeric(lots)
for(i in 1:lots){
  simdat <- LMDataSim2(100,xlims=c(1,10),parms=c(3,3,0.3))
  temp <- LMVerify(simdat, c(3,3,0.3),plot=F)
  difs[i] <- temp$fitted_parms[2]-3
}

hist(difs,xlab="difs between est and true slope")
mean(difs)  ## this test seems relatively robust to violations of homoskedasticity in this case!


##  Exercise 2.2c

## exercise 2.1c

lots=1000
trueparams <- c(3,3,0.3)

inside = replicate(lots, {
  d=LMDataSim2(10,xlims=c(1,10),parms=trueparams)
  LMVerify(d,trueparams,plot=FALSE)$slope_in 
})
sum(inside)/lots     #should be close to 90%



##  goodness of fit visualization...

simdat_uv <- LMDataSim2(500,xlims=c(5,120),parms=c(5,5,.25))
mod = lm(y~x, data=simdat_uv)
nd = data.frame(
  x=seq(min(simdat_uv$x),max(simdat_uv$x),length=100)
)
pred = predict(mod,nd,interval="prediction",level=0.89)
plot(simdat_uv)
lines(nd$x, pred[,"fit"],lwd=2)
lines(nd$x, pred[,"lwr"],lwd=1,lty=2)
lines(nd$x, pred[,"upr"],lwd=1,lty=2)
legend("topleft",lty=2,legend="89% prediction interval")


# Exercise 2.3: power analysis -----------


    ##  first, read in the functions from the "virtual ecologist" lecture. 

## helper function 1 ---------------------------
# function for computing the probability of observing each animal in a single multi-day, multi-observer survey

    # Arguments:
      # N: population abundance
      # people: number of survey participants each day
      # days: survey duration, in days
      # pObs: prob of detection per person per day

SurvProb <- function(N=1000,people=1,days=3,pObs=0.02){
  probPerDay <- 1-(1-pObs)^people      # define the probability of detection per animal per day
  1-(1-probPerDay)^days       # define the probability of detection per animal for the entire survey
}

## function: simulate monitoring data ----------------------------

# develop a function for simulating monitoring data from a declining population

    # Arguments:
      # N0: true initial population abundance
      # trend: proportional change in population size from last year
      # nyears: duration of simulation
      # people: number of survey participants each day
      # days: survey duration, in days
      # survint: survey interval, in years (e.g., 2 means surveys are conducted every other year)

SimulateMonitoringData <- function(N0=1000,trend=-0.03,nyears=25,people=1,days=3,survint=2){
  realAbund <- floor(N0 * (1+trend)^(0:(nyears-1)))   # no fractional individuals
  detected <- sapply(realAbund, function(t)  rbinom(1,t, SurvProb(t,people=people,days=days,pObs=0.02)   ) )
  detected[!c(0:(nyears-1))%%survint==0] <- NA            # if the survey is not performed that year, return a missing value
  detected       # return the number of individuals detected
}

## function: assessing whether or not a decline was detected ------------------------

    # Arguments:
      # monitoringData: simulated results from a long-term monitoring study
      # alpha: define acceptable type-I error rate (acceptable false positive rate)

IsDecline <- function(monitoringData,alpha=0.05){
  time <- 1:length(monitoringData)      # vector of survey years
  model <- lm(monitoringData~time)    # for now, let's use ordinary linear regression (perform linear regression on simulated monitoring data)
  p_value <- summary(model)$coefficients["time","Pr(>|t|)"]      # extract the p-value  
  isdecline <- ifelse(summary(model)$coefficients["time","Estimate"]<0,TRUE,FALSE)     # determine if the simulated monitoring data determined a "significant" decline
  sig_decline <- ifelse((p_value<=alpha)&(isdecline),TRUE,FALSE)    # if declining and significant trend, then the monitoring protocol successfully diagnosed a decline
  return(sig_decline)
}



# Lab exercise 2.3a: solution

GetPower <- function(N0=1000,trend=-0.03,nyears=25,people=1,days=3,survint=2,alpha=0.05){
  lots=1000
  dats <- replicate(lots, SimulateMonitoringData(N0,trend,nyears,people,days,survint) )
  tests <- apply(dats,2,IsDecline,alpha=alpha)
  sum(tests==TRUE)/lots 
}


initabund = 1000
  survints <- c(1:5)
  powers <- numeric(length(survints))
  for(i in 1:length(survints)){
    powers[i] <- GetPower(survint=survints[i])
  }
  
  plot(powers~survints,xlab="Survey interval(years)",ylab="Statistical Power",main="Power to detect trend, by sampling interval") 


# exercise 2.3b

########  survey intervals 

survints <- c(1:5)
powers <- numeric(length(survints))
for(i in 1:length(survints)){
  powers[i] <- GetPower(survint=survints[i])
}

plot(powers~survints,xlab="Survey interval(years)",ylab="Statistical Power",main="Power to detect trend, by sampling interval") 


#######  number of observers

observers <- c(1:5)
powers <- numeric(length(observers))
for(i in 1:length(observers)){
  powers[i] <- GetPower(people=observers[i])
}
plot(powers~observers,xlab="Observers",ylab="Statistical Power",main="") 


####### days per survey bout

days <- c(1:5)
powers <- numeric(length(days))
for(i in 1:length(days)){
  powers[i] <- GetPower(days=days[i])
}
plot(powers~days,xlab="Survey bout length (days)",ylab="Statistical Power",main="") 



########
# exercise 2.3c

## there are many correct answers here!  Here is one way to evaluate multiple scenarios:


scenarios <- expand.grid(survints,observers,days)
names(scenarios) = c("survints", "observers","days")
powers <- numeric(nrow(scenarios))
costs <- numeric(nrow(scenarios))
i=1
for(i in 1:nrow(scenarios)){
  nyears <- 25
  powers[i] <- GetPower(people=scenarios$observers[i],days=scenarios$days[i],survint=scenarios$survints[i])
  nsurveys <- floor(nyears/scenarios$survint[i])
  costs[i] <- nsurveys*2000 + nsurveys*scenarios$observers[i]*scenarios$days[i]*200
}

successful <- which(powers>=0.75)

leastcost <- which.min(costs[successful])

bestscenario <- successful[leastcost]   # find the successful scenario with minimum cost...

scenarios[bestscenario,]    # identify the params for the best scenario!



