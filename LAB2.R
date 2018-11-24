
############################################################
####                                                    ####  
####  NRES 746, Lab 2                                   ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  The virtual ecologist                             ####
############################################################



#############
## exercise 1a

LMDataSim <- function(sample.size,xlimits,coefs,stdev){
  xvals <- runif(sample.size,xlimits[1],xlimits[2])
  yvals <- coefs[1] + coefs[2]*xvals
  yvals <- rnorm(sample.size,yvals,stdev)
  out <- data.frame(x=xvals,y=yvals)
  return(out)
}



LMDataSim(sample.size=10,xlimits=c(40,120),coefs=c(55,-0.09),stdev=6.8)


#############
## exercise 1b

LMVerify_pt1 <- function(sim.dataset,trueparams,plot=T){
  names <- names(sim.dataset)
  form <- as.formula(sprintf("%s ~ %s",names[2],names[1]) )   # note: this is to make sure the function runs for whatever your column names are!
  mod <- lm(form,data=sim.dataset)
  coefs <- coefficients(mod)
  sdest <- sd(residuals(mod))
  difs <- c(coefs,sdest)-trueparams
  names(difs) <- c("intercept",names[1],"stdev")
  if(plot==T){
    xl <- range(sim.dataset[,1])*c(0.75,1.25)
    yl <- range(sim.dataset[,2])*c(0.75,1.25)
    newdata <- data.frame(temp=seq(xl[1],xl[2],0.1))
    names(newdata) <- names[1]
    ci <- predict(mod,newdata,interval="confidence")
    plot(sim.dataset,ylim=yl,xlim=xl)
    abline(mod,col="green",lwd=3)
    lines(newdata[,1],ci[,"lwr"],col="blue",lty=2)
    lines(newdata[,1],ci[,"upr"],col="blue",lty=2)
  }
  cislope <- confint(mod,names[1])
  trueslope <- trueparams[2]
  is.inside <- ifelse((trueslope>=cislope[1])&(trueslope<=cislope[2]),TRUE,FALSE)
  out <- list()
  out$difs <- difs
  out$inside <- is.inside
  return(out)
}


simdat <- LMDataSim(sample.size=10,xlimits=c(40,120),coefs=c(55,-0.09),stdev=6.8)

LMVerify_pt1(sim.dataset = simdat, trueparams = c(55,-0.09,6.8))


#############
## exercise 1c

reps=1000
trueparams <- c(55,-0.09,6.8)
is.inside <- logical(reps)

i=1   # for debugging
for(i in 1:reps){
  simdat <- LMDataSim(sample.size=10,xlimits=c(40,120),coefs=trueparams[c(1:2)],stdev=trueparams[3])
  lmv <- LMVerify_pt1(sim.dataset = simdat, trueparams = trueparams,plot=FALSE)
  is.inside[i] <- lmv[[2]]
}

sum(is.inside)/reps



#############
## exercise 2b

LMDataSim2 <- function(sample.size,xlimits,coefs,sd_coefs){
  xvals <- runif(sample.size,xlimits[1],xlimits[2])
  yvals <- coefs[1] + coefs[2]*xvals
  sds <- sd_coefs[1] + sd_coefs[2]*xvals
  yvals <- rnorm(sample.size,yvals,sds)
  out <- data.frame(x=xvals,y=yvals)
  return(out)
}


test <- LMDataSim2(sample.size=35,xlimits=c(40,120),coefs=c(55,-0.09),sd_coefs=c(-0.2,0.04))
plot(test$x,test$y)


#############
## exercise 2c

LMVerify_pt2 <- function(sim.dataset,trueparams,plot=T){
  names <- names(sim.dataset)
  form <- as.formula(sprintf("%s ~ %s",names[2],names[1]) )   # note: this is to make sure the function runs for whatever your column names are!
  mod <- lm(form,data=sim.dataset)
  coefs <- coefficients(mod)
  difs <- coefs-trueparams
  names(difs) <- c("intercept",names[1])
  if(plot==T){
    xl <- range(sim.dataset[,1])*c(0.75,1.25)
    yl <- range(sim.dataset[,2])*c(0.75,1.25)
    newdata <- data.frame(temp=seq(xl[1],xl[2],0.1))
    names(newdata) <- names[1]
    ci <- predict(mod,newdata,interval="prediction")
    plot(sim.dataset,ylim=yl,xlim=xl)
    abline(mod,col="green",lwd=3)
    lines(newdata[,1],ci[,"lwr"],col="blue",lty=2)
    lines(newdata[,1],ci[,"upr"],col="blue",lty=2)
  }
  cislope <- confint(mod,names[1])
  trueslope <- trueparams[2]
  is.inside <- ifelse((trueslope>=cislope[1])&(trueslope<=cislope[2]),TRUE,FALSE)
  out <- list()
  out$difs <- difs
  out$inside <- is.inside
  return(out)
}


simdat <- LMDataSim2(sample.size=35,xlimits=c(40,120),coefs=c(55,-0.09),sd_coefs=c(-0.3,0.05))
LMVerify_pt2(sim.dataset = simdat, trueparams = c(55,-0.09))


##########
#  Exercise 2d

### pt 1.

reps=1000
difs <- numeric(reps)
i=1
for(i in 1:reps){
  simdat <- LMDataSim2(sample.size=20,xlimits=c(40,120),coefs=c(55,-0.09),sd_coefs=c(-0.3,0.05))
  temp <- LMVerify_pt2(sim.dataset = simdat, trueparams = c(55,-0.09),plot=F)
  difs[i] <- temp$difs["x"]
}

hist(difs,xlab="difs between est and true slope")
mean(difs)  ## this test seems relatively robust to violations of homoskedasticity in this case!


##  Part 2

simdat_lots <- LMDataSim2(sample.size=10000,xlimits=c(40,120),coefs=c(55,-0.09),sd_coefs=c(-0.3,0.05))
simdat_sample <- LMDataSim2(sample.size=10,xlimits=c(40,120),coefs=c(55,-0.09),sd_coefs=c(-0.3,0.05))
temp <- LMVerify_pt2(sim.dataset = simdat_sample, trueparams = c(55,-0.09),plot=T)
points(simdat_lots$x,simdat_lots$y,pch=19,cex=0.5,col=gray(0.5))



############
# Exercise 3: power analysis


####### first, read in the functions from the "virtual ecologist" lecture. 

########
# function for computing the number of observed/detected animals in a single survey

    # Arguments:
      # TrueN: true population abundance
      # surveyors: number of survey participants each day
      # days: survey duration, in days

NumObserved <- function(TrueN=1000,surveyors=1,days=3){
  probPerPersonDay <- 0.02      # define the probability of detection per animal per person-day [hard-coded- potentially bad coding practice!]
  probPerDay <- 1-(1-probPerPersonDay)^surveyors      # define the probability of detection per animal per day (multiple surveyors)(animal must be detected at least once)
  probPerSurvey <- 1-(1-probPerDay)^days       # define the probability of detection per animal for the entire survey
  nobs <- rbinom(1,size=TrueN,prob=probPerSurvey)     # simulate the number of animals detected!
  return(nobs)
}
# NumObserved(TrueN=500,surveyors=2,days=7)   # test the new function


#########
# function for computing expected abundance dynamics of a declining population (deterministic component!)

    # Arguments:
      # LastYearAbund: true population abundance in the previous year
      # trend: proportional change in population size from last year

ThisYearAbund <- function(LastYearAbund=1000,trend=-0.03){
  CurAbund <- LastYearAbund + trend*LastYearAbund    # compute abundance this year
  CurAbund <- floor(CurAbund)  # can't have fractional individuals!
  return(CurAbund)
}


########
# develop a function for simulating monitoring data from a declining population

    # Arguments:
      # initabund: true initial population abundance
      # trend: proportional change in population size from last year
      # years: duration of simulation
      # observers: number of survey participants each day
      # days: survey duration, in days
      # survint: survey interval, in years (e.g., 2 means surveys are conducted every other year)

SimulateMonitoringData <- function(initabund=1000,trend=-0.03,years=25,observers=1,days=3,survint=2){
  prevabund <- initabund        # initialize "previous-year abundance" at initial abundance 
  detected <- numeric(years)     # set up storage variable
  for(y in 1:years){            # for each year of the simulation:
    thisAbund <- ThisYearAbund(prevabund,trend)             # compute the current abundance on the basis of the trend
    detected[y] <- NumObserved(thisAbund,observers,days)     # sample the current population using this monitoring scheme
    prevabund <- thisAbund   # set this years abundance as the previous years abundance (to set up the simulation for next year)
  }
  surveyed <- c(1:years)%%survint==0    # which years were surveys actually performed?
  detected[!surveyed] <- NA            # if the survey is not performed that year, return a missing value
  return(detected)       # return the number of individuals detected
}


#########
# finally, develop a function for assessing whether or not a decline was detected:

    # Arguments:
      # monitoringData: simulated results from a long-term monitoring study
      # alpha: define acceptable type-I error rate (false positive rate)

IsDecline <- function(monitoringData,alpha=0.05){
  time <- 1:length(monitoringData)      # vector of survey years
  model <- lm(monitoringData~time)    # for now, let's use ordinary linear regression (perform linear regression on simulated monitoring data)
  p_value <- summary(model)$coefficients["time","Pr(>|t|)"]      # extract the p-value  
  isdecline <- ifelse(summary(model)$coefficients["time","Estimate"]<0,TRUE,FALSE)     # determine if the simulated monitoring data determined a "significant" decline
  sig_decline <- ifelse((p_value<=alpha)&(isdecline),TRUE,FALSE)    # if declining and significant trend, then the monitoring protocol successfully diagnosed a decline
  return(sig_decline)
}





###########
# Lab exercise 3a: develop a "power" function to return the statistical power to detect a decline under alternative monitoring schemes...

GetPower <- function(observers=1,days=3,alpha=0.05,years=25,survint=2,trend=-0.03,verbose=T){
  nreps <- 1000      # set number of replicate monitoring "experiments"
  initabund <- 1000    # set initial population abundance.
  decline <- logical(nreps)
  for(i in 1:nreps){
    detected <- SimulateMonitoringData(initabund,trend,years,observers,days,survint)
    decline[i] <- IsDecline(detected,alpha)
  }
  Power <- sum(decline)/nreps 
  if(verbose) cat((sprintf("The statistical power to detect a decline for the default parameters is: %s \n",Power)))    # for default
  return(Power)
}


GetPower(survint=3)         # test the new function


########
# exercise 3b

########  survey intervals

initabund = 1000
survints <- c(1:5)
powers <- numeric(length(survints))
for(i in 1:length(survints)){
  powers[i] <- GetPower(survint=survints[i],verbose=F)
}
  
plot(powers~survints,xlab="Survey interval(years)",ylab="Statistical Power",main="Power to detect trend, by sampling interval") 
  
#######  number of observers

initabund = 1000
observers <- c(1:10)
powers <- numeric(length(observers))
for(i in 1:length(observers)){
  powers[i] <- GetPower(observers=observers[i],days=1,survint=3,verbose=F)
}
  
plot(powers~observers,xlab="# observers",ylab="Statistical Power",main="Power to detect trend, by observers") 

####### days per survey bout  

initabund = 1000
days <- c(1:5)
powers <- numeric(length(days))
for(i in 1:length(days)){
  powers[i] <- GetPower(days=days[i],verbose=F)
}
  
plot(powers~days,xlab="days per bout",ylab="Statistical Power",main="Power to detect trend, by days per sampling bout") 
  


########
# exercise 3c

## there are many correct answers here!  Here is one way to evaluate multiple scenarios:


# GetPower(observers=1,days=3,alpha=0.05,years=25,survint=2,trend=-0.03,verbose=T)

scenarios <- data.frame(
  observers = rep(rep(c(1:5),times=5),times=4),
  days = rep(rep(1:5,each=5),times=4),
  survint = rep(4:1,each=25)
)
powers <- numeric(nrow(scenarios))
costs <- numeric(nrow(scenarios))

i=1
for(i in 1:nrow(scenarios)){
  nyears <- 25
  powers[i] <- GetPower(observers=scenarios$observers[i],days=scenarios$days[i],alpha=0.05,years=nyears,survint=scenarios$survint[i],trend=-0.03,verbose=F)
  nsurveys <- floor(nyears/scenarios$survint[i])
  costs[i] <- nsurveys*2000 + nsurveys*scenarios$observers[i]*scenarios$days[i]*200
}

successful <- which(powers>=0.75)

leastcost <- which.min(costs[successful])

bestscenario <- successful[leastcost]   # find the successful scenario with minimum cost... 

scenarios[bestscenario,]    # identify the params for the best scenario!    Here, it's 4 observers surveying for 1 day every 4 years!

##########
#  END LAB 2
##########



