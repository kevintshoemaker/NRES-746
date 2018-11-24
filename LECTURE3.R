
############################################################
####                                                    ####  
####  NRES 746, Lecture 3                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  The virtual ecologist                             ####
####     Building data simulation models                ####
############################################################



###########
# Random number generators (a key component of data simulation models- but usually not the whole story)

runif(1,0,25)   # draw random numbers from various probability distributions
rpois(1,3.4)
rnorm(1,22,5.4)


#############
# Short exercise:

# Generate 50 samples from Normal(mean=10,sd=5) 


# Generate 1000 samples from Poisson(mean=50)


# Generate 10 samples from Beta(shape1=0.1,shape2=0.1)



############
# SIMULATE DATA GENERATION: decompose into deterministic and stochastic components 

##########
# Deterministic component: define function for transforming a predictor variable into an expected response (linear regression)

    # Arguments:
      # x: vector of covariate values
      # a: the intercept of a linear relationship mapping the covariate to an expected response
      # b: the slope of a linear relationship mapping the covariate to an expected response

deterministic_component <- function(x,a,b){
  linear <- a + b*x   # specify a deterministic, linear functional form
  return(linear)
}

xvals = seq(0,100,10)  # define the values of a hypothetical predictor variable (e.g., tree girth)

expected_vals <- deterministic_component(xvals,175,-1.5)   # use the deterministic component to determine the expected response (e.g., tree volume)
expected_vals

plot(xvals,expected_vals)   # plot out the relationship


##########
# Stochastic component: define a function for transforming an expected (deterministic) response and adding a layer of "noise" on top!

    # Arguments:
      # x: vector of expected responses
      # variance: variance of the "noise" component of your data simulation model
stochastic_component <- function(x,variance){     
  sd <- sqrt(variance)       # convert variance to standard deviation       
  stochvals <- rnorm(length(x),x,sd)       # add a layer of "noise" on top of the expected response values
  return(stochvals)
}

    ### Simulate stochastic data!!
sim_vals <- stochastic_component(expected_vals,variance=500)   # try it- run the function to add noise to your expected values. 

plot(xvals,sim_vals)     # plot it- it should look much more "noisy" now!

# ALTERNATIVELY:

sim_vals <- stochastic_component(deterministic_component(xvals,175,-1.5),500)    # stochastic "shell" surrounds a deterministic "core"    


############
# Goodness-of-fit test!

    # Does the data fall into the range of plausble data produced by this fully specified model?

############
# Imagine you have the following "real" data (e.g., tree volumes). 

realdata <- data.frame(Volume=c(125,50,90,110,80,75,100,400,350,290,350),Girth=xvals)
plot(realdata$Girth,realdata$Volume)


#############
# Let's simulate many datasets from our hypothesized data generating model (intercept=10,slope=4,variance=1000):

reps <- 1000    # specify number of replicate datasets to generate
samplesize <- nrow(realdata)    # define the number of data points we should generate for each simulation "experiment"
simresults <- array(0,dim=c(samplesize,reps))   # initialize a storage array for results 
for(i in 1:reps){       # for each independent simulation "experiment":
  exp_vals <- deterministic_component(realdata$Girth,a=10,b=4)          # simulate the expected tree volumes for each measured girth value
  sim_vals <- stochastic_component(exp_vals,1000)  # add stochastic noise
  simresults[,i] <- sim_vals   # store the simulated data for later
}

    # now make a boxplot of the results
boxplot(t(simresults),xaxt="n")    # (repeat) make a boxplot of the simulation results
axis(1,at=c(1:samplesize),labels=realdata$Girth)                          # add x axis labels 


#########
# Now overlay the "real" data
    # how well does the model fit the data?

boxplot(lapply(1:nrow(simresults), function(i) simresults[i,]),xaxt="n")    # (repeat) make a boxplot of the simulation results
axis(1,at=c(1:samplesize),labels=realdata$Girth)                          # add x axis labels 
points(c(1:samplesize),realdata$Volume,pch=20,cex=3,col="red",xaxt="n")     # this time, overlay the "real" data 


###########
# Using data simulation to flesh out sampling distributions for frequentist inference 

    # e.g., the "brute force t test" example:

reps <- 1000                    # number of replicate samples to generate
null_difs <- numeric(reps)              # storage vector for the test statistic for each sample
for(i in 1:reps){
  sampleA <- rnorm(10,10,4)        # sample representing "groups" A and B under the null hypothesis
  sampleB <- rnorm(10,10,4)
  null_difs[i] <- mean(sampleA)-mean(sampleB)         # test statistic (model result)
}

hist(null_difs)           # plot out the sampling distribution
abline(v=3.5,col="green",lwd=3)


###############
# Power analysis example: designing a monitoring program for a rare species


   ### first, let's develop some helper functions:

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
NumObserved(TrueN=500,surveyors=2,days=7)   # test the new function


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
ThisYearAbund(LastYearAbund=500,trend=-0.03)    # test the new function

# NOTE: we could introduce stochastic population dynamics (or density dependence, etc!) for a more realistic model, but we are omitting this here. 


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

SimulateMonitoringData(initabund=1000,trend=-0.03,years=25,observers=1,days=3,survint=2)    # test the new function


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

IsDecline(monitoringData=c(10,20,NA,15,1),alpha=0.05)    # test the function


###########
# Lab exercise: develop a "power" function to return the statistical power to detect a decline under alternative monitoring schemes...

nreps <- 10000      # set number of replicate monitoring "experiments"
initabund <- 1000    # set initial population abundance.

GetPower <- function(nreps=nreps,initabund=initabund,trend=-0.03,years=25,observers=1,days=3,survint=2,alpha=0.05){
     # fill this in!
  return(Power)
}

