
#  NRES 746, Lecture 3  ------------------------
##  University of Nevada, Reno       
##  Data generating models                             
##     Building data simulation models
##     Deterministic and stochastic processes


# Random number generators  -------------------------------

runif(1,0,25)   # draw random numbers from various probability distributions
rpois(1,3.4)
rnorm(1,22,5.4)


# Short exercise:

# Generate 50 samples from Normal(mean=10,sd=4.1) 


# Generate 1000 samples from Poisson(mean=5.4)


# Generate 10 samples from Beta(shape1=0.1,shape2=0.1)


# Try some other distributions and parameters.  NOTE: you can visualize probability densities easily using the "curve" function. or in ggplot:

library(ggplot2)
ggplot() +
  stat_function(fun = dnorm,
                args = list(mean = 10, sd = 2.5),
                color = "red",
                linewidth = 2) +
  xlim(c(0,20)) +
  labs(x="X",y="Prob. Density") +
  theme_classic()
  

# curve(dnorm(x,0,2),-10,10)   # base R version is simpler in this case...

# What happens when you try to use a discrete distribution?


# SIMULATE DATA: ------------------
#  decompose into deterministic and stochastic components (linear regression example)

## Deterministic component  -----------------------------------------

#   define function for transforming a predictor variable into an expected response (linear regression)

    # Arguments:
      # x: vector of covariate values (predictor variable)
      # a: the intercept of a linear relationship mapping the covariate to an expected response
      # b: the slope of a linear relationship mapping the covariate to an expected response

deterministic_component <- function(x,a,b) {a + b*x}   # specify a deterministic, linear functional form

xvals = seq(0,100,5)  # define the values of a hypothetical predictor variable (e.g., tree girth)

expected_vals <- deterministic_component(xvals,175,-1.5)   # use the deterministic component to determine the expected response (e.g., tree volume)
expected_vals

plot(xvals,expected_vals,ylab="response mean", xlab="predictor", type="l")   # plot out the relationship

# plot(xvals,expected_vals,type="l")    # alternatively, plot as a line


## Stochastic component -------------------------------------------- 
##    define a function for transforming an expected (deterministic) response and adding a layer of "noise" on top!

    # Arguments:
      # x: vector of expected responses
      # sd: standard deviation of the "noise" component epsilon
stochastic_component <- function(x,sd){ rnorm(length(x),x,sd)}       # add a layer of "noise" on top of the expected response values

    ### Simulate stochastic data!!
sim_vals <- stochastic_component(expected_vals,sd=10)   # try it- run the function to add noise to your expected values. 

plot(xvals,sim_vals)     # plot it- it should look much more "noisy" now!

# ALTERNATIVELY:

sim_vals <- stochastic_component(deterministic_component(xvals,175,-1.5),10)    # stochastic "shell" surrounds a deterministic "core"    


# Goodness-of-fit test! -------------------------------------

    # Do the data fall into the range of plausible values produced by this model?

# Imagine you have the following "real" data (e.g., tree volumes). 

realdata <- data.frame(Volume=c(0.4, 2.7, 5.3, 13.3, 42.4, 63.5, 45.8, 233.4, 213.7, 383.1),Girth=seq(1,10,length=10))
plot(realdata$Girth,realdata$Volume)


# Simulate many datasets from our hypothesized data generating model (intercept=10,slope=4,variance=1000):

lots <- 1000    # specify number to approximate infinity
N <- nrow(realdata)    # define the number of data points we should generate for each simulation "experiment"

simresults = replicate(lots, rnorm(N,10+realdata$Girth*40,31))

    # now make a boxplot of the results
boxplot(t(simresults),xaxt="n",ylab="Volume",xlab="Girth")    # (repeat) make a boxplot of the simulation results
axis(1,at=c(1:N),labels=realdata$Girth)                          # add x axis labels 


# Now overlay the "real" data
    # how well does the model fit the data?

boxplot(t(simresults),xaxt="n",ylab="Volume",xlab="Girth")    # (repeat) make a boxplot of the simulation results
axis(1,at=c(1:N),labels=realdata$Girth)                          # add x axis labels 
points(c(1:N),realdata$Volume,pch=20,cex=3,col="red",xaxt="n")     # this time, overlay the "real" data 


# Let's simulate many datasets from a 'null' model (intercept=100,slope=0,sd=400):

simresults = replicate(lots, rnorm(N,100,400))
boxplot(t(simresults),xaxt="n",ylab="Volume",xlab="Girth")    # (repeat) make a boxplot of the simulation results
axis(1,at=c(1:N),labels=realdata$Girth)                          # add x axis labels 
points(c(1:N),realdata$Volume,pch=20,cex=3,col="red",xaxt="n")     # this time, overlay the "real" data 



# Power analysis example ---------------------------------
#    designing a monitoring program for a rare species


   ### first, let's develop some helper functions:

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
SurvProb(500,people=2,days=7,pObs=0.02)   # test the new function


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

SimulateMonitoringData(N0=1000,trend=-0.03,nyears=25,people=1,days=3,survint=4)    # test the new function


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

IsDecline(monitoringData=c(10,20,NA,15,1),alpha=0.05)    # test the function


## Review lab exercise (lab 2) ----------------------------
##     develop a "power" function to return the statistical power to detect a decline under alternative monitoring schemes...

GetPower <- function(N0,trend,nyears,people,days,survint,alpha){
     # fill this in!
  return(Power)
}

