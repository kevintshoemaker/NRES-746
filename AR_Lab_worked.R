#~~~~~~~~~~~~~~~~~~~
#Code to accompany the Timeseries & AR model student led lab
#NRES 746 Fall 2023
#Top Level----
#~~~~~~~~~~~~~~~~~~~


#Exercise 1---

 # GOAL: explore ACF and PACF plots

library("MARSS")

# load the data from the "MARSS" package
data(lakeWAplankton, package = "MARSS")

phyto <- lakeWAplanktonTrans #Give it a name to keep the original data structure in tact
phyto <- phyto[phyto[,"Year"] > 1975,] #Subset to years greater than 1975 (missing data in early years)

# Note that our data is monthly!
head(phyto)

#Turn into time series object
phyto.ts <- ts(phyto, start = c(1976,1), freq = 12)

#Plot the ACF
acf(phyto.ts[,"Cryptomonas"], na.action = na.pass, las = 1)
#Plot the PCF
pacf(phyto.ts[,"Cryptomonas"], na.action = na.pass, las = 1)

plot(phyto.ts[,"Cryptomonas"])

#Exercise 2---

 # GOAL: basic frequentist autoregressive model
     # (try to incorporate a seasonal covariate.. )

head(phyto)

phyto.lag <- as.data.frame(phyto)
phyto.lag$lag1 <- NA
phyto.lag$lag1[2:nrow(phyto.lag)] <- phyto.lag$Cryptomonas[1:(nrow(phyto.lag)-1)]
library(lubridate)
phyto.lag$DATE <- ymd(paste(phyto.lag$Year,phyto.lag$Month,"1",sep="_"))

tokeep <- c("Cryptomonas","lag1","DATE","Month")
  
names(phyto.lag)
phyto.lag <- phyto.lag[,tokeep]
head(phyto.lag)
phyto.lag <- na.omit(phyto.lag)

fit.lm <- lm(Cryptomonas~lag1 + sin(2*pi*Month/12),data=phyto.lag)

summary(fit.lm)

nd <- data.frame(Month=1:12,lag1=0)
plot(nd$Month, predict(fit.lm,nd),type="l")

plot(Cryptomonas~DATE,data=phyto.lag,type="l")

lines(phyto.lag$DATE,predict(fit.lm),col="green",lwd=2)


# access the fitted series (for plotting)
fit <- fitted(fit.lm)  

# find predictions for original time series
pred <- predict(fit.lm, newdata=data.frame(Time=Time))    

plot(temperature ~ Time, data= Data, xlim=c(1, 900))
lines(fit, col="red")
lines(Time, pred, col="blue")


 # did

#Exercise 3----

##Building Jags Model----

###Data Management and Exploration----

library(jagsUI)

data <-  read.csv('portal_timeseries.csv') #Read in the data, will need to edit to wherever you have it stored

head(data)

plot.ts(data$NDVI) #Plot the time series

#Plot AC plots, anything interesting? :)
acf(data$NDVI)
pacf(data$NDVI) 

#For secret (forecasting) reasons we're going to cut the last 10 observations out
oos <- data[-(1:(nrow(data)-10)),] #Take the last 100 rows of the data frame and save it as a new dataframe
data.training <- data[(1:(nrow(data)-10)),] #Take out the last 100 rows and the rest of the data is what we will use to make our model



###Build our jags model!----

filename <-  "NDVI_Model.txt"

cat("

model{

#Likelihood 

for(t in 2:time){ #Loop through the number of observations, we have to start at t=2 because at t=1 there is no t-1 data to pull from!
  NDVI[t] ~ dnorm(expNDVI[t],precNDVI)
  expNDVI[t] <- b0 + b1*NDVI[t-1] + b2*rain[t-1]
}

#Priors
b0 ~ dnorm(0,0.1) #intercept term - dnorm uses mu and precision (1/var). So prec of 0.1 is a var of 10
b1 ~ dnorm(0,0.1) #slope term for previous month's NDVI
b2 ~ dnorm(0,0.1) #slope term for previous month's RAIN

sd ~ dnorm(0,0.1)T(0,) #We're truncating the prior to be positive

#Derived terms
precNDVI <- pow(sd,-2) #compute precision from stdv. 1/var

}
",
file = filename)

###Package Data for Jags----

datalist <- list(time = nrow(data.training), NDVI = data.training$NDVI, rain = data.training$rain)

#Specify initial values
#   Each chain needs it's own initial value for the parameters 
#   So we'll create a function to generate initial values for us

inits <- function(){
  vals <- list(
    b0 = runif(1,-10,10),
    b1 = runif(1,-10,10), 
    b2 = runif(1,-10,10),
    sd = runif(1,1,10)
  )
  return(vals)
}

inits() #testing
datalist

#Specify what parameters we're interested in (what do we obtain a joint posterior dist for)
#   What are our free params
params <- c("b0", "b1", "b2","sd")

###Run jags----

nc <- 3 #number of chains
niter <- 10000 #Number of iterations, includes burnin
nb <- 2500 #burnin
thin <- 5 #thining 

#Inits should be function itself, not just one run of it
#N.adapt is the tuning parameter 
mod <- jags(data = datalist, inits = inits, parameters.to.save = params, 
            model.file = filename,
            n.chains = nc, n.adapt=1000, n.iter = niter, n.burnin=nb, n.thin=thin,
            parallel=F)

###Visualize posteriors and assess convergence----
mod #Summary of the model
plot(mod) #Plot trace plots and posterior distributions of our parameters

mod
sims <-  mod$sims.list #Pull out each iteration parameter estimates


##Forecasting----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#If you feel up to it, let's forecast
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#We're going to try and see how well we can forecast those last 10 samples we kept out!

#Step 1: Rebuild the likelikehood

Predictions <- matrix(NA, length(sims$b0), nrow(oos)) #Create a blank matrix to store our predictions
dim(Predictions) # nrow = # of iterations, ncol = years we're making predictions for (withheld data)

#We have to keep parameter estimates together for each iteration so we need a double for loop

for(p in 1:length(sims$b0)){ #Loops through the number of iterations
  NDVI <- data$NDVI[1] #give it a value to start at
  for(t in 1:nrow(oos)){
    NDVI <- sims$b0[p]+sims$b1[p]*NDVI+sims$b2[p]*oos$rain[t] #Copy the process from what you did in the model
    Predictions[p,t] <- NDVI
  }
}

#Take some summary statistics from our predictions to plot
Mean <- apply(Predictions,2, mean)
UCI <- apply(Predictions,2,quantile, prob = .975)
LCI <- apply(Predictions,2,quantile, prob = .025)

#Now let's plot our forecast
par(mfrow = c(1,1))
plot(Mean, type = 'l', ylim = c(0.1,0.4))
lines(UCI, lty = 2, col = 'steelblue')
lines(LCI, lty = 2, col = 'steelblue')
points(oos$NDVI, col = 'red') #Actual data

#Then you can do things like RMSE to get a number on how far our predictions are to the observed data







