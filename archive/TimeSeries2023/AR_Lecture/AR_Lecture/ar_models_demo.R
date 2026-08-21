#~~~~~~~~~~~~~~~~~~~
#Code to accompany the Timeseries & AR model student led lectures
#NRES 746 Fall 2023
#Top Level----
#~~~~~~~~~~~~~~~~~~~

set.seed(123) #Set seed to get reproduceable results

#install.packages("MARSS") #if you don't have it!

###
#Examples of time series----
###

##Nile Example----

?ts

#Load the Nile dataset
data(Nile) 

class(Nile) #ts object already!

#It's easy to turn a ts object back into a regular dataframe
Nile <- data.frame(Y=as.matrix(Nile), date=time(Nile))

#Look at the data
head(Nile) 

#Looks like our response variable is yearly!

#Now let's turn it back into a timeseires... for the example

#Turn it into a timeseries
Nile.ts <- ts(data = Nile$Y, start = 1871, end = 1970, frequency = 1) #our data will be column 'Y', the discharge!And frequency is 1 since it is measured once a year! 

ts.plot(Nile.ts) #ts.plot() plots time series objects

#Plot our original dataframe that wasn't a timeseries object
ts.plot(Nile$Y) #remember we plot our response variable, Y

##Beaver Example----

#Load the beavers dataset from R
data(beavers) #We'll work with beaver1 :)

#Look at the data
head(beaver1)

summary(beaver1) 

#Subset out only the first day
beaver1 <- beaver1[which(beaver1$day == 346),] #Take all the data on day 346

#We'll be looking at beaver temperatures over time!

#Turn it into a timeseries
beaver1.ts <- ts(data = beaver1$temp) #our data will be the temperature. 

ts.plot(beaver1.ts) #ts.plot plots time series objects

###
#Correlograms (ACF and PACF Plots)----
###

#Example ACF plots (simulated data)----

##Building blank ACF plots----

#Build out the blank ACF plot with no values for ACF
par(mai=c(1,1,0,0), omi=c(0.1,0.1,0.1,0.1))
plot(NA, NA, type="n", xlim=c(0,15), ylim=c(-1,1),
     xlab="", xaxt="n", ylab="", las = 1)
abline(h=0)
axis(side=1, at=seq(15), labels=FALSE)
axis(side=1, at=seq(0,15,5))
mtext(expression(paste("Lag ", (italic(k)))), side=1, line=3, cex=1.2)
mtext(expression(paste("ACF ", (italic(r[k])))), side=2, line=3, cex=1.2)

#Build out the ACF plot at lag = 0 (where ACF always = 1! Always 100% correlated to itself!)
par(mai=c(1,1,0,0), omi=c(0.1,0.1,0.1,0.1))
plot(NA, NA, type="n", xlim=c(0,15), ylim=c(-1,1),
     xlab="", xaxt="n", ylab="", las = 1)
abline(h=0)
axis(side=1, at=seq(15), labels=FALSE)
axis(side=1, at=seq(0,15,5))
mtext(expression(paste("Lag ", (italic(k)))), side=1, line=3, cex=1.2)
mtext(expression(paste("ACF ", (italic(r[k])))), side=2, line=3, cex=1.2)

lines(c(0,0), c(0,1), lwd=2, col="darkred")
text(x=1, y =1, expression(italic(r)[0] == 1), col="darkred")

#Add approximate 95% confidence intervals
par(mai=c(1,1,0,0), omi=c(0.1,0.1,0.1,0.1))
plot(NA, NA, type="n", xlim=c(0,15), ylim=c(-1,1),
     xlab="", xaxt="n", ylab="", las = 1)
abline(h=0)
axis(side=1, at=seq(15), labels=FALSE)
axis(side=1, at=seq(0,15,5))
mtext(expression(paste("Lag ", (italic(k)))), side=1, line=3, cex=1.2)
mtext(expression(paste("ACF ", (italic(r[k])))), side=2, line=3, cex=1.2)

lines(c(0,0), c(0,1), lwd=2, col="darkred")
text(x=1, y =1, expression(italic(r)[0] == 1), col="darkred")

# add 95% CI's
nn <- 30
alpha <- 0.05
ts.SD <- qnorm(1-alpha/2, 0, 1)/sqrt(nn)
abline(h=-ts.SD, lty="dashed", col="blue")
abline(h=ts.SD, lty="dashed", col="blue")

##Trend Only ACF Plot----
## length of ts
nn <- 100

## trend only
par(mfrow=c(1,2), mai=c(1,1,0,0), omi=c(0.1,0.1,0.6,0.1))
tt <- seq(nn)
plot.ts(tt, ylab=expression(italic(x[t])))
acf(tt)
mtext("Linear trend {1,2,3,...,100}", outer=TRUE, line=1, cex=1.5)

##Seasonal ACF Plot----
par(mfrow=c(1,2), mai=c(1,1,0,0), omi=c(0.1,0.1,0.6,0.1))
## compute the 2 predictor variables
tt <- sin(2*pi*seq(nn)/12)
plot.ts(tt, ylab=expression(italic(x[t])))
acf(tt)
mtext("Discrete (monthly) sine wave", outer=TRUE, line=1, cex=1.5)

##Trend & Seasonal ACF Plot----
par(mfrow=c(1,2), mai=c(1,1,0,0), omi=c(0.1,0.1,0.6,0.1))
## compute the 2 predictor variables
tt <- sin(2*pi*seq(nn)/12) - seq(nn)/50
plot.ts(tt, ylab=expression(italic(x[t])))
acf(tt, lag.max=30)
mtext("Linear trend + seasonal effect", outer=TRUE, line=1, cex=1.5)

#Beaver ACF & PACF Example ----
##ACF Plots (Beaver data)----

ts.plot(beaver1$temp) #Plot the time series of beaver temperature again

acf(beaver1$temp)

##PACF Plots (Beaver data)----
pacf(beaver1$temp)

#Phytoplankton ACF and PACF Example----
library(MARSS) #load up the MARSS package which is home to the phytoplankton data

data(lakeWAplankton) #read in the data

#Data Management
phyto <- lakeWAplanktonTrans #Give it a name to keep the original data structure in tact
phyto <- phyto[phyto[,"Year"] > 1975,] #Subset to years greater than 1975 (missing data in early years)

head(phyto) #Our data is in monthly increments now!

#Turn into time series object
phyto.ts <- ts(phyto, start = c(1975,1), freq = 12)

#Plot the time series again
par(mai=c(1,1,0,0), omi=c(0.1,0.1,0.1,0.1)) #changes margins
plot.ts(phyto.ts[,"Cryptomonas"], ylab=expression(log(italic(Cryptomonus))), las = 1)

#Plot the ACF
acf(phyto.ts[,"Cryptomonas"], na.action = na.pass, las = 1)

#Plot the PCF
pacf(phyto.ts[,"Cryptomonas"], na.action = na.pass, las = 1)

###
#Simple Time Series Models----
###

##White Noise time series model and ACF----
par(mfrow = c(1,2), mai = c(1.5,0.9,0.1,0.1), omi = c(0,0,0,0))
tt <- rnorm(100)
plot.ts(tt, ylab = expression(italic(w[t])))
acf(tt)

##Random Walk time series model and ACF----
par(mfrow = c(1,2), mai = c(1.5,0.9,0.1,0.1), omi = c(0,0,0,0))
tt <- cumsum(rnorm(100))
plot.ts(tt, ylab = expression(italic(x[t])))
acf(tt)

##Autoregresive Model w/ Phytoplankton----

#Replot the time series
plot.ts(phyto.ts[,"Cryptomonas"], ylab=expression(log(italic(Cryptomonus))), las = 1)


par(mfrow = c(1,2))
#Plot the ACF
acf(phyto[,"Cryptomonas"], na.action = na.pass, las = 1)

#Plot the PACF
pacf(phyto[,"Cryptomonas"], na.action = na.pass, las = 1)


colnames(phyto)


class(phyto) #not a ts object and that is ok!! 

#Let's data manage to make things more readable

#Pick out the columns we're interested in and turn it into a dataframe so we can call them nicely
phyto.df <- as.data.frame(phyto[,c("Year", "Month", "Temp", "TP", "Cryptomonas")])


#Sometimes working with previous times requires clever indexing here let's test out how we can do that 


#proof of concept that indexing is ok on our full dataset
test <- cbind(phyto.df $Cryptomonas[-1], phyto.df $Cryptomonas[-nrow(phyto)]) #offsets the data so t and t-1 line up
colnames(test) <- c("t", "t minus 1")

head(test)

#Now lets test some models!!

#No covariates model
ar1 <- lm(phyto.df $Cryptomonas[-1]~phyto.df $Cryptomonas[-nrow(phyto.df )])
summary(ar1)

#Covariates model
ar1.co <- lm(phyto.df $Cryptomonas[-1]~phyto.df $Cryptomonas[-nrow(phyto.df )]+phyto.df $Temp[-1]+phyto.df $TP[-1])
summary(ar1.co)

#Covariates model with phosphorus at t-1
ar1.co.prevp <- lm(phyto.df $Cryptomonas[-1]~phyto.df $Cryptomonas[-nrow(phyto)]+phyto.df $Temp[-1]+phyto.df $TP[-nrow(phyto.df )])
summary(ar1.co.prevp)

##Adding previous lag steps

#Find the matches of cryptomonas at t and t-1
test1 <- cbind(phyto.df[-1,], phyto.df[-nrow(phyto), c("Cryptomonas")]) 

colnames(test1) <- c(colnames(test1)[1:5], "t min 1")

#Find the matches of cryptomonas at t and t-5
test2 <- cbind(phyto.df[-(1:5),],
               phyto.df[-((nrow(phyto)-4):nrow(phyto)), c("Cryptomonas")]
)
colnames(test2) <- c(colnames(test2)[1:5], "t min 5")

#Merge together
complete <- merge(test1, test2)

complete <- complete[order(complete$Year, complete$Month),] #Reorder our columns

head(complete) #Make sure it's correct

#Run the regression
ar1.co.season <- lm(complete$Cryptomonas~complete$`t min 1`+ complete$`t min 5`)
summary(ar1.co.season)
