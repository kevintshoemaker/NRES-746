
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #1                    ####
####                                                    ####
############################################################


############################################################
####  Model selection and multi-model inference         ####
############################################################



#########
# Load packages

suppressWarnings(library(forecast))
suppressWarnings(library(tseries))
suppressWarnings(library(dtw))


## Set up the time series

year<-seq(1968, 2017, by=1)
winter<-c(33.4, 36.3, 30.3, 33.5, 26.1, 29.5, 29.7, 28.8, 27.6, 27.8, 26.7, 29.0, 27.7,
          24.9, 27.6, 29.0, 22.0, 27.3, 26.9, 23.2, 23.6, 26.7,
          25.8, 21.6, 22.1, 25.8, 21.8, 22.9, 26.9, 20.0, 23.2, 24.7, 21.5, 23.7, 23.9, 21.6, 25.4, 
          24.5, 23.4, 25.1, 26.0, 24.8, 22.2, 25.9, 26.9, 23.7, 24.1, 21.8, 25.4, 23.3)
summer<-c(28.7, 22.8, 28.5, 26.3, 27.8, 22.9, 25.3, 23.7, 25.8,
          28.3, 25.0, 24.9, 22.8, 29.8, 19.7, 17.4, 22.7, 22.1,
          22.6, 26.1, 28.0, 23.0, 23.0, 22.2, 25.2, 19.9, 23.7, 17.7, 21.1,
          19.1, 18.2, 19.2, 19.5, 22.2, 24.7, 21.1, 22.3, 20.4, 17.5, 19.9, 15.4, 18.0, 
          15.8, 15.7, 19.7, 19.4, 23.4, 22.3, 17.2, 16.3)

tahoedatasummer<-data.frame(year, summer)
tahoedatasummer$season<-2
tahoedatawinter<-data.frame(year, winter)
tahoedatawinter$season<-1

names(tahoedatasummer) <- c("year", "clarity","season")
names(tahoedatawinter) <- c("year", "clarity", "season")
tahoedata<-as.data.frame(rbind(tahoedatasummer, tahoedatawinter))
tahoedata <- tahoedata[order(tahoedata$year,tahoedata$season),] 
rownames(tahoedata) <- c()
head(tahoedata)


#Create a time searies using the ts() command 
tahoetimeseries <- ts(tahoedata$clarity, frequency=2, start=c(1968,1))
class(tahoetimeseries) #This confirms that the dataset is now a timeseries object

#Create a time series plot
plot.ts(tahoetimeseries, main="Declining Lake Tahoe Clarity Over Time",
        ylab="Depth of H2O Clarity (meters)", xlab="Year", lwd=1)


# Decompose the time series into trend, seasonal, and random components

tahoetimeseriescomponents <- decompose(tahoetimeseries, type="additive") # specify an additive decomposition 

tahoetimeseriescomponents$trend #portion of variation in H20 Clarity due to the trend
tahoetimeseriescomponents$seasonal #portion of variation in H20 Clarity due to seasonal patterns
tahoetimeseriescomponents$random #portion of variation in H20 Clarity due to random effects
# We can visualize each of the components

plot(tahoetimeseriescomponents)


par(mfrow=c(1,3))
plot.ts(tahoetimeseries, main="Declining Lake Tahoe \nClarity", xlab="Year", ylab="Clarity (meters)")

seasonallyadjusted <- tahoetimeseries-tahoetimeseriescomponents$seasonal #Remove the seasonal component
plot(seasonallyadjusted, main="Seasonal Component \nRemoved", ylab="Clarity (meters)", xlab="Year")

trend<- seasonallyadjusted-tahoetimeseriescomponents$random #With both seasonal and random components removed, only trend remains.
plot(trend, main="Seasonal and Random \nComponents Removed", xlab="Year", ylab="Clarity (meters)")


forecasted <- HoltWinters(tahoetimeseries)
forecasted


plot(forecasted, xlab="Year")


future<-forecast(HoltWinters(tahoetimeseries), h=40)
plot(future)


Box.test(future$residuals, lag=20, type="Ljung-Box")


acf(future$residuals, lag.max=20, na.action = na.pass)


## A function for building correlograms and testing errors

ErrorDistribution<-function(forecastedvalues, lags){
  errors<-na.omit(forecastedvalues$residuals) 
  par(mfrow=c(1,2))
  hist(errors, main="Distribution of Forecast \nErrors", freq=FALSE, col="royalblue1", 
       xlab="Error Magnitude", ylim=c(0,(dnorm(0, 0, sd(errors))+1))) #Histogram of the residual        distribution
  curve(dnorm(x, 0, sd(errors)), type='l', add=TRUE, lwd=3) #Histogram of residuals should follow this curve.
  legend("topleft", lwd=3, c("Ideal Normal Distribution"), cex=0.5) #Create a legend
  acf(errors, main="Corellogram", lag=lags, na.action=na.pass, ylab="Correlation Coefficient") #Plot the corellogram 
  legend("topright", c("95% Confidence \nIntervals\n"), lty=2, col=4, cex=0.6) 
  shapiro<-shapiro.test(errors) 
  ljung.box<-Box.test(errors, lag=lags, type="Ljung-Box") #Ljung-Box test: tests temporal autocorrelation
  matrix<-matrix(data=c(shapiro$statistic, shapiro$p.value, 
                        ljung.box$statistic, ljung.box$p.value), nrow=2, ncol=2)
  colnames(matrix)<-c("Shapiro-Wilk Test", "Ljung-Box Test")
  rownames(matrix)<-c("Test Statistic", "P-Value") #Makes a matrix of the test results.
  return(matrix)
  
}

ErrorDistribution(forecastedvalues=future, lags=20)


# Testing our forecasting ability

plot.ts(window(tahoetimeseries, start=c(1990,1), end=c(2017, 2)), main="Declining Lake Tahoe Clarity Over Time", ylab="Depth of H2O Clarity (meters)", xlab="Year", ylim=c(0, 40), lwd=2)
shortforecast<-data.frame(forecast(HoltWinters(window(tahoetimeseries, start=c(1968, 1), end=c(2009, 2))), h=16)) #Shortened forecast (using 1968-2009 data)
predictionyears<-seq(2010, 2017.5, by=0.5) #Prediction years

#Plot Forecast with 95% prediction intervals 
points(predictionyears, shortforecast$Point.Forecast, type='l', col="red", lwd=2) #Forecast
points(predictionyears, shortforecast$Lo.95, type='l', lty=2, col="red", lwd=1) #lower 95% PI
points(predictionyears, shortforecast$Hi.95, type='l', lty=2, col="red", lwd=1) #Upper 95% PI
legend("topleft", c("Observed Values (2000-2017)","Forecasted Values (2010-2017)", "95% Prediction Intervals"), lwd=c(2,2,1), col=c(1,2,2), lty=c(1,1,2), cex=0.75)


AirPassengers
class(AirPassengers) 
## plot of Air Passengers between 1948 to 1960

plot.ts(AirPassengers)


## Create a Holt-Winters projection of the last 

forecastAirPassengers<-forecast(HoltWinters(AirPassengers, seasonal="mult", gamma=), h=40)
plot.ts(AirPassengers, main="Air Passenger Data", ylab="# of Airline Passengers",xlab="Year", xlim=c(1950, 1960), ylim=c(100, 555), lwd=2)
shortforecast<-data.frame(forecast(HoltWinters(window(AirPassengers, start=c(1950, 1), end=c(1955, 12))),h=48))
predictionyears<-seq(1956, 1959+(11/12), by=1/12)

points(predictionyears, shortforecast$Point.Forecast, type='l', col="red", lwd=2)
points(predictionyears, shortforecast$Lo.95, type='l', lty=2, col="red", lwd=1)
points(predictionyears, shortforecast$Hi.95, type='l', lty=2, col="red", lwd=1)
legend("topleft", c("Observed Values (Jan,1949-Dec,1960)", "Forecasted Values (Jan,1956-Jan,1960)", "95% Confidence Intervals"), lwd=c(2,2,1), col=c(1,2,2), lty=c(1,1,2), cex=0.75)

##Ljung-Box test for autocorrelation of residuals
Box.test(forecastAirPassengers$residuals, lag=20, type="Ljung-Box")


## Correlogram and Ljung-Box test for the first 48 months (= lag of 4 years)

acf(AirPassengers, lag.max=48, na.action = na.pass,xlab="lag (years)")
pacf(AirPassengers, lag.max=48, na.action = na.pass,xlab="lag (years)")
adf.test(AirPassengers) #Augmented Dickey-Fuller Test
kpss.test(AirPassengers) #Kwiatkowski-Phillips-Schmidt-Shin: null hypothesis: x is level or trend stationary


# The diff() command de-trends the timeseries to make it more stationary
# This plot gives: 
# 1. overal shape of data
# 2. ACF (autocorrelation function: represents the correlation between consecutive data points in the time series)
# 3. pACF (partial autocorrelation fuction: partial correlation coefficients between the series and lags of itself) 

dAP<-diff(AirPassengers, differences = 1)
tsdisplay(ts(dAP, freq=1),lag.max= 48, main="First order Difference", xlab="Months")
adf.test(dAP)
kpss.test(dAP)


ndiffs(AirPassengers)


# Use the diff() command on the logged data and create ACF and pACF plots

ldAP<-diff(log(AirPassengers), differences = 1) #use diff function to take first order diffrence of logged time series
tsdisplay(ts(ldAP, freq=1),lag.max= 48, main="First order Difference",xlab="Months")
adf.test(ldAP)
kpss.test(ldAP)


# Compare differenced and logged time series and correlograms

par(mfrow=c(2,3),mar=c(4, 4, 1, 1) + 0.1)
acf(AirPassengers, lag.max=48, na.action = na.pass, main="Original Time Series",xlab="Lag (years)")
acf(dAP, lag.max=48, na.action = na.pass, main="Differenced Time Series",xlab="Lag (years)")
acf(ldAP, lag.max=48, na.action = na.pass, main="Differenced and Logged Time Series",xlab="Lag (years)")

plot.ts(AirPassengers, xlab="Year")
plot.ts(dAP, xlab="Year", ylab="logged # of Air Passengers")
plot.ts(ldAP, xlab="Year", ylab="Differenced and logged # of Air Passengers")


sldAP<-diff(ldAP,lag=12, differences= 1) # Take first order difference of seasonal cycle (12mo)
suppressWarnings(tsdisplay(ts(sldAP,freq=1),main="Logged Raw Data \nDifferenced at lags 1 and 12",lag.max=40))
par(mfrow=c(1,2))
acf(sldAP, lag.max=48, na.action = na.pass, main="Differenced Time Series",xlab="Lag (years)")
plot.ts(sldAP, xlab="Year", ylab="log ( # of Air Passengers )")
adf.test(sldAP)
kpss.test(sldAP)


arima_model <- auto.arima(log(AirPassengers))
summary(arima_model)

forecastAP <- forecast(arima_model, level = c(95), h = 40)
plot(forecastAP)


## Test for normality and autocorrelation of residuals

ErrorDistribution(forecastedvalues=forecastAP, lags=20)


dat<-seq(0,6.28,len=100)
query<-sin(dat)+runif(100)/10 #create somewhat noisy data

## Create a template or baseline time series to compare our query time series to
template<-cos(dat)
plot(template); lines(query,col="blue")


## Find the best match with the canonical recursion formula
alignment<-dtw(query,template,keep=TRUE) #will need the dtw package for this

plot(dtw(query,template,keep=TRUE, step=rabinerJuangStepPattern(6,"c")),type="twoway")
plot(dtw(query,template,keep=TRUE, step=rabinerJuangStepPattern(6,"c")),type="threeway")


#View ONI data
onidat=read.csv("onidata.csv", header=T)
head(onidat)
onits <- ts(onidat$ONI, frequency=2, start=c(1968,1)) #convert to time series
plot.ts(onits, main="ONI Over Time",
        ylab="Oceanic Niño Index (ONI)", xlab="Year", lwd=1)

