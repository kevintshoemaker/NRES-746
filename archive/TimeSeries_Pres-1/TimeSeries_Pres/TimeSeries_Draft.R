## ------------------------------------------------------------------------
#install.packages(c("forecast","tseries", "lubridate", "tidyverse", "GeneCycle"))
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)


## ------------------------------------------------------------------------
#import data
df <- read_csv("CENMET.csv")
str(df)


## ------------------------------------------------------------------------
df$date <- ymd(df$date)
str(df)


## ------------------------------------------------------------------------
ggplot(df, aes(date, tmax)) + geom_line() + scale_x_date('month') +
  ylab("Tmax") + 
  xlab("") +
  theme_bw()
  


## ------------------------------------------------------------------------
df1 <- df[c(1:1100),]
ggplot(df1, aes(date, tmax)) + geom_line() + scale_x_date('month') +
  ylab("Tmax [K]") + 
  xlab("") +
  theme_bw()


## ------------------------------------------------------------------------
#add month again for TS analysis
df1$month<- month(df1$date)
df1$year <- year(df1$date)


## ------------------------------------------------------------------------
#plot month over month DO to see the range and outliers -- We see some variation within months between years
ggplot(df1, aes(date, tmax)) + geom_point(color = "navyblue") +
  facet_wrap(~month) + scale_x_date('month') + 
  ylab("Tmax [K]") +
  xlab("Wrap by month for all years")+
  theme_bw()



## ------------------------------------------------------------------------
#Create a time series object based on Tmax to pass to tsclean()
TS.Tmax <- ts(df1[,c("tmax")])


## ------------------------------------------------------------------------
#tsclean() function to ID and replace outliers and input missing values if they exist
df1$clean_count <- tsclean(TS.Tmax)

#Graph cleaned data - ignore type as scale warning
ggplot() + 
  geom_line(data = df1, aes(x = date, y = clean_count)) +
  ylab('Clean Count') +
  theme_bw()


## ------------------------------------------------------------------------
#Get monthly and weekly moving averages (MA) and compare to cleaned daily data which
#still has a lot of variance and volatility in it.
df1$plotdata.ma7 <- ma(df1$clean_count, order=7)
df1$plotdata.ma30 <- ma(df1$clean_count, order=30)


## ----fig.height = 9, fig.width=12, fig.align="center"--------------------
ts1 <- ggplot() +
  geom_line(data = df1, aes(x = date, y = clean_count, colour = "Counts")) +
  geom_line(data = df1, aes(x = date, y = plotdata.ma7, colour = "Weekly Moving Average")) +
  geom_line(data = df1, aes(x = date, y = plotdata.ma30, colour = "Monthly Moving Average")) + 
  ylab('Tmax [K]') +
  xlab('Date') + 
  theme_bw()

ts1


## ------------------------------------------------------------------------
count_ma = ts(na.omit(df1$plotdata.ma7), frequency = 52)
decomp = stl(count_ma, s.window = "periodic")
deseasonal_cnt <- seasadj(decomp) #This function removes the seasonal component of the data
plot(decomp)


## ------------------------------------------------------------------------
df2 <- decompose(count_ma, type = "additive")
plot(df2)


## ------------------------------------------------------------------------
deseasonal_df <- count_ma - df2$seasonal
detrend_df <- count_ma - df2$trend
derandomize_df <- count_ma - df2$random
par(mfrow=c(1,3))
plot(count_ma, main = "Temp time series")
plot(detrend_df, main = "Deseasonalized")
plot(derandomize_df, main = "Random component removed")


## ------------------------------------------------------------------------
df3 <- deseasonal_df - df2$random
plot(df3)


## ------------------------------------------------------------------------
adf.test(count_ma, alternative = "stationary") 


## ------------------------------------------------------------------------
#ACF plots display correlation between a series and its lags
Acf(count_ma, main = '')

#PACF plots display correlation between a series and its lags that explained by previous lags
Pacf(count_ma, main = '')


## ------------------------------------------------------------------------
count_d1 = diff(deseasonal_cnt, differences = 1) #difference of 1 is sufficient
plot(count_d1)


## ------------------------------------------------------------------------
adf.test(count_d1, alternative = "stationary")
#with difference of 1 data, we bring the 
#Dickey-Fuller test to -10.201, 
#Lag order to 10 and 
#p value to 0.01 (or less) 


## ------------------------------------------------------------------------
#Look for spikes at specific lag points of the differenced series
par(mfrow=c(2,1))
Acf(count_d1, main='ACF for differenced Series')
Pacf(count_d1, main = 'PACF for differenced Series')


## ------------------------------------------------------------------------
#Evaluate and iterate - does the model make sense?
fit <- auto.arima(deseasonal_cnt, seasonal = FALSE)
tsdisplay(residuals(fit), lag.max = 60, main = "(1,1,1) Model Residuals") #high lag max to make sure we're seeing all of the lags


## ------------------------------------------------------------------------
#graph shows serious lags at 7, so modify for q = 7
fit2 = arima(deseasonal_cnt, order=c(1,1,7))
tsdisplay(residuals(fit2), lag.max = 20, main = "Model Residuals")
#Not perfect yet, but getting better


## ------------------------------------------------------------------------
#Test model performance with a holdout set
hold <- window(ts(deseasonal_cnt), start = 991)
fit_no_holdout = arima(ts(deseasonal_cnt[-c(991:1090)]), order = c(1,1,7))
fcast_no_holdout <- forecast(fit_no_holdout,h=100)
plot(fcast_no_holdout, main='ARIMA Forecast')
lines(ts(deseasonal_cnt))


## ------------------------------------------------------------------------
fit3 = arima(deseasonal_cnt, order=c(1,1,12))
tsdisplay(residuals(fit3), lag.max = 20, main = "Model Residuals")
fcast1 <- forecast(fit3, h=60)
plot(fcast1)


## ------------------------------------------------------------------------
hold <- window(ts(deseasonal_cnt), start = 991)
fit_no_holdout = arima(ts(deseasonal_cnt[-c(991:1090)]), order = c(1,1,12))
fcast_no_holdout <- forecast(fit_no_holdout,h=100)
plot(fcast_no_holdout, main='ARIMA Forecast')
lines(ts(deseasonal_cnt)) #a bit better, but not ideal


## ------------------------------------------------------------------------
#get auto fit p,d,q values
auto.arima(deseasonal_cnt, seasonal = FALSE) #ideal ARIMA is (p=1,d=0,q=4)
fit4 = arima(deseasonal_cnt, order=c(1,0,4))
tsdisplay(residuals(fit4), lag.max = 20, main = "Model Residuals")
fcast2 <- forecast(fit4, h=60)
plot(fcast2)


## ------------------------------------------------------------------------
hold <- window(ts(deseasonal_cnt), start = 991)
fit_no_holdout = arima(ts(deseasonal_cnt[-c(991:1090)]), order = c(1,0,4))
fcast_no_holdout <- forecast(fit_no_holdout,h=100)
plot(fcast_no_holdout, main='ARIMA Forecast')
lines(ts(deseasonal_cnt)) 


## ------------------------------------------------------------------------
xs <- seq(0,4*pi,pi/100)
wave.1 <- sin(xs)
plot(xs,wave.1,type="l"); abline(h=0,lty=3)



## ------------------------------------------------------------------------

#Simple function to plot a wave
wave.func <- function(amp = 1, freq = 4, phase = 0, trans = 0, length = 30) {
  xs <- seq(0,length,1/100)
  wave <-  trans + amp*sin(freq*xs+phase)
  plot(xs,wave,type="l")
  abline(h=0,lty=3)
  return(wave)
}


## ------------------------------------------------------------------------
x =  wave.func(amp = 3, freq = 4, phase = 1, trans = 1)

## ------------------------------------------------------------------------
y = wave.func(amp = 3, freq = .5, phase = 1, trans = 1)

## ------------------------------------------------------------------------
z = wave.func(amp = 5, freq = 1, phase = 3, trans = 2)


## ------------------------------------------------------------------------
xs <- seq(0,30,1/100)
wave.total = x+y+z
plot(xs,wave.total,type="l")


## ------------------------------------------------------------------------
plot.fourier <- function(fourier.series, f.0, ts) {
  w <- 2*pi*f.0
  trajectory <- sapply(ts, function(t) fourier.series(t,w))
  plot(ts, trajectory, type="l", xlab="time", ylab="f(t)"); abline(h=0,lty=3)
}


## ------------------------------------------------------------------------
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}



## ------------------------------------------------------------------------

acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time

dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components

f   <- function(t,w) { 
  dc.component + 
    sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}

plot.fourier(f,f.0,ts=ts)



## ------------------------------------------------------------------------
w <- 2*pi*f.0
trajectory <- sapply(ts, function(t) f(t,w))
head(trajectory,n=30)

X.k <- fft(trajectory)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits=c(0,20))


## ------------------------------------------------------------------------
set.seed(101)
acq.freq <- 200
time     <- 1
w        <- 2*pi/time
ts       <- seq(0,time,1/acq.freq)
trajectory <- 3*rnorm(101) + 3*sin(3*w*ts)
plot(trajectory, type="l")

X.k <- fft(trajectory)
plot.frequency.spectrum(X.k,xlimits=c(0,acq.freq/2))



## ------------------------------------------------------------------------
trajectory1 <- trajectory + 25*ts # let's create a linear trend 
plot(trajectory1, type="l")


## ------------------------------------------------------------------------
X.k <- fft(trajectory1)
plot.frequency.spectrum(X.k,xlimits=c(0,acq.freq/2))



## ------------------------------------------------------------------------
trend <- lm(trajectory1 ~ts)
detrended.trajectory <- trend$residuals
plot(detrended.trajectory, type="l")


## ------------------------------------------------------------------------
X.k <- fft(detrended.trajectory)
plot.frequency.spectrum(X.k,xlimits=c(0,acq.freq/2))



## ------------------------------------------------------------------------
library(zoo) #this is used to turn date in to an index


weather <- read.csv("CENMET.csv")       # daily conditions (1 Hz = 1 day)
weather <- weather[order(nrow(weather):1),]

plot(index(weather$date), weather$tmax, type="l")
trend <- lm(tmax ~ index(date), data = weather)
abline(trend, col="red")


## ------------------------------------------------------------------------
#Although this does not have a trend, lets detrend it anyway
detrended.trajectory <- trend$residuals
plot(detrended.trajectory, type="l", main="detrended time series")

#There is also a package called GeneCycle which will perform the fft
library(GeneCycle)

f.data <- GeneCycle::periodogram(detrended.trajectory)
harmonics <- 1:50 
plot(f.data$freq[harmonics]*length(detrended.trajectory), 
     f.data$spec[harmonics]/sum(f.data$spec), 
     xlab="Harmonics (Hz)", ylab="Amplitute Density", type="h")

