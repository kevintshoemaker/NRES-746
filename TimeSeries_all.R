library(rnoaa) # package for downloading NOAA, NCDC, and GHCN data
library(forecast) # time series forecasting
library(tseries) # useful functions and tests 
library(tsutils) # more tests
library(tidyverse)
library(lubridate) # date management
library(ggplot2)
library(dataRetrieval) # packaage for downloading USGS stream flow data
library(dplyr) # for data formatting 
library(trend) # for MK test and sen's slope
library(zyp) # y-intercept calc for sen's slope
library(changepoint)
stationID <- 'USW00023185' # Reno airport (1937-03-01 to present)
dataset <- 'GSOM'   # Global Summary of the Month
climateVars <- c('TMIN', 'TMAX', 'PRCP', 'EMXP') # select climate variables to download
startDate <- as.Date('1940-01-01')
endDate <- as.Date('2020-12-31')
GHCN <- getGHCN(dataset, stationID, climateVars, startDate, endDate) # Call to download data using a wrapper function

ggplot(GHCN, aes(as.Date(date), value, color = datatype)) +
  geom_point() + 
  facet_wrap(~datatype, scales = 'free') +
  scale_x_date(date_breaks = '20 years', date_labels = "%Y") +
  labs(x = "", y = "")

tmin <- GHCN %>% 
  filter(datatype == 'TMIN') %>%  # select only 'TMIN' data from downloaded data frame
  select(value) # select only the data value column, we don't need the stationID, date, or datatype columns

# the ts() function converts to an r time series object
# input data is a vector or matrix of data
tmin.ts <- ts(tmin, frequency = 12, start = c(1940, 1))

head(tmin.ts, n = 60) # look at the first 5 years of data
plot.ts(tmin.ts, ylab = "Min. Temp. (C)") # plot the time series

# explore the time series object and other information it contains
start(tmin.ts)
end(tmin.ts)
frequency(tmin.ts)

plot.ts(tmin.ts, ylab = "Min. Temp. (C)") # plot the time series
tsutils::mseastest(tmin.ts, outplot = 1)
tminDecompose <- decompose(tmin.ts, type = "multi")
plot(tminDecompose)

tminSeasonAdj <- tmin.ts - tminDecompose$seasonal

plot.ts(tminSeasonAdj, 
        main = "Seasonally-Adjusted Minimum Monthly Temperature",
        ylab = "Min. Temp. (C)")

tminTrendAdj <- tmin.ts - tminDecompose$trend

plot.ts(tminTrendAdj, 
        main = "Detrended Minimum Monthly Temperature",
        ylab = "Min. Temp. (C)")

tminForecast <- forecast::stlf(tmin.ts, method = "naive")

plot(tminForecast, 
        main = "Minimum Monthly Temperature with Forecasting through 2022",
        ylab = "Min. Temp. (C)")
# identifying parameters for data download
siteNumber <- "11413000" # USGS gauge number
parameterCd <- "00060"  # mean daily discharge in cfs
startDate <- "1940-10-01" # starting date bc missing few years in the 30's
endDate <- "2021-11-01" # date today
# download data using readNWISdv function
yuba_raw <- dataRetrieval::readNWISdv(siteNumber, parameterCd, startDate, endDate)
head(yuba_raw) # check
colnames(yuba_raw)[4] <-"discharge_cfs" # rename column 4
yuba_raw <-yuba_raw[,-5] # delete unused col
class(yuba_raw$Date) # check it's actually a data not a character string
yuba_raw$year <-lubridate::year(yuba_raw$Date) # add year col
yuba_raw$month <-lubridate::month(yuba_raw$Date) # add month col
yuba_raw$day <-lubridate::day(yuba_raw$Date) # add day of month col
yuba_raw$doy <-lubridate::yday(yuba_raw$Date) # add day of year col
yuba_raw <-dataRetrieval::addWaterYear(yuba_raw)
head(yuba_raw) # check
yuba <- as.data.frame(dplyr::group_by(yuba_raw,waterYear) 
                      %>% mutate(dowy = seq(1:n())))
head(yuba)
tail(yuba)
# set theme for plot, changes background to white
ggplot2::theme_set(theme_light(base_size =11))

# plot
ggplot(yuba, aes(y = discharge_cfs, x = Date))+
  geom_line(size = .3) + # change line width
  labs(title="Yuba near Goodyears Bar,CA Discharge", 
       y="Discharge (cfs)", x="Date")
# filter the date down to the 2021 + 1 month water year and plot
yuba_wy2021 <-dplyr::filter(yuba, Date >= as.Date("2020-10-01") & Date <= as.Date("2021-11-01"))

# plot wy 2021
ggplot(yuba_wy2021, aes(y = discharge_cfs, x = Date))+
  geom_line(size = .7, col = "darkblue") + # change width, set color
  labs(title="Yuba near Goodyears Bar,CA Discharge WY 2021", 
       y="Discharge (cfs)", x="Date")
# create a dataframe of maximum daily discharge per year
wy_max <-as.data.frame(yuba %>%
  group_by(waterYear) %>%
  filter(discharge_cfs == max(discharge_cfs, na.rm=TRUE)))

head(wy_max)            

# max day of water year
ggplot(wy_max, aes(y = dowy, x = Date)) +
  geom_point(col = "goldenrod")+
  geom_smooth(method = "lm", se = FALSE)+ # add lm trend line
  labs(title="Yuba Max Annual Discharge DOWY", 
       y="Day of Water Year", x="Date")

# max discharge cfs
ggplot(wy_max, aes(y = discharge_cfs, x = Date)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+ # add lm trend line
  labs(title="Yuba Max Annual Discharge DOWY Discharge", 
       y="Discharge (cfs)", x="Date")
max <-dplyr::filter(yuba, discharge_cfs == max(discharge_cfs, na.rm=TRUE))
print(max)
melt_onset_dowy <-function(df){
  
  # crop days until until spring and early summer  
  wy_spring <-filter(df, dowy > 100 & dowy < 250) # crop days until until spring and early summer
  mean_spring <-mean(wy_spring$discharge_cfs) # find meann discharge of spring date
  above_mean_days <-wy_spring$dowy[which(wy_spring$discharge_cfs > mean_spring)] # find dates which discharge > mean
  diff_vect <-diff(above_mean_days)
  results <-rle(diff_vect)
  meets <-as.numeric(min(which(results$lengths == max(results$lengths))))
  values_to_cut <-results$lengths[1:(meets-1)]
  #sum the lengths
  index_values_until_start <-as.numeric(sum(values_to_cut))
  #subtract that many values off to get your SDD
  onset_dowy <-as.integer(min(above_mean_days[-c(1:index_values_until_start)]))
  return(onset_dowy)
}
# create sequence of available years
years <-as.matrix(seq(min(yuba$waterYear),2021,1))

# create an empty matrix to read the data into
loop_mat <-matrix(nrow = length(years))

# filter data to stop at WY 2021
yuba_loop <-dplyr::filter(yuba, Date >= as.Date("1940-10-01") & Date <= as.Date("2021-09-30"))
tail(yuba_loop)

# loop the function for calculating melt onset through all the way years
for (i in 1:length(years)){
  
  yearly_df <-subset(yuba_loop, waterYear == years[i]) # sequence through years to create yearly df
  loop_mat[i] <-melt_onset_dowy(yearly_df) # apply function and sequentially save in empty matrix

}

# bind years vector and new melt onset dates
melt_onset_df <-as.data.frame(cbind(years, loop_mat))
colnames(melt_onset_df)[1] <-"waterYear"
colnames(melt_onset_df)[2] <-"spring_flow_dowy"
tail(melt_onset_df) # check
# plot the new data and add a "trend line", just visually
# seems to trend towards early flow onset
ggplot(melt_onset_df, aes(y = spring_flow_dowy, x = waterYear)) +
  geom_point(col = "red")+
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title="Yuba Spring Flow Onset Date", 
       y="Day of Water year", x="Water Year")
yuba_wy2021 <-filter(yuba, waterYear == 2021) # subset wy 2021
melt_date <-melt_onset_df$spring_flow_dowy[which(melt_onset_df$waterYear == "2021")]

# for wy2021
ggplot(yuba_wy2021, aes(y = discharge_cfs, x = dowy))+
  geom_line(size = .7, col = "green") + # change width, set color
  geom_vline(xintercept = melt_date, col = "red") + # date calculated
  labs(title="Yuba near Goodyears Bar,CA Discharge WY 2021", 
       y="Discharge (cfs)", x="Date")

# 1970 has the earliest calculated date of 105 or January 15th!
min <-min(melt_onset_df$spring_flow_dowy) 
yuba_wy1970 <-filter(yuba, waterYear == 1970) # low year

# plot to check our function is being reasonable here
# and it is, shows the variability in the hydroclimatic system
ggplot(yuba_wy1970, aes(y = discharge_cfs, x = dowy))+
  geom_line(size = .7, col = "cyan") + # change width, set color
  geom_vline(xintercept = min, col = "red") + # date calculated
  labs(title="Yuba near Goodyears Bar,CA Discharge WY 1970", 
       y="Discharge (cfs)", x="Date")
# test the data for normality
hist(melt_onset_df$spring_flow_dowy, breaks = 30)

# but let's use a Shapiro test to check for real
shapiro.test(melt_onset_df$spring_flow_dowy)  
# We can turned to the Mann-Kendall trend test, which is non-parametric, not requiring a normal dist
# this only tells us if there is a significant trend, not the estimated slope
mk_results <-trend::mk.test(melt_onset_df$spring_flow_dowy)
print(mk_results) # significant

# for a slope estimate we uses the Sen's slope function
sens_result <-trend::sens.slope(melt_onset_df$spring_flow_dowy)
print(sens_result)
slope_est <-as.numeric(sens_result[[1]]) # slope estimate

# the trend package only gives us slope, not a y-int estimate for plotting
# we'll use the zyp package for that
y_int_results <-zyp::zyp.sen(spring_flow_dowy ~ waterYear, melt_onset_df) # y comes first!!
print(y_int_results) # inspect
y_int <-y_int_results$coefficients[[1]] # pull out y-int estimate for ploting
# not let's plot the data with our MK estimate trend, instead of lm one from ggplot
ggplot(melt_onset_df, aes(y = spring_flow_dowy, x = waterYear)) +
  geom_point(col = "red")+
  # geom_smooth(method = "lm", se = FALSE) + 
  # use info from our MK and Sen's test to graph the trend line
  geom_abline(intercept = y_int, slope = slope_est, color="blue", size=1.2) +
  labs(title="Yuba Spring Flow Onset Date MK Trend", 
       y="Day of Water year", x="Water Year")

# quick test plot to visualize the time series
plot(martin$Date, martin$discharge_cfs, type = "l",
     main = "Martin Creek Discharge,NV",
     ylab = "Discharge (cfs)",
     xlab = "Date")
df_day_ts <- ts(data = martin$discharge_cfs, frequency = 365.25, start = c(1921,10), end= c(2021,10))
summary(df_day_ts)
#decompose the data both ways
ddata <- decompose(df_day_ts, type = "additive")
ddatam <- decompose(df_day_ts, type = "multiplicative")
plot(ddata) #additive
plot(ddatam) #multiplicative
#compare which model has smaller sums of squares for the ACF of the residuals
add_info <- sum(acf(ddata$random, lag.max = length(ddata$random), na.action = na.omit)$acf^2)
mult_info <- sum(acf(ddatam$random, lag.max = length(ddatam$random), na.action = na.omit)$acf^2)
add_info
#mult_info
#we will use ddata for the remainder of the demo
#Now we visualize our decomposed data and our trend component
#plot(ddata)
plot(ddata$trend)
abline(reg=lm(ddata$trend~time(ddata$trend)), col="red")
# examine the data for normality
hist(ddata$trend, breaks = 30)

# try autocorrelation function (ACF) (Autoregressive process) and the partial autocorrelation function (PACF) (Moving Average process) to examine the data
acf(df_day_ts,lag.max = length(df_day_ts),
    xlab = "lag #", ylab = 'ACF',main='ACF results in very gradual attenuation')  #possible non-stationarry
pacf(df_day_ts, lag.max = length(df_day_ts), xlab = "lag #", ylab = 'PACF',main='PACF results in rapid dampening') #possible stationary

#Ljung-Box test for independence
#the Ljung-Box test examines whether there is significant evidence for non-zero 
#correlations at given lags. If p-value >0.05, we accept the null hypothesis of independence.
lag.length = 25
options(warn=-1)
Box.test(df_day_ts, lag=lag.length, type="Ljung-Box")
#p-value suggests we should reject the null and assume data is not independent

#Another test we can conduct is the Augmented Dickey–Fuller (ADF) t-statistic 
#test to find if the series has a unit root. Non-stationary data will have a large p-value.
tseries::adf.test(df_day_ts)
#this test suggests that we might have stationary data

#A fourth test we can use to check assumptions is to check for stationarity or trend is the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test.
tseries::kpss.test(df_day_ts, null="Level") #results in p-value = 0.01, reject null
#this test suggests that our data might not be stationary
forecast::ndiffs(df_day_ts)  # number of differences need to make it stationary
forecast::nsdiffs(df_day_ts) # number for seasonal differencing needed
#if we needed to difference (as some of our tests tell us we might):
stationaryCFS <- diff(df_day_ts, differences=1)
#examine differenced data
plot(stationaryCFS, type="l", main="Differenced and Stationary")  # appears to be stationary
hist(stationaryCFS)
lag.length = 25
options(warn=-1)
Box.test(stationaryCFS, lag=lag.length, type="Ljung-Box") #still autocorrelated
tseries::adf.test(stationaryCFS) #stationary
tseries::kpss.test(stationaryCFS) #stationary
#Our data might not meet the normality assumption and it is still autocorrelated. We may need to mean-center our data to reduce problems with forecasting and reduce problems due to temporal autocorrelation
options(warn=-1)
mymodel <- forecast::auto.arima(stationaryCFS, lambda = "auto", parallel = TRUE, num.cores = 8)
mymodel
#let's try forecasting now that we have an ARIMA model built
fit <- forecast::Arima(df_day_ts, c(3,0,1)) #ARIMA model shows three coefficients for the AR component, 0 seasonal differences, and 1 coefficient for the MA component
plot(forecast::forecast(fit))
#we write a function that will help us look for the correct penalty value to set in our change point analysis
cptfn <- function(data, pen) {
  ans <- changepoint::cpt.mean(data, test.stat="Normal", method = "PELT", penalty = "Manual", pen.value = pen) 
  length(cpts(ans)) +1
}
# run our cptfn function expecting there is a signal with a known change point
pen.vals <- seq(0, 12,.2) #set some penalty values to use in the elbow plot
elbowplotData <- unlist(lapply(pen.vals, function(p) 
  cptfn(data = df_day_ts, pen = p))) #apply our function on the data
#and now we plot the penalty parameters as a function of the time series
plot(pen.vals,elbowplotData, 
     xlab = "PELT penalty parameter",
     ylab = " ",
     main = " ")

#look at the graphs, especially the PELT Penalty Parameter to determine
#where the change point might be

penalty.val <- 6 # this value is determined from elbow plots

#next we use cpt.mean and the SegNeigh method with a max number of change points set to 4 to detect change points in the data
options(warn=-1)
cptm_data <- changepoint::cpt.mean(df_day_ts, penalty='Manual',pen.value=penalty.val,method='SegNeigh', Q=4)

cpts_dat <- changepoint::cpts(cptm_data) # change point time points
cpts_dat #gives us the time points

trend_narm <- na.omit(ddata$trend)  #decomposed trend data without the NA's
options(warn=-1)
cptm_trend <- changepoint::cpt.mean(trend_narm, penalty='Manual',pen.value=penalty.val,method='SegNeigh', Q=4) 
cpts_trnd <- changepoint::cpts(cptm_trend) # change point time points
cpts_trnd #gives us the time points

plot(cptm_data,
     xlab = "Time",
     ylab = "Avg. Daily Discharge (cfs)",
     main = "Change in Mean Signal")

plot(cptm_trend,
     xlab = "Time",
     ylab = "Decomposed Trend of CFS",
     main = "Change in Mean Signal")
