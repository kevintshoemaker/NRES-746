
###############################################################
####                                                       ####  
####  Quantile Regression                                  ####
####  Cheyenne Acevado, Christine Albano, Ashley Eustis    #### 
####  NRES 746                                             ####
####                                                       #### 
###############################################################



# install.packages("quantreg")
# install.packages("Qtools")
# install.packages("ggplot2")
library(quantreg)
library(Qtools)
library(ggplot2)


# Simulate data
# generate some data to represent a dataset with a non-constant variance
x <- rep(1:195, 2) # our independent variable
b_0 <- 0 # intercept
b_1 <- 1 # slope
sigma <- 0.1 + 0.5*x # our non-constant variance
eps <- rnorm(x, mean = , sd = sigma)
y <- b_0 + b_1*x + eps


#plot the data to show heteroscedascity
plot(y~x)



# OLS regression
# The standard OLS (Ordinary Least Squares) model explains the relationship between independent variables and the conditional mean of the dependent variable.

# Run OLS regression and get the output:
OLSreg <- lm(y~x)
summary(OLSreg)

# Now we can build some quantile regression models for different quantiles

# estimate the model at the median (this is also known as the 2nd quantile)
qr1 <- rq(y ~ x) # our default tau is 0.50 
summary(qr1) # get the output

# Now lets fit a model for the remaining quantiles

# estimate the model for the first quartile
qr2 <- rq(y ~ x, tau = 0.25) # tau = 0.25 is the first quantile
summary(qr2) # get the output

# estimate the model for the third quartile
qr3 <- rq( y~ x, tau = 0.75) # third quantile
summary(qr3) # get the output

# estimate the model for the 95th quantile
Qreg95 <- rq(y~x, tau=0.95)
summary(Qreg95, se = "rank") # get the output



# plot the data
plot(x,y, type = "n")
points(x,y,cex=.5,col="blue")  # add the data points
taus <- c(.05,.1,.25,.75,.9,.95) # choose the quantiles you want to plot
f <- coef(rq((y)~(x),tau=taus)) # make your coefficents from your quantile model
xx <- seq(min(x),max(x),190) # sequence the min and max values
yy <- cbind(1,xx)%*%f 
for(i in 1:length(taus)){  # loop through your quantile values to make a line for each quantile
  lines(xx,yy[,i],col = "gray") # this will show the slope of th eline for each quantile
}
abline(lm(y ~ x),col="red",lty = 2) # plot the mean Least Squares Estimate fit line
abline(rq(y ~ x), col="blue") # plot the meadian
legend("topleft",legend = c("mean (LSE) fit", "median (LAE) fit"), col = c("red","blue"),lty = c(2,1)) # make the legend!


#Example of plotting of coefficients and their confidence bands
plot(summary(rq(y~x,tau = 1:49/50)))


# use anova to compare the 1st and 3rd quantiles
anova(qr2, qr3)

# compare the 1st, 2nd, and 3rd quantiles
anova(qr1, qr2, qr3)


# Let's use the Barro Data to explore quantile regression models with a larger dataset with more covarites. 
# This is a regression data set consisting of 161 observations on determinants of cross country GDP growth rates with 14 covariates.

data(barro) # load the Barro Data
# ?barro # for more info on this data

names(barro)
# our dependent variable is "y.net" the Annal Change Per Capita GDP
# our independent variables are everything else

# Now let's build a nested model where we compare the effect of Initial Per Capita GDP (lgdp2), Female Secondary Education (fse2), and Education/GDP (gedy2) with our dependent variable (Annual Change Per Capita GDP)

# estimate this model for the 2nd quantile (or the meadian) 
fit0 <- rq(y.net ~  lgdp2 + fse2 + gedy2 , data = barro)  
summary(fit0) # the default is the meadian since we didn't specify tau

# Let's build another nested model where we compare the effect of Initial Per Capita GDP (lgdp2), Female Secondary Education (fse2), Education/GDP (gedy2), Investment GDP (Iy2), and Public Consumption/GDP (gcony2) with our dependent variable (Annual Change Per Capita GDP)

fit1 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro) 
summary(fit1) # this model will also output the 2nd quantile

# Now compare the two models
anova(fit1,fit0) 


# Still using the barro data
# Let's examine different quantiles from our second nested model above


# estimate the model for the third quantile
fit2 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.75)
summary(fit2) # this model gives us the output for the 3rd quantile because we specified tau as 0.75

# estimate the model for the first quantile
fit3 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.25) 
summary(fit3) # this model outpus the 1st quantile

# Test whether the coefficients are significantly different for the quantiles using ANOVA
# remember: use the anova function after you have computed test statistics for two or more quantile regression fits

# Now compare the two quantiles
anova(fit2,fit3) 

# lets compare the 1st, 2nd and 3rd quantiles
anova(fit1,fit2,fit3,joint=FALSE) 


# Alternatively, fitting can be done in one call:
fit <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, method = "fn", tau = 1:4/5, data = barro)
summary(fit) # this looks that the 20th, 40th, 60th, and 80th percentiles
plot(summary(fit)) # see what it looks like



# iid - conditional distributions are independent and identically distributed (iid).

# generate data with constant variance 
x <- seq(0,100,length.out = 100)        # independent variable
int <- 6                                # true intercept of mean
slope <- 0.1                            # true slope
set.seed(1)                             # make the next line reproducible
e <- rnorm(100,mean = 0, sd = 5)        # normal random error with constant variance
y <- int + slope*x + e                    # dependent variable
dat <- data.frame(x,y)                  # make into datafame
plot(dat$x, dat$y)                      # plot data

# note that summary gives confidence intervals but not standard errors
summary(rq(y~x, tau=c(0.1, 0.5, 0.9)))

# add SE option to statement
#iid - assumes iid distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="iid"))

#nid - assumes non-identical conditional distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="nid"))

#boot - no distribution assumption. Note that this is most similar to the iid, which is the appropriate choice for data with this distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="boot"))


# generate data with non-constant variance (adapted from https://data.library.virginia.edu/getting-started-with-quantile-regression/)
#generate data with nonconstant variance
x <- seq(0,100,length.out = 100)        # independent variable
varfunc <- 0.1 + 0.1*x                  # non-constant variance
int <- 6                                # true intercept
slope <- 0.1                            # true slope
set.seed(1)                             # make the next line reproducible
e <- rnorm(100,mean = 0, sd = varfunc)  # normal random error with non-constant variance
y <- int + slope*x + e                    # dependent variable
dat <- data.frame(x,y)                  # make into datafame
plot(dat$x, dat$y)                      # plot data

# add SE option to statement
# note that summary gives confidence intervals but not standard errors
summary(rq(y~x, tau=c(0.1, 0.5, 0.9)))

# add SE option to statement
#iid - assumes iid distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="iid"))

#nid - assumes non-identical conditional distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="nid"))

#boot - no distribution assumption. Note that this is most similar to the nid, which is the appropriate choice for data with this distribution
plot(summary(rq(y~x, tau=c(0.1, 0.5, 0.9)), se="boot"))


####
####  Evaluating Goodness-of-Fit
####

# generate data with non-constant variance (adapted from https://data.library.virginia.edu/getting-started-with-quantile-regression/)

x <- seq(0,100,length.out = 100)        # independent variable
varfunc <- 0.1 + 0.1*x                  # non-constant variance
int <- 6                                # true intercept
slope <- 0.1                            # true slope
set.seed(1)                             # make the next line reproducible
e <- rnorm(100,mean = 0, sd = varfunc)  # normal random error with non-constant variance
y <- int + slope*x + e                    # dependent variable
dat <- data.frame(x,y)                  # make into datafame
plot(dat$x, dat$y)                      # plot data


####  R1 metric

# fitted model for median quantile
fitmod5 <- rq(y ~ x, tau = .5, data = dat)


# null model for median quantile
nullmod5<-rq(y~1, data=dat, tau=0.5)

#weighted residuals
nullres5<-resid(nullmod5)
fitres5 <- resid(fitmod5)

# calculate R1 metric as 1 minus ratio of summed absolute residuals from fitted and null model
1-(sum(abs(fitres5))/sum(abs(nullres5)))

# or more simply, pull these sums from the model object
1 - fitmod5$rho/nullmod5$rho

# compare to R2
summary(lm(y ~ x, data = dat))$r.squared

# plot data with mean and median regression lines. note in this case they are the same line but goodness of fit differs (QR < OLS)
plot(dat$x, dat$y)
abline(fitmod5, col="red")
abline(lm(y ~ x, data = dat), col="blue")


####  difference in AIC

# delta AIC method for QR model
AIC(nullmod5)[1]-AIC(fitmod5)[1]

#compare to OLS
AIC(lm(y ~ 1, data = dat))[1]-AIC(lm(y ~ x, data = dat))[1]

# notice the dAIC for the QR model vs. null is larger than that for the OLS vs. null


# cumulative sum of the gradient vector (eigenvector) (Qtools package)

# specify QR model object, significance level, B= number of montecarlo samples
GOFTest(fitmod5, alpha = 0.05, B = 1000)

#try with dataset with smaller variance
x <- seq(0,100,length.out = 100)        # independent variable
varfunc <- 0.1 + 0.01*x                  # non-constant variance
int <- 6                                # true intercept
slope <- 0.1                            # true slope
set.seed(1)                             # make the next line reproducible
e <- rnorm(100 ,mean = 0, sd = varfunc)  # normal random error with non-constant variance
y <- int + slope*x + e                    # dependent variable
dat <- data.frame(x,y)                  # make into datafame
plot(x,y)

fitmod5b <- rq(y ~ x, tau = .5, data = dat)
GOFTest(fitmod5b, alpha = 0.05, B = 1000)

#try with dataset with smaller variance and larger sample size
x <- seq(0,1000,length.out = 1000)        # independent variable
varfunc <- 0.1 + 0.01*x                  # non-constant variance
int <- 6                                # true intercept
slope <- 0.1                            # true slope
set.seed(1)                             # make the next line reproducible
e <- rnorm(1000 ,mean = 0, sd = varfunc)  # normal random error with non-constant variance
y <- int + slope*x + e                    # dependent variable
dat <- data.frame(x,y)                  # make into datafame
plot(x,y)

fitmod5c <- rq(y ~ x, tau = .5, data = dat)
GOFTest(fitmod5c, alpha = 0.05, B = 1000)



####
#### Plotting Quantile Regression in ggplot
####

data(engel)
m <- ggplot(engel, aes(income, foodexp)) + geom_point()
m + geom_quantile()
m + geom_quantile(quantiles = 0.5)

q10 <- seq(0.05, 0.95, by = 0.05)
m + geom_quantile(aes(colour = ..quantile..),quantiles = q10)


#  Nonlinear Quantile Regression  (adapted from: https://rpubs.com/MarkusLoew/10676)

data(Mammals) #Observations on the maximal running speed of mammal species and their body mass.
attach(Mammals)
x <- log(weight)
y <- log(speed)
mamdat<-as.data.frame(cbind(x,y))

# plot data
plot(x,y, xlab="Weight in log(Kg)", ylab="Speed in log(Km/hour)")

# start with just a simple quadratic approach
# y ~ a * x^2 + b * x  + c

my.equation <- y ~ a * x^2 + b * x  + c

# fit the equation to the data via "non-linear least squares"
# choose some good starting values for parameter estimation
nls.fit <- nls(my.equation,
               data = mamdat,
               start = list(a = 2, b = 3, c = 5))

# look at the result
summary(nls.fit)

# create a dummy range of that we use to predict speed from our fitted model
predict_range <- data.frame(x = seq(-5, 9, length = 250))

# calculate for each x-range value the corresponding y-range
my.line <- within(predict_range, y <- predict(nls.fit, newdata = predict_range))

# add the line to the existing graph
# This line represents the "mean" fit, no quantile regression involved
# plot data
plot(x,y, xlab="Weight in log(Kg)", ylab="Speed in log(Km/hour)")
lines(y ~ x, data = my.line, col = "red")


# Non-linear quantile regression
# aiming for the upper 99% quantile
my.rq <- nlrq(my.equation,
           data = mamdat,
           start = list(a = 2, b = 2, c = 5),
           tau = .99)
summary(my.rq)

# calculating the values from the model
my.line95 <- within(predict_range, 
             y <- predict(my.rq, 
                            newdata = predict_range))
plot(x,y, xlab="Weight in log(Kg)", ylab="Speed in log(Km/hour)")
lines(y ~ x, data = my.line, col = "red")
lines(y ~ x, data = my.line95, col = "blue")


# Esterase data
data(esterase)

# Fit quantiles 0.25 and 0.75 - number of dithered samples
fit1 <- rq.counts(Count ~ Esterase, tau = 0.25, data = esterase)
coef(fit1)
fit2 <- rq.counts(Count ~ Esterase, tau = 0.75, data = esterase)
coef(fit2)


# Plot
with(esterase, plot(Count ~ Esterase))
lines(esterase$Esterase, fit1$fitted.values, col = "blue")
lines(esterase$Esterase, fit2$fitted.values, col = "red")
legend(8, 1000, lty = c(1,1), col = c("blue", "red"), legend = c("tau = 0.25","tau = 0.75"))

