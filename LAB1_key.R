
#  NRES 746, Lab 1  ----------------------
##  University of Nevada, Reno                        
##  Computational algorithms and review of statistics  


# LAB 1 answer key  --------------------

# this is just one solution- there are lots of other possibilities

## exercise 1a  ----------------------

CoefVar <- function(x){    
  cv <- sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
  return(cv)     
}
#CoefVar(c(2,3,4,3,2,3,4))   # test your function


# Testing your code ------------------------ 
#   Explore the "trees" dataset 

#?trees

summary(trees)    # learn more about the data

trees$Height    # extract the "Height" column from the trees dataset.

CoefVar(trees$Height)    # run your new function!


library(ggplot2)

## exercise 1b ---------------------

    # ggplot version
DrawLine <- function(x,y){ 
  coefs <- coef(lm(y~x))
  df = data.frame(x=x,y=y)
  plot = ggplot(df,aes(x,y)) +
    geom_point() +
    geom_smooth(method="lm") +
    theme_classic()
  print(plot)
  return(coefs)
}

      # base R version
DrawLine_base <- function(x,y){ 
  plot(y~x,xlab="predictor",ylab="response")
  mod = lm(y~x)
  coefs <- coef(mod)
  abline(mod)
  return(coefs)
}


#DrawLine(trees$Height,trees$Volume)

# ?faithful
# summary(faithful)

## ggplot version:
DrawLine(faithful$waiting,faithful$eruptions)    # test your function using the old faithful eruptions data

## base R version:
DrawLine_base(faithful$waiting,faithful$eruptions)


## exercise 1c --------------------

    # ggplot version
DrawLine2 <- function(x,y,smooth=TRUE,span=1){ 
  if(smooth) coefs <- loess(y~x,span=span) else coefs <- coef(lm(y~x))
  df = data.frame(x=x,y=y)
  plot = ggplot(df,aes(x,y)) +
    geom_point() +
    theme_classic()
  if(smooth) plot = plot + geom_smooth(method="loess",span=span) else plot = plot + geom_smooth(method="lm")
  print(plot)
  return(coefs)
}

      # base R version
DrawLine2_base <- function(x,y,smooth=TRUE,span=1){  
  if(smooth==FALSE){
    plot(y~x)
    mod <- lm(y~x)
    coefs <- coef(mod)
    abline(coefs,lwd=2)
  }else{
    coefs <- loess(y~x,span=span)
    scatter.smooth(y~x,span=span,lpars=list(lwd=2))
  }
  return(coefs)     
}

xvec <- c(1:10)
yvec <- rnorm(length(xvec),c(2:6,7:3),2)
DrawLine2(xvec,yvec,smooth=T,span=1)    # run your new function!

DrawLine2(x=trees$Height,y=trees$Volume,smooth=F)

DrawLine2(faithful$waiting,faithful$eruptions,smooth=T,span=1)    # test using the old faithful eruptions data

DrawLine2(faithful$waiting,faithful$eruptions,smooth=T,span=1.9)


## exercise 1d: CLT function  ------------------

CLTdemo <- function(N=10,min=0,max=1,visualize=T){
  lots <- 1000       # placeholder representing infinity 
  many_samples = replicate(lots,runif(N,min,max))  # generate lots of samples!
  many_xbars = colMeans(many_samples)  # compute lots of sample means
  if(visualize){
    layout(matrix(1:2,nrow=1))
    hist(many_xbars,freq=F)    # distribution of sample means.
    qqnorm(many_xbars)   # normal q-q plot
  }
  return(shapiro.test(many_xbars)$p.value)
}


CLTdemo(3,10,20)    # run your new function!

 ## you might want to replicate this many times...

replicate(5,CLTdemo(3,10,20,visualize=F))   # for example, replicate 5 times



## Exercise 1e  ----------------------------

sizes <- seq(2,15)

evidence <- sapply(sizes, function(this) replicate(100,CLTdemo(this,10,20,visualize = F) ) )
pass_rate <- apply(evidence,2,function(this) sum(this>0.1)/length(this)  )

par(mfrow=c(1,1))
names(pass_rate) <- sizes
barplot(pass_rate,xlab="sample size",ylab="Proportion passed tests")


## exercise 2a: predict ozone concentration  ------------------

# summary(model2)
# ?airquality

CtoF <- function(c){ c*9/5 + 32 }
nd <- data.frame(
  Temp = CtoF(c(26,30)),
  Solar.R = 200,
  Wind = 9
)
diff(predict(model2,nd))   # difference of 12.88 ppb Ozone 


## exercise 2b: predict ozone concentration  ------------------

##  written answer:  the null hypothesis is that the true regression coefficient for the population of interest is zero. 


## exercise 2c: diagnostic  ------------------

summary(model2)

par(mfrow=c(3,2),mai=c(1,1,0,0))
plot(model2)

plot(air_cleaned$Wind,residuals(model2))
plot(air_cleaned$Temp,residuals(model2))

### Short answer: the diagnostic plots look okay except there seems to be some evidence for nonlinear response of Ozone to Wind



## exercise 2d: variable importance  ------------------

# summary(model2)
# ?airquality

library(ggeffects)
library(cowplot)

r1 = predict_response(model2, terms = "Wind")
r2 = predict_response(model2, terms = c("Temp","Solar.R"))

a = plot(r1,show_data = T, jitter = 0.1)
b = plot(r2,show_data = T, jitter = 0.1)

plot_grid(a,b)

anova(model2)
drop1(model2,test="Chisq")

model3 = lm(Ozone~scale(Wind)+scale(Solar.R)*scale(Temp),data=air_cleaned)
summary(model3)

## short answer: all variables are important and it is difficult to tell which is 'most important'. Standardized coefficients indicate that temperature may be the strongest effect, followed by wind. Interestingly, the main effect of the 'Temp' variable was 'non-significant' in model 2 but the effect of "Temp acted primarily through the interaction term. Note: when interaction terms are included in a model it does NOT usually make sense to interpret the main effect without also considering the interaction term. 


# Exercise 3: brute force z-tests! --------------------


#   modified from lecture...

ztest_bruteforce <- function(d, mu, sigma){
  return_list <- list() # bundle results to return using a list object
  return_list$Xbar <- mean(d)  # Compute the sample statistic
  N <- length(d)   # compute sample size
  lots <- 1000       # stand-in for infinity [compute lots of replicate samples under the null hypothesis] 
  return_list$nulldist <- replicate(lots, mean(rnorm(N,mu,sigma)) )     # Simulate random sample and generate sampling distribution
  return_list$p_value <- sum(return_list$nulldist<=return_list$Xbar)/lots   # how many of these are more extreme than the sample statistic?
  return(return_list)
}

# mu = 3.8     # testing code (commented out)
# sigma = 0.9
# mysample = c(3.14,3.27,2.56,3.77,3.34,4.32,3.84,2.19,5.24,3.09)
# ztest_bruteforce(d = mysample, mu, sigma )

# zstat = (mean(mysample)-mu)/(sd(mysample)/sqrt(length(mysample)))
# pnorm(zstat)


## Question 3a  -----------------------

library(tidyverse)

my_ztest <- function(d, mu, sigma, alternative='b'){
  return_list <- list() # bundle results to return using a list object
  return_list$Xbar <- mean(d)  # Compute the sample statistic
  N <- length(d)   # compute sample size
  lots <- 10000       # stand-in for infinity 
  return_list$nulldist <- replicate(lots, mean(rnorm(N,mu,sigma)) )     # generate sampling distribution
  if(!(alternative%in%c('b','g','l'))) ~ stop("alternative must have one of the values 'l', 'g', or 'b'")
  dist_from_null <- return_list$nulldist-mu
  return_list$p_value <- case_when(    # you can use if-then-else, but the case_when method is a little cleaner
    alternative == 'l' ~ sum(dist_from_null<=(return_list$Xbar-mu))/lots,
    alternative == 'g' ~ sum(dist_from_null>=(return_list$Xbar-mu))/lots,
    alternative == 'b' ~ sum(abs(dist_from_null)>=abs(return_list$Xbar-mu))/lots
  )
  return(return_list)
}


mu = 4     # testing code
sigma = 2
mysample = c(3.14,3.27,2.56,3.77,3.34,4.32,3.84,2.19,5.24,3.09)
test = my_ztest(mysample, mu, sigma,alternative = 'l' )
str(test)
test$p_value

   # also try using the code from lecture to compare against a "real" z-test!


x_vs_null <- function(x, nulldata){
  to_return <- list()   # initialize object to return
  to_return$xbar <- mean(x)   # sample mean
  N <- length(x)   # compute sample size
  lots <- 1000                 # set the number of replicate samples
  to_return$nulldist <- replicate(lots,mean(sample(nulldata,N,replace=T)))    # sampling distribution under null
  to_return$p_value <- sum(to_return$nulldist<=to_return$xbar)/lots  # how many of these samples are more extreme than the sample statistic?
  to_return
}


nulldata=c(2.2,3.86,6.39,4.6,3.43,5.16,4.36,4.22,6.31,4.61,5.13,4.12,4.64,4.03,5.01,7.33,5.35,4.7,2.82,4.87,3.87,5.95,5.28,4.02,3.58,4.03,5.38,5.5,3.07,3.29,3.45,5.25,5.7,1.26,5.28,4.19,4.76,4.2,4.81,2.5)
mysample = c(3.14,3.27,2.56,3.77,3.34,4.32,3.84,2.19,5.24,3.09)
test = x_vs_null(mysample, nulldata )
test$p_value


# Exercise 4: Bootstrapping regression coefficients! ---------------

# first, grab the code from the lecture (bootstrapping R-squared)
    
Rsquared <- function(df,responsevar="Volume"){    # univariate models only- interaction and multiple regression not implemented here
  response <- df[,responsevar]       # extract the response variable
  varnames <- setdiff(names(df),responsevar)      # extract all column names that are not the response    
  rsq <- numeric(length(varnames))        # named storage vector
  names(rsq) <- varnames               
  for(i in names(rsq)){         # loop through predictors
      model <- lm(as.formula(paste0(responsevar,"~",i)),df)       # regress response on predictor
      rsq[i] <- summary(model)$r.square       # extract R-squared statistic   1-RSS/TSS
  }
  return(rsq)     
}

boot_sample <- function(df,fun,nboot,y){
  t(replicate(nboot,fun(df[sample(1:nrow(df),replace=T),],responsevar=y) ))  #  randomly sample observations with replacement and generate statistics from the bootstrapped sample
}


## exercise 4a  -------------------

RegressionCoefs <- function(df,responsevar){    # univariate models only- interaction and multiple regression not implemented here
  response <- df[,responsevar]       # extract the response variable
  names <- names(df)                  
  coefs <- numeric(length(names))        # named storage vector
  names(coefs) <- names(df)               
  coefs <- coefs[names(coefs)!=responsevar]           # assume that all columns that are not the response variable are possible predictor variables
  for(i in names(coefs)){         # loop through predictors
      predictor <- df[,i]                  # extract this predictor
      model <- lm(response~predictor)       # regress response on predictor
      coefs[i] <- coefficients(model)["predictor"]       # extract slope term
  }
  return(coefs)     
}



RegressionCoefs(trees,"Volume")   # should return two regression coefficients


## exercise 4b  -----------------

 # note: this is copied from the lecture code...
boot_sample <- function(df,fun,nboot,y){
  t(replicate(nboot,fun(df[sample(1:nrow(df),replace=T),],responsevar=y) ))  #  randomly sample observations with replacement and generate statistics from the bootstrapped sample
}


BootCoefs <- function(df,fun,n_samples,responsevar){
  boot <- boot_sample(df,fun,n_samples,responsevar)   # generate test statistics (coefs) for 1000 bootstrap samples
  apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))       # summarize the quantiles to generate confidence intervals for each predictor variable
}



BootCoefs(df=trees,fun=RegressionCoefs,1000,responsevar="Volume")

df <- mtcars[,c(1,3,4,6)]
responsevar="mpg"
BootCoefs(df,RegressionCoefs,1000,responsevar)


## Exercise 4c  ---------------------

BootCoefs(trees,RegressionCoefs,responsevar="Volume")

data.frame(
  Girth=confint(lm(Volume~Girth,trees))["Girth",],    # compare with lm() 
  Height=confint(lm(Volume~Height,trees))["Height",]    # compare with lm() 
)

# There are slight differences due to stochasticity, and the bootstrapped coefficient intervals may be slightly wider, but the answers are essentially the same. 

# end of lab 1


