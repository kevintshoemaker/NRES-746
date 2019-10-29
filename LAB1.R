
############################################################
####                                                    ####  
####  NRES 746, Lab 1                                   ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Computational algorithms and standard statistics  ####
############################################################


######################
# LAB 1 (possible) answers
######################

#############
## exercise 1a

CoefVar <- function(vector){    
  cv <- sd(vector,na.rm=TRUE)/mean(vector,na.rm=TRUE)
  return(cv)     
}



answer1a <- CoefVar(trees$Height)    # run your new function!
answer1a


check_result(pass_if(answer1a),grader_args=list(user_quo = quo({test <- trees$Height;CoefVar <-            function(vector){sd(vector)/mean(vector)};CoefVar(vector=test)})))

test_result(grader_args=list(user_quo = quo({a=10;a}),solution_quo=quo(10)))    # use this syntax!




#############
## exercise 1b

DrawLine <- function(x,y){  
  plot(y~x)
  mod <- lm(y~x)
  coefs <- coef(mod)
  abline(mod)
  return(coefs)
  
}



answer1b <- DrawLine(trees$Height,trees$Volume)    # run your new function!

DrawLine(faithful$waiting,faithful$eruptions)    # test using the old faithful eruptions data




#############
## exercise 1c

DrawLine2 <- function(x,y,smooth=TRUE,span=1){  
  plot(y~x)
  if(smooth==FALSE){
    mod <- lm(y~x)
    coefs <- coef(mod)
    abline(coefs)
  }else{
    coefs <- loess(y~x,span=span)
    scatter.smooth(y~x,span=span)
  }
  return(coefs)     
}

xvec <- c(1:10)
set.seed(100) 
yvec <- rnorm(length(xvec),c(2:6,7:3),2)
answer1c <- DrawLine2(xvec,yvec,smooth=T,span=0.5)    # run your new function!

DrawLine2(x=trees$Height,y=trees$Volume,smooth=F,span=NA)

DrawLine2(faithful$waiting,faithful$eruptions,span=.5)    # test using the old faithful eruptions data

DrawLine2(faithful$waiting,faithful$eruptions,span=.1)


#############
## exercise 1d

####################
# CENTRAL LIMIT THEOREM function
####################


CLTdemo <- function(n.samples=1000,sample.size=10,min=10,max=20){
  lots <- 100000       # number approximating infinity      
  datafountain <- runif(lots,min,max)      # here we define the full set of possible random numbers to draw random samples from (the population of potential data)
  
  #######
  # Draw multiple samples from the pool (population) of possible data.  
  
  samplemean <- numeric(n.samples)     # set up storage vector
  i=1
  for(i in 1:n.samples){      # for each replicate (independent random sample)
    sample <- sample(datafountain,sample.size)   # draw an independent random sample from the population of interest
    samplemean[i] <- mean(sample)    # compute and record the sample mean
  }
  
  
  par(mfrow=c(1,2))
  hist(datafountain,freq=F,ylim=c(0,1),main="",xlab="Value")      # plot out the distribution of sample means
  hist(samplemean,freq=F,add=T,col="red")    # overlay the distribution of the underlying data from which we are drawing samples.
  legend("topleft",fill=c("white","red"),legend=c("data population","sample means"),bty="n")
  qqnorm(samplemean)
  
  
  out <- shapiro.test(samplemean)
  
  return(out)
}




answer1d <- CLTdemo(n.samples=5000,sample.size=4,min=10,max=20)    # run your new function!



##############
# Exercise 1e

reps <- 50
sizes <- seq(3,10,1)
prop <- numeric(length(sizes))
i=sizes[1]
counter=1
for(i in sizes){
  okay <- logical(reps)
  j=1
  for(j in 1:reps){
    temp <- CLTdemo(n.samples=1000,sample.size=i,min=10,max=20,plot=FALSE)    # run your new function!
    okay[j] <- temp$p.value<=0.05
  }
  num <- length(which(okay))
  prop[counter] <- num/reps   # proportion rejected
  counter=counter+1
}

par(mfrow=c(1,1))
names(prop) <- sizes
barplot(prop,xlab="sample size",ylab="Normality rejected (proportion)")





#############
# CHALLENGE 4: brute force t-tests!
#############

##### from lecture...

t.test.algorithm <- function(dat = reshape_df, group = "Treatment", value = "Mass" ){
  
  #############
  # Compute the sample statistic
  #############
  
  indexA <- which(dat[,group]=="A")     # rows representing treatment A
  indexB <- which(dat[,group]=="B")     # rows representing treatment B
  observed_dif <- mean(dat[indexA,value]) - mean(dat[indexB,value])
  
  sample.size <- length(indexA)
  
  #############
  # Simulate the STATISTICAL POPULATION under the null hypothesis
  #############
  
  lots <- 1000000  # large number approximating infinity 
  
  popMean_null <- mean(dat[,value])           # assume groups A and B come from a population with common mean 
  popSD_null <- sd(dat[,value])                      # and common standard deviation... 
  popData_null <- rnorm(n=lots,mean=popMean_null,sd=popSD_null)    # the statistical "population" of interest (under null model w no treatment effect)

  #################
  # Repeat sampling process (sampling from population) using a FOR loop
  #################
  
  reps <- 1000                 # set the number of replicates
  null_difs <- numeric(reps)       # initialize a storage structure to hold one anomaly (sampling error) per replicate
  
  for(i in 1:reps){            # for each replicate... 
    sampleA <- sample(popData_null,size=sample.size)      # draw a sample assuming no treatment effect       
    sampleB <- sample(popData_null,size=sample.size)      # draw a sample assuming no treatment effect (again!)
    null_difs[i] <- mean(sampleA)-mean(sampleB)           # compute and store the sampling error produced under the null hypothesis
  }
  
  ordered_difs <- sort(abs(null_difs))       # sort the vector of sampling errors 
  higher_anomaly <- length(which(ordered_difs>=abs(observed_dif)))       # how many of these sampling errors equal or exceed the sample statistic?
  p_value <- higher_anomaly/reps
  
  to_return <- list()   # initialize object to return
  
  to_return$null_difs <- null_difs
  to_return$p_value <- p_value
  to_return$observed_dif <- observed_dif
  
  return(to_return)

}



###############
# Question 4a

t.test.onetail <- function(dat = reshape_df, group = "Treatment", value = "Mass" ){
  
  #############
  # Compute the sample statistic
  #############
  
  indexA <- which(dat[,group]=="A")     # rows representing treatment A
  indexB <- which(dat[,group]=="B")     # rows representing treatment B
  observed_dif <- mean(dat[indexB,value]) - mean(dat[indexA,value])
  
  sample.size <- length(indexA)
  
  #############
  # Simulate the STATISTICAL POPULATION under the null hypothesis
  #############
  
  lots <- 1000000  # large number approximating infinity 
  
  popMean_null <- mean(dat[,value])           # assume groups A and B come from a population with common mean 
  popSD_null <- sd(dat[,value])                      # and common standard deviation... 
  popData_null <- rnorm(n=lots,mean=popMean_null,sd=popSD_null)    # the statistical "population" of interest (under null model w no treatment effect)

  #################
  # Repeat sampling process (sampling from population) using a FOR loop
  #################
  
  reps <- 1000                 # set the number of replicates
  null_difs <- numeric(reps)       # initialize a storage structure to hold one anomaly (sampling error) per replicate
  
  for(i in 1:reps){            # for each replicate... 
    sampleA <- sample(popData_null,size=sample.size)      # draw a sample assuming no treatment effect       
    sampleB <- sample(popData_null,size=sample.size)      # draw a sample assuming no treatment effect (again!)
    null_difs[i] <- mean(sampleB)-mean(sampleA)           # compute and store the sampling error produced under the null hypothesis
  }
  
  ordered_difs <- sort(null_difs)       # sort the vector of sampling errors 
  higher_anomaly <- length(which(ordered_difs>=observed_dif))       # how many of these sampling errors equal or exceed the sample statistic?
  p_value <- higher_anomaly/reps
  
  to_return <- list()   # initialize object to return
  
  to_return$null_difs <- null_difs
  to_return$p_value <- p_value
  to_return$observed_dif <- observed_dif
  
  return(to_return)

}



df <- data.frame(
  TreatmentA = c(175, 168, 168, 190, 156, 181, 182, 175, 174, 179),
  TreatmentB = c(185, 169, 173, 173, 188, 186, 175, 174, 179, 180) 
)


reshapefunc <- function(df,varname){
  sample.size <- nrow(df)     # determine sample size 

  reshape_df <- data.frame(                # "reshape" the data frame so each observation gets its own row (standard format)
    Treatment = rep(c("A","B"),each=sample.size),
    Value = c(df[,1],df[,2])
  )
  
  names(reshape_df)[2] <- varname   # rename the value column to whatever the user specifies
  return(reshape_df)
}


dfnew <- reshapefunc(df,"Mass")

ttest_onetail <- t.test.onetail(dat=dfnew,group = "Treatment", value = "Mass" )
ttest_onetail$p_value

ttest_twotail <- t.test.algorithm(dat=dfnew,group = "Treatment", value = "Mass" )    # try two-tailed test
ttest_twotail$p_value

t.test(df$TreatmentA,df$TreatmentB,alternative = "less")          # compare with true t-test
t.test(df$TreatmentA,df$TreatmentB,alternative = "two.sided")



df <- data.frame(
      TreatmentA = c(105, 118, 119, 112, 109),
      TreatmentB = c(138, 130, 150, 123, 144)
    )
df2<-reshapefunc(df,varname="Mass")

answer4a <- t.test.onetail(dat=df2,group = "Treatment", value = "Mass")



###############
# A possible answer to question 4b

t.test.vardif <- function(dat = df.vardif2, group = "Treatment", value = "Mass" ){
  
  #############
  # Compute the sample statistic
  #############
  
  indexA <- which(dat[,group]=="A")     # rows representing treatment A
  indexB <- which(dat[,group]=="B")     # rows representing treatment B
  observed_dif <- mean(dat[indexA,value]) - mean(dat[indexB,value])
  
  sample.size <- length(indexA)
  
  varA <- sd(dat[indexA,value])    # estimate variance of each group separately
  varB <- sd(dat[indexB,value])
  
  #############
  # Simulate the STATISTICAL POPULATION under the null hypothesis
  #############
  
  lots <- 1000000  # large number approximating infinity 
  
  popMean_null <- mean(dat[,value])           # assume groups A and B come from a population with common mean 
  popData_nullA <- rnorm(n=lots,mean=popMean_null,sd=varA)    # the statistical "population" of interest (under null model w no treatment effect but allowing variances to differ)
  popData_nullB <- rnorm(n=lots,mean=popMean_null,sd=varB)    # the statistical "population" of interest (under null model w no treatment effect but allowing variances to differ)
  #################
  # Repeat sampling process (sampling from population) using a FOR loop
  #################
  
  reps <- 1000                 # set the number of replicates
  null_difs <- numeric(reps)       # initialize a storage structure to hold one anomaly (sampling error) per replicate
  
  for(i in 1:reps){            # for each replicate... 
    sampleA <- sample(popData_nullA,size=sample.size)      # draw a sample assuming no treatment effect       
    sampleB <- sample(popData_nullB,size=sample.size)      # draw a sample assuming no treatment effect (again!)
    null_difs[i] <- mean(sampleB)-mean(sampleA)           # compute and store the sampling error produced under the null hypothesis
  }
  
  ordered_difs <- sort(abs(null_difs))       # sort the vector of sampling errors 
  higher_anomaly <- length(which(ordered_difs>=abs(observed_dif)))       # how many of these sampling errors equal or exceed the sample statistic?
  p_value <- higher_anomaly/reps
  
  to_return <- list()   # initialize object to return
  
  to_return$null_difs <- null_difs
  to_return$p_value <- p_value
  to_return$observed_dif <- observed_dif
  
  return(to_return)

}



###########
# Test data for unequal variances (prior to reshaping)...

df.vardif <- data.frame(
  TreatmentA = c(135, 128, 139, 122, 126, 121, 128, 135, 134, 129),
  TreatmentB = c(215, 69, 143, 153, 218, 186, 125, 98, 271, 340)
)

summary(df.vardif)    # summarize!


df.vardif2 <- reshapefunc(df.vardif,"Mass")

ttest_vardif <- t.test.vardif(dat=df.vardif2, group = "Treatment", value = "Mass")    # test new function
ttest_vardif$p_value

ttest_varequal <- t.test.algorithm(dat=df.vardif2)    # versus original function
ttest_varequal$p_value

t.test(df.vardif$TreatmentA,df.vardif$TreatmentB, var.equal = FALSE, alternative="two.sided")   # versus built-in function




df <- data.frame(
      TreatmentA = c(1135, 1128, 1139, 1122, 1126),
      TreatmentB = c(1915, 69, 3143, 53, 1818)
    )
df2<-reshapefunc(df,varname="Mass")
    
answer4b <- t.test.vardif(dat=df2,group = "Treatment", value = "Mass")



###############
# A possible answer to question 4c

t.test.ndif <- function(dat = df.ndif2, group = "Treatment", value = "Mass" ){
  
  #############
  # Compute the sample statistic
  #############
  
  levs <- levels(dat[,group])
  
  indexA <- which(dat[,group]==levs[1])     # rows representing treatment A
  indexB <- which(dat[,group]==levs[2])     # rows representing treatment B
  observed_dif <- mean(dat[indexA,value],na.rm=TRUE) - mean(dat[indexB,value],na.rm=T)
  
  sample.size <- tapply(dat[,value],dat[,group],function(t) length(which(!is.na(t))))

  
  varA <- sd(dat[indexA,value],na.rm = T)    # estimate variance of each group separately
  varB <- sd(dat[indexB,value],na.rm=T)
  
  #############
  # Simulate the STATISTICAL POPULATION under the null hypothesis
  #############
  
  lots <- 1000000  # large number approximating infinity 
  
  popMean_null <- mean(dat[,value],na.rm=TRUE)           # assume groups A and B come from a population with common mean 
  popData_nullA <- rnorm(n=lots,mean=popMean_null,sd=varA)    # the statistical "population" of interest (under null model w no treatment effect but allowing variances to differ)
  popData_nullB <- rnorm(n=lots,mean=popMean_null,sd=varB)    # the statistical "population" of interest (under null model w no treatment effect but allowing variances to differ)
  #################
  # Repeat sampling process (sampling from population) using a FOR loop
  #################
  
  reps <- 1000                 # set the number of replicates
  null_difs <- numeric(reps)       # initialize a storage structure to hold one anomaly (sampling error) per replicate
  
  for(i in 1:reps){            # for each replicate... 
    sampleA <- sample(popData_nullA,size=sample.size[1])      # draw a sample assuming no treatment effect       
    sampleB <- sample(popData_nullB,size=sample.size[2])      # draw a sample assuming no treatment effect (again!)
    null_difs[i] <- mean(sampleB)-mean(sampleA)           # compute and store the sampling error produced under the null hypothesis
  }
  
  ordered_difs <- sort(abs(null_difs))       # sort the vector of sampling errors 
  higher_anomaly <- length(which(ordered_difs>=abs(observed_dif)))       # how many of these sampling errors equal or exceed the sample statistic?
  p_value <- higher_anomaly/reps
  
  to_return <- list()   # initialize object to return
  
  to_return$null_difs <- null_difs
  to_return$p_value <- p_value
  to_return$observed_dif <- observed_dif
  
  return(to_return)

}



###########
# Test data for unequal sample sizes (before reshaping)...

df.ndif <- data.frame(
  TreatmentA = c(135, 128, 139, 122, 126, 121, 128, 135, 134, 129, 134, 125, 130, 132, 125),
  TreatmentB = c(98, 271, 340, rep(NA,times=12))
)

summary(df.ndif)    # summarize!


df.ndif2 <- reshapefunc(df.ndif,"Mass")

ttest_ndif <- t.test.ndif(dat=df.ndif2)    # test new function
ttest_ndif$p_value


df <- data.frame(
      TreatmentA = c(135, 128, 139, 122, 126,rep(NA,3)),
      TreatmentB = c(315, 69, 143, 53, 818,4,55,190)
    )
df2<-reshapefunc(df,varname="Mass")
answer4c <- t.test.ndif(dat=df2,group = "Treatment", value = "Mass")


#####################
# CHALLENGE 5: Bootstrapping regression coefficients!
#####################

# first, grab the code from the lecture (bootstrapping R-squared)

#########
# Function for returning a vector of R-squared statistics from models regressing a response variable on multiple possible predictor variables
   # here we assume that all columns in the input data frame that are NOT the response variable are potential predictor variables.

Rsquared <- function(df,responsevar="Volume"){    # univariate models only- interaction and multiple regression not implemented here
  response <- df[,responsevar]       # extract the response variable
  names <- names(df)                  
  rsq <- numeric(length(names))        # named storage vector
  names(rsq) <- names(df)               
  rsq <- rsq[names(rsq)!=responsevar]           # assume that all columns that are not the response variable are possible predictor variables
  for(i in names(rsq)){         # loop through predictors
      predictor <- df[,i]                  # extract this predictor
      model <- lm(response~predictor)       # regress response on predictor
      rsq[i] <- summary(model)$r.square       # extract R-squared statistic
  }
  return(rsq)     
}

boot_sample <- function(df,statfunc,n_samples,n_stats,responsevar="Volume"){
  indices <- c(1:nrow(df))
  output <- matrix(NA,nrow=n_samples,ncol=n_stats)        # storage object- to store a single bootstrapped sample from the original data
  
  for(i in 1:n_samples){              # for each bootstrap replicate:
    boot_rows <- sample(indices,size=nrow(df),replace=T)         # randomly sample observations with replacement
    newdf <- df[boot_rows,]                       # dataframe of bootstrapped observations
    output[i,] <- statfunc(newdf,responsevar)                 # generate statistics from the bootstrapped sample  (e.g., compute Rsquared after regressing y on all possible x variables)
  }
  return(output)
}

# boot <- boot_sample(df=trees,statfunc=Rsquared,n_samples=1000,n_stats=2)   # generate test statistics (Rsquared vals) for 1000 bootstrap samples
# confint <- apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))       # summarize the quantiles to generate confidence intervals for each predictor variable
# colnames(confint) <- names(stat)
# t(confint)



#############
## exercise 5a

RegressionCoefs <- function(df=trees,responsevar="Volume"){    # univariate models only- interaction and multiple regression not implemented here
  response <- df[,responsevar]       # extract the response variable
  names <- names(df)                  
  coefs <- numeric(length(names))        # named storage vector
  names(coefs) <- names(df)               
  coefs <- coefs[names(coefs)!=responsevar]           # assume that all columns that are not the response variable are possible predictor variables
  i=names(coefs)[1]
  for(i in names(coefs)){         # loop through predictors
      predictor <- df[,i]                  # extract this predictor
      model <- lm(response~predictor)       # regress response on predictor
      coefs[i] <- coefficients(model)["predictor"]       # extract slope term
  }
  return(coefs)     
}



RegressionCoefs(df=trees,responsevar="Volume")


df <- mtcars[,c(1,3,4,6)]
answer5a <- coefficients(lm(df$mpg~df$hp))[2]


#############
## exercise 5b

boot_sample <- function(df=trees,statfunc=Rsquared,n_samples=1000,n_stats=2,responsevar="Volume"){
  indices <- c(1:nrow(df))
  output <- matrix(NA,nrow=n_samples,ncol=n_stats)        # storage object- to store a single bootstrapped sample from the original data
  
  for(i in 1:n_samples){              # for each bootstrap replicate:
    boot_rows <- sample(indices,size=nrow(df),replace=T)         # randomly sample observations with replacement
    newdf <- df[boot_rows,]                       # dataframe of bootstrapped observations
    output[i,] <- statfunc(newdf,responsevar)                 # generate statistics from the bootstrapped sample  (e.g., compute Rsquared after regressing y on all possible x variables)
  }
  return(output)
}


BootCoefs <- function(df=trees,statfunc=RegressionCoefs,n_samples=1000,n_stats=2,responsevar="Volume"){    # univariate models only- interaction and multiple regression not implemented here
  boot <- boot_sample(df=df,statfunc=statfunc,n_samples=n_samples,n_stats=ncol(df)-1,responsevar=responsevar)   # generate test statistics (Rsquared vals) for 1000 bootstrap samples
  confint <- apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))       # summarize the quantiles to generate confidence intervals for each predictor variable
  colnames(confint) <- paste("stat",c(1:n_stats),sep="")
  return(t(confint))     
}



BootCoefs(df=trees,statfunc=RegressionCoefs,n_samples=1000,n_stats=2,responsevar="Volume")


df <- mtcars[,c(1,3,4,6)]
responsevar="mpg"
answer5b <- BootCoefs(
              df=df,
              statfunc=RegressionCoefs,
              n_samples=1000,
              n_stats=3,
              responsevar=responsevar
            )

#answer5b <- round(confint(lm(df$mpg~df$disp))[2,],2)[1]


############
# Exercise 5c

BootCoefs(df=trees,statfunc=RegressionCoefs,n_samples=1000,n_stats=2,responsevar="Volume")

confint.lm(lm(Volume~Girth,trees))    # compare with lm() 



###########
# end of lab 1
###########

