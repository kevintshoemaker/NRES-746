
#  NRES 746, Lecture 1                          
#   University of Nevada, Reno                        
#   Bespoke algorithms for inference   


# SALMON EXAMPLE (made-up!) ------------------

pop_mean = 4.5
pop_sd = 0.9

mysample = c(3.14,3.27,2.56,3.77,3.34,4.32,3.84,2.19,5.24,3.09)

myN <- length(mysample)     # determine sample size   

mysamplemean = mean(mysample)     # note the equal sign as alternative assignment operator

## visualize the population of conventional-raised salmon  -------------------

curve(dnorm(x,pop_mean,pop_sd),0,10,
      xlab="Body mass (kg)",ylab="Probability density",lwd=2)

### now overlay this on the observed data  --------------------

hist(mysample,freq=F,
     xlab="Body mass (kg)",ylab="Probability density",main="",
     xlim=c(0,10))
curve(dnorm(x,pop_mean,pop_sd),0,10,
      col="red",lwd=2,add=T)
abline(v=mysamplemean,col="blue",lwd=3)


# Perform "canned" z-test  ----------------------------

library(BSDA)
z.test(x=mysample,mu=pop_mean, sigma.x=pop_sd,alternative = "less")


# alternative z-test in base R (no packages)  -----------------------

pop_se = pop_sd/sqrt(myN)   # standard deviation for sample means drawn from the null population

curve(dnorm(x,pop_mean,pop_se),0,10,     # visualize the sampling distribution under null hypothesis
      xlab="Body mass (kg)",ylab="Probability density")     # versus the observed sample mean
abline(v=mysamplemean,col="blue",lwd=3)

p.val = pnorm(mysamplemean,pop_mean,pop_se)    # note that neither pop_mean or pop_se are random variables- they are known with certainty. Therefore we can use a normal distribution (the known data distribution under the null hypothesis, as specified above) to define the sampling error. 
p.val     # this is the same as the p value from the z-test above...


# ALTERNATIVE ALGORITHMIC Z-TEST! ----------------------

## Simulate the STATISTICAL POPULATION under the null hypothesis -----------------

lots <- 1000000  # large number filling in for infinity 

null_population <- rnorm(n=lots,mean=pop_mean,sd=pop_sd)    # the statistical "population" of interest (under null model w no 'treatment' effect)


## Draw a SAMPLE from the null population ----------------

null_sample <- sample(null_population,size=myN)    # use R's native "sample()" function to sample randomly from the null distribution

round(null_sample,2)
null_statistic <- mean(null_sample)  
null_statistic    # here is one sample mean that we can generate under the null hypothesis


## Repeat this process using a FOR loop ----------------------

null_replicates <- 1000                 # set the number of replicate samples to generate  (yet another number intended to approximate infinity!)
null_statistics <- numeric(null_replicates)       # initialize a storage vector for sample means under the null hypothesis

for(i in 1:null_replicates){            # for each replicate... 
  null_sample <- sample(null_population,size=myN)      # draw a random sample of body masses assuming no treatment effect       
  null_statistics[i] <- mean(null_sample)           # compute and store the sampling distribution produced under the null hypothesis
}

hist(null_statistics,xlim=c(0,10))       # plot out the sampling distribution using base R plotting
abline(v=mysamplemean,col="green",lwd=3)     # overlay the observed sample statistic  


## Generate a p-value  --------------------------

more_extreme <- length(which(null_statistics<=mysamplemean))       # how many of these sampling errors equal or exceed the "extremeness" of the observed statistic?
p_value <- more_extreme/null_replicates       # compute a p-value! 
p_value    


# Develop a function that wraps up all the above steps into one! ------------------

ztest_bruteforce <- function(d, mu, sigma){
  return_list <- list() # bundle results to return using a list object
  return_list$Xbar <- mean(d)  # Compute the sample statistic
  N <- length(d)   # compute sample size
  reps <- 1000       # stand-in for infinity [compute lots of replicate samples under the null hypothesis] 
  return_list$nulldist <- numeric(reps)       # initialize a storage structure for sampling distribution
  for(i in 1:reps){            # for each replicate... 
    return_list$nulldist[i] <- mean(rnorm(N,mu,sigma))           # Simulate random sample and generate sampling distribution
  }
  return_list$p_value <- sum(return_list$nulldist<=return_list$Xbar)/reps   # how many of these are more extreme than the sample statistic?
  return(return_list)
}

ztest <- ztest_bruteforce(d = mysample, mu=4.5, sigma=0.9 )   # try to run the new function

ztest$p_value     # get the p_value

hist(ztest$nulldist,xlim=range(c(ztest$Xbar,ztest$nulldist)))       # plot out all the samples under the null hypothesis as a histogram
abline(v=ztest$Xbar,col="green",lwd=3)     # indicate the observed sample statistic. 


# Nonparametric t-test (permutation test) ------------------------

## Start with a made-up data frame ---------------------

df <- data.frame(
  trtA = c(175, 168, 168, 190, 156, 181, 182, 175, 174, 179),
  Control = c(185, 169, 173, 173, 188, 186, 175, 174, 179, 180) 
)

summary(df)    # summarize! 

N <- nrow(df)     # determine sample size N   

# Get data in proper format

reshape_df <- data.frame(                # "reshape" the data frame so each observation gets its own row (standard 'tidy' format)
  Treatment = rep(factor(c("trtA","Control"),levels=c("Control","trtA")),each=N),
  Mass = c(df$trtA,df$Control),
  stringsAsFactors = T
)

plot(Mass~Treatment, data=reshape_df)    # explore/visualize the data

## Compute the observed difference between group means  -----------------

observed_dif <-  diff(with(reshape_df,tapply(Mass,Treatment,mean)))



## Run permutation t-test ----------------

reps <- 5000            # Number of replicates - another number representing infinity!
null_difs <- numeric(reps)   # initialize storage variable
for (i in 1:reps){			# For each replicate:		
  newGroup <- reshape_df$Treatment[sample(c(1:nrow(reshape_df)))]			   # assign each observation a random treatment group
	null_difs[i] <- mean(reshape_df$Mass[newGroup=="trtA"])	- mean(reshape_df$Mass[newGroup=="Control"])	   #  compute the difference between the group means after reshuffling the data
}
hist(null_difs)    # Plot a histogram of null differences between group A and group B under the null hypothesis (sampling errors)
abline(v=observed_dif,col="green",lwd=3)   # Add a vertical line to the plot to indicate the observed difference


## Compute a p-value based on the permutation test -------------------

  #just like we did before (except now 2-tailed)!

p_value <- sum(abs(null_difs)>=abs(observed_dif))/reps 
p_value


## Develop a function that performs a permutation-t-test! -----------------------

ttest_permutation <- function(dat = reshape_df, group = "Treatment", value = "Mass" ){
  to_return <- list()   # initialize object to return
  to_return$observed_dif = diff(tapply(dat[[value]],dat[[group]],mean)  )            # Compute the sample statistic
  to_return$null_difs <- numeric(reps)   # initialize storage variable
  for (i in 1:reps){			# For each replicate:		
  	to_return$null_difs[i] <- diff(tapply(dat[[value]],sample(dat[[group]]),mean)  ) 	    # compute and store sample stat after permuting the treatments
  }
  
  to_return$p_value <- sum(abs(null_difs)>=abs(observed_dif))/reps  
  return(to_return)
}

mytest <- ttest_permutation(reshape_df)   # use default values for all function arguments

mytest$p_value

hist(mytest$null_difs)    # Plot a histogram of null differences between group A and group B under the null hypothesis (sampling errors)
abline(v=mytest$observed_dif,col="green",lwd=3)   # Add a vertical line to the plot to indicate the observed difference



# Demonstration: bootstrapping a confidence interval! ---------------------

## use the "trees" dataset in R:
head(trees)   # use help(trees) for more information


## Data exploration  --------------------

plot(trees$Volume~trees$Height, main = 'Black Cherry Tree Height/Volume Relationship', xlab = 'Height', ylab = 'Volume', pch = 16, col ='blue')
plot(trees$Volume~trees$Girth, main = 'Black Cherry Tree Girth/Volume Relationship', xlab = 'Girth', ylab = 'Volume', pch = 16, col ='red')


## Function for returning a vector of R-squared statistics  -------------
     # Function for returning a vector of R-squared statistics from models regressing a response variable on multiple possible predictor variables
   # here we assume that all columns in the input data frame that are NOT the response variable are potential predictor variables.

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


# test the function to see if it works!  ----------------------

stat <- Rsquared(trees,"Volume")
stat


Rsquared(mtcars,"mpg")


# new function to generate multiple "bootstrap" estimates of a test statistic  ----------------

boot_sample <- function(df,fun,nboot,y){

  t(replicate(nboot,fun(df[sample(1:nrow(df),replace=T),],responsevar=y) ))  #  randomly sample observations with replacement and generate statistics from the bootstrapped sample
}


# Generate a few bootstrapped samples!  ------------------

boot <- boot_sample(trees,Rsquared,500, "Volume")       # generate test stats from lots of bootstrapped samples

summary(boot)
stat


# use bootstrapping to generate confidence intervals for R-squared statistic!  ----------------

boot <- boot_sample(trees,Rsquared,1000, "Volume")    # generate test statistics (Rsquared vals) for 1000 bootstrap samples
confint <- apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))       # summarize the quantiles to generate confidence intervals for each predictor variable
confint


