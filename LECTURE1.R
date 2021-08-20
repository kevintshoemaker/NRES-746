
############################################################
####                                                    ####  
####  NRES 746, Lecture 1                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Computational algorithms vs standard statistics   ####
############################################################



###################
# SALMON EXAMPLE (made-up!)

population.mean = 4.5
population.sd = 0.9

my.sample = c(3.14,3.27,2.56,3.77,3.34,4.32,3.84,2.19,5.24,3.09)

sample.size <- length(my.sample)     # determine sample size   

obs.samplemean = mean(my.sample)     # note the equal sign as assignment operator

## visualize the population of conventional-raised salmon

curve(dnorm(x,population.mean,population.sd),0,10,
      xlab="Body mass (kg)",ylab="Probability density")

## now overlay this on the observed data

hist(my.sample,freq=F,
     xlab="Body mass (kg)",ylab="Probability density",main="",
     xlim=c(0,10))
curve(dnorm(x,population.mean,population.sd),0,10,
      col="red",lwd=2,add=T)
abline(v=obs.samplemean,col="blue",lwd=3)


################
# Perform standard z-test
################

library(BSDA)
z.test(x=my.sample,mu=population.mean, sigma.x=population.sd,alternative = "less")


############
# alternative z-test

std.err = population.sd/sqrt(sample.size)

curve(dnorm(x,population.mean,std.err),0,10,     # visualize the sampling distribution under null hypothesis
      xlab="Body mass (kg)",ylab="Probability density")     # versus the observed sample mean
abline(v=obs.samplemean,col="blue",lwd=3)

p.val = pnorm(obs.samplemean,population.mean,std.err)
p.val     # this is the same as the p value from the z-test above...


######################
   # ALTERNATIVE ALGORITHMIC APPROACH!
######################

#############
# Simulate the STATISTICAL POPULATION under the null hypothesis
#############

infinity <- 1000000  # large number approximating infinity 

popData_null <- rnorm(n=infinity,mean=population.mean,sd=population.sd)    # the statistical "population" of interest (under null model w no 'treatment' effect)


#############
# Draw a SAMPLE from that null data
#############

null.sample <- sample(popData_null,size=sample.size)    # use R's native "sample()" function to sample from the null distribution

round(null.sample,2)
null.samplemean <- mean(null.sample)  
null.samplemean    # here is one sample mean that we can generate under the null hypothesis


#################
# Repeat this process using a FOR loop
#################

n.samples <- 1000                 # set the number of replicate samples to generate
null.samplemeans <- numeric(n.samples)       # initialize a storage vector for sample means under the null hypothesis

for(i in 1:n.samples){            # for each replicate... 
  this.nullsample <- sample(popData_null,size=sample.size)      # draw a sample of body masses assuming no treatment effect       
  null.samplemeans[i] <- mean(this.nullsample)           # compute and store the sampling distribution produced under the null hypothesis
}

hist(null.samplemeans,xlim=c(0,10))       # plot out the sampling distribution
abline(v=obs.samplemean,col="green",lwd=3)     # overlay the observed sample statistic. 


############
# Generate a p-value algorithmically!!
############

ordered_means <- sort(null.samplemeans)       # sort the vector of null sample means
more_extreme <- length(which(ordered_means<=obs.samplemean))       # how many of these sampling errors equal or exceed the "error" represented by the observed statistic?
p_value <- more_extreme/n.samples       # compute a p-value! 
p_value    


#############
# Develop a function that wraps up all the above steps into one!
#############

z.test.algorithm <- function(sample, pop.mean, pop.sd){
  
  #############
  # Compute the sample statistic
  #############
  
  observed_mean <- mean(sample)
  
  sample.size <- length(observed_mean)   # compute sample size

  #################
  # Generate SAMPLING DISTRIBUTION
  #################
  
  reps <- 1000                 # set the number of replicate samples
  null_dist <- numeric(reps)       # initialize a storage structure for sampling distribution
  
  for(i in 1:reps){            # for each replicate... 
    nullsamp <- rnorm(10,pop.mean,pop.sd)      # draw a sample assuming no treatment effect       
    null_dist[i] <- mean(nullsamp)           # compute and store the sampling error produced under the null hypothesis
  }
  
  more.extreme <- length(which(null_dist<=observed_mean))       # how many of these sampling errors equal or exceed the sample statistic?
  p_value <- more.extreme/reps
  
  to_return <- list()   # initialize object to return
  
  to_return$null_dist <- null_dist
  to_return$p_value <- p_value
  to_return$observed_mean <- observed_mean
  
  return(to_return)

}

ztest <- z.test.algorithm(sample = my.sample, pop.mean=population.mean, pop.sd=population.sd )   # try to run the new function

ztest$p_value     # get the p_value

hist(ztest$null_dist)       # plot out all the sampling errors under the null hypothesis as a histogram
abline(v=ztest$observed_mean,col="green",lwd=3)     # indicate the observed sample statistic. 


#############
# Start with a made-up data frame!
#############

df <- data.frame(
  A = c(175, 168, 168, 190, 156, 181, 182, 175, 174, 179),
  B = c(185, 169, 173, 173, 188, 186, 175, 174, 179, 180) 
)

summary(df)    # summarize! 

sample.size <- length(df$A)     # determine sample size    

#######
# Get data in proper format

reshape_df <- data.frame(                # "reshape" the data frame so each observation gets its own row (standard 'tidy' format)
  Treatment = rep(c("A","B"),each=sample.size),
  Mass = c(df$A,df$B),
  stringsAsFactors = T
)


########
# Alternative (commented out)- using the 'tidyverse'

# library(tidyr)
# reshape_df <- pivot_longer(df,everything(),names_to = "Treatment",values_to="Mass")


plot(Mass~Treatment, data=reshape_df)    # explore/visualize the data

#######
# Compute the observed difference between group means

observed_dif <- mean(reshape_df$Mass[reshape_df$Treatment=="A"])	- mean(reshape_df$Mass[reshape_df$Treatment=="B"])



##################
# NON-PARAMETRIC T-TEST -- PERMUTATION TEST
##################

reps <- 5000            # Define the number of permutations to run (number of replicates)
null_difs <- numeric(reps)   # initialize storage variable
for (i in 1:reps){			# For each replicate:		
  newGroup <- reshape_df$Treatment[sample(c(1:nrow(reshape_df)))]			   # randomly shuffle the observed data with respect to treatment group
	dif <- mean(reshape_df$Mass[newGroup=="A"])	- mean(reshape_df$Mass[newGroup=="B"])	   #  compute the difference between the group means after reshuffling the data
	null_difs[i] <- dif	    # store this value in a vector
}
hist(null_difs)    # Plot a histogram of null differences between group A and group B under the null hypothesis (sampling errors)
abline(v=observed_dif,col="green",lwd=3)   # Add a vertical line to the plot to indicate the observed difference


########
# Compute a p-value based on the permutation test, just like we did before (except now 2-tailed)!
########

more_extreme <- length(which(abs(null_difs)>=abs(observed_dif)))
p_value <- more_extreme/reps  
p_value


#############
# Develop a function that performs a permutation-t-test!
#############

t.test.permutation <- function(dat = reshape_df, group = "Treatment", value = "Mass" ){
  
  #############
  # Compute the sample statistic
  #############
  
  indexA <- which(dat[,group]=="A")     # rows representing treatment A
  indexB <- which(dat[,group]=="B")     # rows representing treatment B
  observed_dif <- mean(dat[indexA,value]) - mean(dat[indexB,value])
  
  reps <- 5000            # Define the number of permutations to run (number of replicates)
  null_difs <- numeric(reps)   # initialize storage variable
  for (i in 1:reps){			# For each replicate:		
    newGroup <- reshape_df$Treatment[sample(c(1:nrow(reshape_df)))]			   # randomly shuffle the observed data with respect to treatment group
  	dif <- mean(reshape_df$Mass[newGroup=="A"])	- mean(reshape_df$Mass[newGroup=="B"])	   #  compute the difference between the group means after reshuffling the data
  	null_difs[i] <- dif	    # store this value in a vector
  }
  
  more_extreme <- length(which(abs(null_difs)>=abs(observed_dif)))
  p_value <- more_extreme/reps  
  
  to_return <- list()   # initialize object to return
  
  to_return$null_difs <- null_difs
  to_return$p_value <- p_value
  to_return$observed_dif <- observed_dif
  
  return(to_return)
  
}

my.ttest <- t.test.permutation()   # use default values for all function arguments

my.ttest$p_value

hist(my.ttest$null_difs)    # Plot a histogram of null differences between group A and group B under the null hypothesis (sampling errors)
abline(v=my.ttest$observed_dif,col="green",lwd=3)   # Add a vertical line to the plot to indicate the observed difference



##############
# Demonstration: bootstrapping a confidence interval!

## use the "trees" dataset in R:

head(trees)   # use help(trees) for more information


#########
# Basic data exploration

plot(trees$Volume~trees$Height, main = 'Black Cherry Tree Height/Volume Relationship', xlab = 'Height', ylab = 'Volume', pch = 16, col ='blue')
plot(trees$Volume~trees$Girth, main = 'Black Cherry Tree Girth/Volume Relationship', xlab = 'Girth', ylab = 'Volume', pch = 16, col ='red')


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


#########
# test the function to see if it works!

stat <- Rsquared(trees,"Volume")
stat


############
# new function to generate "bootstrap" samples from a data frame

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


##########
# Generate a few bootstrapped samples!

boot <- boot_sample(df=trees,statfunc=Rsquared,n_samples=10,n_stats=2)       # generate test stats from lots of bootstrapped samples
colnames(boot) <- names(stat)         # name the columns to recall which predictor variables they represent

boot
stat


#############
# use bootstrapping to generate confidence intervals for R-squared statistic!

boot <- boot_sample(df=trees,statfunc=Rsquared,n_samples=1000,n_stats=2)   # generate test statistics (Rsquared vals) for 1000 bootstrap samples
confint <- apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))       # summarize the quantiles to generate confidence intervals for each predictor variable
colnames(confint) <- names(stat)
t(confint)


