
############################################################
####                                                    ####  
####  NRES 746, Lecture 4                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Likelihood                                        ####
####     Assessing the probability of the data          ####
####     under a known data-generating model            ####
############################################################



#########
# try an exponential model

Deterministic_component <- function(xvals,a,b){
  yexp <- a*exp(b*xvals)        # deterministic exponential decay (assuming b is negative)
  return(yexp)
}

DataGenerator_exp <- function(xvals,params){
  yexp <- Deterministic_component(xvals,params$a,params$b)  # get signal
  yvals <- rnorm(length(yexp),yexp,sqrt(params$c))     # add noise
  return(yvals)
}



###########
# generate data under an assumed process model

xvals=mtcars$disp    # xvals same as data (this is a "fixed effect", so there is no random component here- we can't really "sample" x values)
params <- list()  
params$a=30             # set model parameters arbitrarily (eyeballing to the data)
params$b=-0.005   # = 1/200
params$c=1

yvals <- DataGenerator_exp(xvals,params)

plot(yvals~xvals)      # plot the simulated data



##########
# assess goodness-of-fit of a known data-generating model

PlotRangeOfPlausibleData <- function(xvals,params,reps){ 
  samplesize <- length(xvals)
  results <- array(0,dim=c(samplesize,reps))   # storage array for results
  for(i in 1:reps){
    yvals <- DataGenerator_exp(xvals,params)
    results[,i] <- yvals
  }
      # now make a boxplot of the results
  boxplot(lapply(1:nrow(results), function(i) results[i,]),at=xvals, xaxt="n",main="Plausible data under this model",ylab="mpg",xlab="Displacement",boxwex=6)
  cleanseq <- (seq(0,max(round(xvals/100)),length=(max(round(xvals/100)))+1))*100
  axis(1,at=cleanseq,labels = cleanseq)    # label the x axis properly
  
}


reps <- 1000    # number of replicate datasets to generate

PlotRangeOfPlausibleData(xvals,params,reps)    # run the function to visualize the range of data that could be produced under this model


############
# finally, overlay the real data to evaluate goodness of fit!

real_yvals <- mtcars$mpg
PlotRangeOfPlausibleData(xvals,params,reps)
points(xvals,real_yvals,pch=20,cex=3,col="green")


##########
# now change the parameters and see if the data fit to the model

params$a=40       # was 30
params$b=-0.001   # was 0.005

    
PlotRangeOfPlausibleData(xvals,params,reps)
points(xvals,real_yvals,pch=20,cex=3,col="green")    # overlay the real data


#######
# try again- select a new set of parameters

params$a=33       # was 40
params$b=-0.002   # was 0.001

    
PlotRangeOfPlausibleData(xvals,params,reps)
points(xvals,real_yvals,pch=20,cex=3,col="green")    # overlay the real data


#############
# Work with likelihood!

obs.data <- mtcars[1,c("mpg","disp")]    # for simplicity, consider only the first observation
obs.data


############
# "best fit" parameters from above
############

params <- list()    # set up empty list to store parameters
params$a=33            # fill the list with the "best fit" parameter set from above    
params$b=-0.002   
params$c=1

params

expected_val <- Deterministic_component(obs.data$disp,params$a,params$b)   
expected_val      # expected mpg for the first observation in the "mtcars" dataset


###########
# Visualize the likelihood of this single observation. 

mean = expected_val   # expected (mean) value for this observation, given the "known" data generating model
stdev = sqrt(params$c)    # standard deviation

curve(dnorm(x,mean,stdev),10,30,xlab="possible vals for mpg under the specified model",ylab="probability density")   # probability of all plausible mpg values under the data generating model.  
abline(v=obs.data$mpg,col="red",lwd=2)    # overlay the observed data


###########
# compute the likelihood of the first observation

likelihood = dnorm(obs.data$mpg,mean,stdev)
likelihood


###########
# Visualize the likelihood of two observations. 

obs.data <- mtcars[c(1,3),c("mpg","disp")]
obs.data

par(mfrow=c(1,2))  # set up graphics!

for(i in 1:nrow(obs.data)){
  curve(dnorm(x,Deterministic_component(obs.data$disp[i],params$a,params$b),sqrt(params$c)),10,30,xlab="mpg",ylab="probability density")   # probability density
  abline(v=obs.data$mpg[i],col="red",lwd=2)
}


##########
# compute the likelihood of observing BOTH data points

Likelihood <- dnorm(obs.data$mpg[1],Deterministic_component(obs.data$disp[1],params$a,params$b),sqrt(params$c)) *
              dnorm(obs.data$mpg[2],Deterministic_component(obs.data$disp[2],params$a,params$b),sqrt(params$c))  
Likelihood


##############
# and now... four observations!

obs.data <- mtcars[c(1,3,4,5),c("mpg","disp")]
obs.data

par(mfrow=c(2,2))  # set up graphics!

for(i in 1:nrow(obs.data)){
  curve(dnorm(x,Deterministic_component(obs.data$disp[i],params$a,params$b),sqrt(params$c)),10,30,xlab="mpg",ylab="probability density")   # probability density
  abline(v=obs.data$mpg[i],col="red",lwd=2)
}



##########
# compute the likelihood of observing all four data points

Likelihood <- 1     # initialize the likelihood
for(i in 1:nrow(obs.data)){
  Likelihood <- Likelihood * dnorm(obs.data$mpg[i],Deterministic_component(obs.data$disp[i],params$a,params$b),sqrt(params$c))
}
Likelihood


# Alternatively, we can use the "prod" function in R

Likelihood <- prod(dnorm(obs.data$mpg,Deterministic_component(obs.data$disp,params$a,params$b),sqrt(params$c)))
Likelihood


##########
# Finally, compute the likelihood of ALL data points in the entire data set, using the "prod()" function

full.data <- mtcars[,c("mpg","disp")]
Likelihood <- prod(dnorm(full.data$mpg,Deterministic_component(full.data$disp,params$a,params$b),sqrt(params$c)))
Likelihood


##########
# Compute the log-likelihood (much easier to work with!)

Log.Likelihood <- sum(dnorm(full.data$mpg,Deterministic_component(full.data$disp,params$a,params$b),sqrt(params$c),log=TRUE)) 
Log.Likelihood  
exp(Log.Likelihood)   # we can convert back to likelihood if we want...


##########
# Example likelihood function!

# Arguments:
#   params: bundled vector of free parameters for the known data-generating model
#   df: a data frame that holds the observed data
#   yvar: the name of the response variable (ancillary)
#   xvar: the name of the predictor variable (ancillary)

LogLikFunction <- function(params,df,yvar,xvar){
  LogLik <- sum(dnorm(df[,yvar],Deterministic_component(df[,xvar],params['a'],params['b']),sqrt(params['c']),log=TRUE))
  return(LogLik)
}
LogLikFunction(unlist(params),df=mtcars,yvar="mpg",xvar="disp")


############
# Use numerical optimization methods to identify the maximum likelihood estimate (and the likelihood at the MLE)

MLE <- optim(fn=LogLikFunction,par=unlist(params),df=mtcars,yvar="mpg",xvar="disp",control=list(fnscale=-1))  # note, the control param is set so that "optim" maximizes rather than minimizes the Log-likelihood. 


MLE$par   # maximum likelihood parameter estimates


MLE$value   # log likelihood for the best model


##############
# visualize goodness-of-fit for the best model

bestParams <- as.list(MLE$par)   # extract the MLE parameter vals

xvals <- mtcars$disp
yvals <- mtcars$mpg
PlotRangeOfPlausibleData(xvals,bestParams,1000)
points(xvals,yvals,pch=20,cex=3,col="green")


##############
# Visualize a "slice" of the likeihood function

upperval <- -1/1000
lowerval <- -1/200
allvals <- seq(lowerval,upperval,length=1000)
likelihood_slice <- numeric(1000)   # set up storage vector! 
newParams <- bestParams 
for(i in c(1:length(allvals))){
  newParams$b <- allvals[i]
  likelihood_slice[i] <- exp(LogLikFunction(unlist(newParams),mtcars,"mpg","disp"))    # get the data likelihood across slice of parameter space
}

plot(allvals,likelihood_slice,type="l",main="Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Likelihood")


#################
# Work with log-likelihood instead...

upperval <- -1/1000
lowerval <- -1/200
allvals <- seq(lowerval,upperval,length=1000)
loglikelihood_slice <- numeric(1000)   # set up storage vector! 
newParams <- bestParams 
for(i in c(1:length(allvals))){
  newParams$b <- allvals[i]
  loglikelihood_slice[i] <- LogLikFunction(unlist(newParams),mtcars,"mpg","disp")    # get the data likelihood across slice of parameter space
}

plot(allvals,loglikelihood_slice,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")


###########
# zoom in closer to the MLE

upperval <- -1/550
lowerval <- -1/350
allvals <- seq(lowerval,upperval,length=1000)
loglikelihood_slice <- numeric(1000)   # set up storage vector! 
newParams <- bestParams 
for(i in c(1:length(allvals))){
  newParams$b <- allvals[i]
  loglikelihood_slice[i] <- LogLikFunction(unlist(newParams),mtcars,"mpg","disp")    # get the data likelihood across slice of parameter space
}

plot(allvals,loglikelihood_slice,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")


##############
# what parameter values are within 2 log likelihood units of the best value?

bestVal <- bestParams$b
bestVal



plot(allvals,loglikelihood_slice,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")
abline(v=bestVal,lwd=3,col="blue")
abline(h=(MLE$value-2))


##############
# Generate an approximate 95% confidence interval for the "b" parameter

reasonable_parameter_values <- allvals[loglikelihood_slice>=(MLE$value-2)]
min(reasonable_parameter_values)
max(reasonable_parameter_values)
plot(allvals,loglikelihood_slice,type="l",main="Log Likelihood slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")
abline(v=bestVal,lwd=3,col="blue")
abline(h=(MLE$value-2),lty=2)
abline(v=min(reasonable_parameter_values),lwd=1,col="blue")
abline(v=max(reasonable_parameter_values),lwd=1,col="blue")


###############
# A better confidence interval, using the likelihood "profile"


##########
# first, visualize the likelihood surface in 2 dimensions

upperval_b <- -1/800
lowerval_b <- -1/300

upperval_a <- 50
lowerval_a <- 5

allvals_a <- seq(lowerval_a,upperval_a,length=500)
allvals_b <- seq(lowerval_b,upperval_b,length=500)

loglikelihood_surface <- matrix(0,nrow=500,ncol=500)   # set up storage matrix! 

newParams <- bestParams 
for(i in 1:length(allvals_a)){  # loop through possible a params
  newParams$a <- allvals_a[i]
  for(j in 1:length(allvals_b)){    # loop through possible b params
    newParams$b <- allvals_b[j]
    loglikelihood_surface[i,j] <- LogLikFunction(unlist(newParams),mtcars,"mpg","disp")    # get the data likelihood across slice of parameter space
  }
}

image(x=allvals_a,y=allvals_b,z=loglikelihood_surface,zlim=c(-100,-75),col=topo.colors(12))


##########
# add a contour line, assuming deviances follow a chi-squared distribution

conf95 <- qchisq(0.95,2)/2  # this evaluates to around 3. Since we are varying freely across 2 dimensions, we use chisq with 2 degrees of freedom
image(x=allvals_a,y=allvals_b,z=loglikelihood_surface,zlim=c(-100,-75),col=topo.colors(12))
contour(x=allvals_a,y=allvals_b,z=loglikelihood_surface,levels=(MLE$value-conf95),add=TRUE,lwd=3,col=gray(0.3))


#############
# visualize likelhood profile!

              ### A parameter
profile_A <- apply(loglikelihood_surface,1,max)
reasonable_parameter_values_A <- allvals_a[profile_A >=(MLE$value-qchisq(0.95,1)/2)]
min(reasonable_parameter_values_A)
max(reasonable_parameter_values_A)
plot(allvals_a,profile_A,type="l",main="Log Likelihood profile",xlab="Parameter Slice for \'a\'",ylab="Log-Likelihood")
abline(v=MLE$par["a"],lwd=3,col="blue")
abline(v=min(reasonable_parameter_values_A),lwd=1,col="blue")
abline(v=max(reasonable_parameter_values_A),lwd=1,col="blue")


############
# profile for the b parameter... 

profile_B <- apply(loglikelihood_surface,2,max)
reasonable_parameter_values_B <- allvals_b[profile_B >=(MLE$value-qchisq(0.95,1)/2)]
min(reasonable_parameter_values_B)
max(reasonable_parameter_values_B)
plot(allvals_b,profile_B,type="l",main="Log Likelihood profile",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")
abline(v=MLE$par["b"],lwd=3,col="blue")
abline(v=min(reasonable_parameter_values_B),lwd=1,col="blue")
abline(v=max(reasonable_parameter_values_B),lwd=1,col="blue")


###############
# demo: likelihood ratio test

curve(dchisq(x,2),0,10,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=2")


curve(dchisq(x,2),0,10,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=2")
abline(v=qchisq(0.95,2),col="red",lwd=2)


curve(dchisq(x,1),0,5,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=1")
abline(v=qchisq(0.95,1),col="red",lwd=2)

