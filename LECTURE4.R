
#  NRES 746, Lecture 4     -------------------------------                  
#   University of Nevada, Reno                       
#   Likelihood                                       
#     Assessing the probability of the data          
#     under a fully specified data-generating model        


# Demo: using data simulation to make inferences ----------------

data(mtcars)    # use the 'mtcars' data set as an example 

# ?mtcars

plot(mpg~disp, data = mtcars, las = 1, pch = 16, xlab = "Displacement (cu. in.)", ylab = "Miles/Gallon")   # visualize the relationship


# try an exponential model

mu_func <- function(x,a,b){
  a*exp(b*x)        # deterministic exponential decline (assuming b is negative)
}

DataGenerator <- function(x,params){
  a=params[1]; b=params[2]; sigma=params[3]
  rnorm(length(x),mu_func(x,a,b),sigma)     
}


## generate data under an assumed process model -----------------

xvals=mtcars$disp    # xvals same as data (there is no random component here- we can't really "sample" x values)
params <- c(  
  a = 30,             # set model parameters arbitrarily (eyeballing to the data) (see Bolker book)
  b = -0.005,   # = 1/200
  sigma=2
)

yvals <- DataGenerator(xvals,params)

plot(yvals~xvals,las = 1, pch = 16, xlab = "Displacement (cu. in.)", ylab = "Miles/Gallon")      # plot the simulated data



## assess goodness-of-fit of a known data-generating model --------------------

VisualizeModelWithData <- function(x,params){ 
  lots=1000
  plotdat = data.frame(xseq=round(seq(min(x),max(x),length=25)))
  results = replicate(lots,DataGenerator(plotdat$xseq,params))
      # now make a boxplot of the results
  bounds=sapply(1:nrow(results), function(t) c(lb=quantile(results[t,],0.025,names=F), ub=quantile(results[t,],0.975,names=F), mean=mean(results[t,]))  )
  plotdat = cbind(plotdat,t(bounds))
  plot = ggplot(plotdat,aes(x=xseq,y=mean)) +
    geom_ribbon(aes(ymin=lb,ymax=ub),fill="gray") +
    geom_path(lwd=2) + 
    geom_point(data=mtcars,aes(x=disp,y=mpg),col="darkgreen") +
    labs(x="Displacement",y="mpg",title ="Plausible data under this model" ) +
    theme_classic()
  print(plot)
}


VisualizeModelWithData(xvals,params)    # run the function to visualize the range of data that could be produced under this model


# now change the parameters and see if the data fit to the model

params["a"]=45       # was 30
params["b"]=-0.005   
    
VisualizeModelWithData(xvals,params)   


# try again- select a new set of parameters

params["b"]=-0.002   # was -0.005
params["sigma"]=3    # was 2
    
VisualizeModelWithData(xvals,params)   

# Work with likelihood! ---------------------

obs.data <- mtcars[1,c("mpg","disp")]    # for simplicity, consider only the first observation
obs.data


## "best fit" parameters from above   ----------------

params["a"]=35            # fill the list with the "best fit" parameter set from above (this is still just an educated guess) 
params["b"]= -0.0029   
params["sigma"]=0.95

expected_val <- mu_func(obs.data$disp,params["a"],params["b"])   
expected_val      # expected mpg for the first observation in the "mtcars" dataset


## Visualize the likelihood of this single observation.  ------------------------

mu_star = expected_val   # expected (mean) value for this observation, given the data generating model
sigma_star = params['sigma']    # standard deviation

curve(dnorm(x,mu_star,sigma_star),10,30,xlab="Response outcome under the specified model",ylab="probability density")   # probability of all plausible mpg values under the data generating model.  
abline(v=obs.data$mpg,col="red",lwd=2)    # overlay the observed data


## compute the likelihood of the first observation  ---------------------

likelihood1 = dnorm(obs.data$mpg,mu_star,sigma_star)
likelihood1


## Visualize the likelihood of two observations. ---------------

obs.data <- mtcars[c(1,3),c("mpg","disp")]
obs.data

par(mfrow=c(1,2))  # set up graphics!

for(i in 1:nrow(obs.data)){
  curve(dnorm(x,mu_func(obs.data$disp[i],params['a'],params['b']),params['sigma']),10,30,xlab="mpg",ylab="probability density")   # probability density
  abline(v=obs.data$mpg[i],col="red",lwd=2)
}


## compute the likelihood of observing BOTH data points  ----------------

Likelihood <- dnorm(obs.data$mpg[1],mu_func(obs.data$disp[1],params['a'],params['b']),params['sigma']) *
              dnorm(obs.data$mpg[2],mu_func(obs.data$disp[2],params['a'],params['b']),params['sigma'])  
Likelihood

prod(dnorm(obs.data$mpg,mu_func(obs.data$disp,params['a'],params['b']),params['sigma']))

# and now... four observations!

obs.data <- mtcars[c(1,3,4,5),c("mpg","disp")]
obs.data

par(mfrow=c(2,2))  # set up graphics!

for(i in 1:nrow(obs.data)){
  curve(dnorm(x,mu_func(obs.data$disp[i],params['a'],params['b']),params['sigma']),10,30,xlab="mpg",ylab="probability density")   # probability density
  abline(v=obs.data$mpg[i],col="red",lwd=2)
}



    # compute the likelihood of observing all four data points
prod(dnorm(obs.data$mpg,mu_func(obs.data$disp,params['a'],params['b']),params['sigma']))


       # Finally, compute the likelihood of ALL data points in the entire data set, using the "prod()" function
full.data <- mtcars[,c("mpg","disp")]

prod(dnorm(full.data$mpg,mu_func(full.data$disp,params['a'],params['b']),params['sigma']))


## Compute the log-likelihood (easier to work with!)  -----------------------

l <- sum(dnorm(full.data$mpg,mu_func(full.data$disp,params['a'],params['b']),params['sigma'],log=TRUE)) 
l
exp(l)   # we can convert back to likelihood if we want...


# Example likelihood function!  --------------------------------

# Arguments:
#   params: bundled vector of free parameters for the known data-generating model
#   df: a data frame that holds the observed response variable and covariates
#   yvar: the name of the response variable (ancillary)
#   xvar: the name of the predictor variable (ancillary)

mtcars_LL <- function(params,df=mtcars,yvar="mpg",xvar="disp"){
  sum(dnorm(df$mpg,mu_func(df$disp,params['a'],params['b']),params['sigma'],log=TRUE)) 
}
mtcars_LL(unlist(params),df=mtcars,yvar="mpg",xvar="disp")


# Use numerical optimization methods to identify the maximum likelihood estimate (and the likelihood at the MLE)

optimizedLik <- optim(fn=mtcars_LL,par=params,control=list(fnscale=-1),hessian = T)  # note, the control param is set so that "optim" maximizes rather than minimizes the Log-likelihood. 


MLE = optimizedLik$par   # maximum likelihood parameter estimates
MLE

LogLik = optimizedLik$value   # log likelihood for the best model
LogLik


# visualize goodness-of-fit for the best model  ----------------------

xvals <- mtcars$disp
yvals <- mtcars$mpg
VisualizeModelWithData(xvals,MLE)


# Estimating parameter uncertainty -------------------------

# Visualize a "slice" of the likelihood function

allvals_b <- seq(-1/1000,-1/200,length=200)
paramslist <- lapply(allvals_b, function(t){MLE['b']=t;MLE } )
slice_b <- sapply(paramslist, function(t) exp(mtcars_LL(t))) 

plot(allvals_b,slice_b,type="l",main="Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Likelihood")


# Work with log-likelihood instead...

slice_b <- sapply(paramslist, function(t) mtcars_LL(t))    # 

plot(allvals_b,slice_b,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")


# zoom in closer to the MLE

allvals_b <- seq(-1/550,-1/350,length=200)
paramslist <- lapply(allvals_b, function(t){MLE['b']=t;MLE } )
slice_b <- sapply(paramslist, function(t) mtcars_LL(t)) 

plot(allvals_b,slice_b,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")

# what parameter values are within 2 log likelihood units of the best value?  -------------

plot(allvals_b,slice_b,type="l",main="Log Likelihood Slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")
abline(v=MLE['b'],lwd=3,col="blue")
abline(h=(LogLik-2))


# Generate an approximate 95% confidence interval for the "b" parameter -----------------

reasonable_b_slice <- allvals_b[slice_b>=(LogLik-2)]
reasonable_b_limits <- c(min(reasonable_b_slice),max(reasonable_b_slice))

plot(allvals_b,slice_b,type="l",main="Log Likelihood slice",xlab="Parameter Slice for \'b\'",ylab="Log-Likelihood")
abline(v=MLE['b'],lwd=3,col="blue")
abline(h=(LogLik-2))
abline(v=reasonable_b_limits,lwd=1,col="blue")


# A better confidence interval, using the likelihood "profile" -----------------

# first, visualize the likelihood surface in 2 dimensions

a <- seq(5,50,length=200)
b <- seq(-1/300,-1/800,length=200)

LLsurface <- expand.grid(a,b)
colnames(LLsurface) <- c("a","b")
paramslist <- lapply(1:nrow(LLsurface), function(t){MLE['a']=LLsurface$a[t];MLE['b']=LLsurface$b[t];MLE } )
LLsurface$LL <- sapply(paramslist, function(t) mtcars_LL(t)) 

summary(LLsurface)

ggplot(LLsurface,mapping =aes(x=a,y=b,z=LL)) +
  geom_raster(aes(fill=LL)) +
  geom_contour(breaks=seq(-100,-75,3) ,lwd=1.2) +
  scale_fill_gradient(limits=c(-125,-75)) 
  


# add a contour line, assuming deviances follow a chi-squared distribution

conf95 <- qchisq(0.95,2)/2  # this evaluates to around 3. Since we are varying freely across 2 dimensions, we use chisq with 2 degrees of freedom

ggplot(LLsurface,mapping =aes(x=a,y=b,z=LL)) +
  geom_raster(aes(fill=LL)) +
  geom_contour(breaks=seq(-100,-75,3) ,lwd=1.2) +
  geom_contour(breaks=LogLik-conf95 ,lwd=2,col="black") +
  scale_fill_gradient(limits=c(-125,-75)) +
  labs(title = "Log Likelihood Surface with 2D Confidence Region")


# visualize likelihood profiles!

profile_a <- LLsurface |> group_by(a) |> summarize(LL=max(LL)) 
profile_b <- LLsurface |> group_by(b) |> summarize(LL=max(LL)) 

reasonable_a <- profile_a$a[profile_a$LL >=(LogLik-qchisq(0.95,1)/2)]

ggplot(profile_a,aes(a,LL)) + geom_path(lwd=2) + 
  geom_vline(xintercept = c(min(reasonable_a),max(reasonable_a) )) +
  xlim(c(25,40)) + ylim(c(-110,-75))



# profile for the b parameter... 

reasonable_b <- profile_b$b[profile_b$LL >=(LogLik-qchisq(0.95,1)/2)]

ggplot(profile_b,aes(b,LL)) + geom_path(lwd=2) + 
  geom_vline(xintercept = c(min(reasonable_b),max(reasonable_b) )) 


# Compare profile and slice intervals

compdf <- rbind(profile_b,data.frame(b=allvals_b,LL=slice_b))
compdf$Method <- rep(c("Profile","Slice"),each=nrow(profile_b))

ggplot(compdf,aes(x=b,y=LL,color=Method)) + geom_path(lwd=2) + 
  geom_vline(xintercept = c(min(reasonable_b),max(reasonable_b) ),color="darkgreen",lty=2) +
  geom_vline(xintercept = reasonable_b_limits ,color="purple",lty=2) +
  scale_color_manual(values=c("darkgreen","purple"))


## use the normal approximation to estimate confidence intervals!

varcov <- solve(-optimizedLik$hessian)   # compute the variance covariance matrix for the coefficients from the hessian matrix

# approximate confidence interval as 2 standard errors from the MLE
lb = MLE - 2*sqrt(diag(varcov))
ub = MLE + 2*sqrt(diag(varcov))

cbind(MLE,lb,ub)

# compare with profile likelihood method for b
c(min(reasonable_b),max(reasonable_b) )    # close enough!?


###### Practice exercise: develop a likelihood function for estimating the probability of detection of a rare frog species

ncaps <- c(3,2,6)   # number of times out of 10 that a frog is detected at 3 known-occupided wetland sites

#### construct a likelihood function


#### find the MLE using 'optim()' in R


#### find the approximate 95% confidence interval using the "rule of 2"




# demo: likelihood ratio test  -------------------------------

curve(dchisq(x,2),0,10,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=2")


curve(dchisq(x,2),0,10,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=2")
abline(v=qchisq(0.95,2),col="red",lwd=2)


curve(dchisq(x,1),0,5,ylab="probability density",xlab="x", main="Chi-Squared distribution, df=1")
abline(v=qchisq(0.95,1),col="red",lwd=2)

