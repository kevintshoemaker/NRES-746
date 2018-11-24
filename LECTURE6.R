
############################################################
####                                                    ####  
####  NRES 746, Lecture 6                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Bayesian analysis #1: concepts                    ####
############################################################



###################
### Bayesian analysis example using Binomial distribution


######
# first visualize the "prior" for the probability *p* as a uniform distribution

curve(dunif(x),ylim=c(0,2),col="red", ylab = "probability", xlab="parameter \"p\"")

#hist(runif(10000),freq=F,ylim=c(0,2),col="red")


#############
# Alternative prior: beta distribution (conjugate prior)

curve(dbeta(x,1,1),ylim=c(0,2),col="red",ylab = "probability", xlab="parameter \"p\"")

#hist(rbeta(10000,1,1),freq=F,ylim=c(0,2),col="red")   # histogram of random numbers from a flat beta distribution


#############
# frog call example: imagine we detected the frog in 3 of 10 visits to a known-occupied wetland

######
# visualize the data likelihood alongside the prior probability

#   # recall that likelhood surface is not a probability distribution (thus, the 2 y axes)

data = 3
param.space <- seq(0,1,by=0.001)
likelihood <- dbinom(data,size=10,prob=param.space)
par(mai=c(1,1,0,1))
curve(dbeta(x,1,1),ylim=c(0,2),col="blue",lty=2,ylab="Probability density (prior)",xlab="param.space")
points(param.space,likelihood*5,type="l",col="red",lwd=2,lty=2)
axis(4,at=seq(0,2,by=0.4),labels = seq(0,0.5,by=.1))
mtext("Likelihood (relative plausibility)", side=4, col="red",line=3)



###########
# Brute-force Bayes

####
# prior across parameter space 

prior <- dbeta(param.space,shape1=1,shape2=1)    # flat prior
#prior

######
## Numerator for Bayes rule: weight the data likelihood by the prior

weighted.likelihood <- likelihood*prior      # Numerator for Bayes rule

######
## Denominator for Bayes rule: compute normalization constant

normalization.constant <- sum(weighted.likelihood)   # "evidence" is a single constant number: the sum of the weighted likelihoods!

######
## Posterior!! (numerator/denominator)

posterior <- weighted.likelihood/normalization.constant   # this is Bayes' rule!

#######
## Plot it out!

par(mai=c(1,1,0,1))
plot(param.space,prior,ylim=c(0,5),type="l",lwd=1,lty=2,col="blue",ylab="Probability Density",xlab="param.space")
points(param.space,posterior*length(param.space),type="l",col="blue",lwd=2,lty=1)  # convert posterior to probability density
points(param.space,likelihood*5,type="l",col="red",lwd=1,lty=2)
axis(4,at=seq(0,2,by=0.4),labels = seq(0,0.5,by=.1))
mtext("Likelihood", side=4, col="red",line=3)


#############
# Try an informative prior!

prior <- dbeta(param.space,shape1=15,shape2=5)
#prior

## weight the data likelihood by the prior

weighted.likelihood <- likelihood*prior

## compute normalization constant

normalization.constant <- sum(weighted.likelihood)

## Posterior!!

posterior <- weighted.likelihood/normalization.constant

## Plot it out!
par(mai=c(1,1,0,1))
plot(param.space,prior,ylim=c(0,5),type="l",lwd=1,lty=2,col="blue",ylab="Probability Density",xlab="param.space")
points(param.space,posterior*length(param.space),type="l",col="blue",lwd=2,lty=1)
points(param.space,likelihood*5,type="l",col="red",lwd=1,lty=2)
axis(4,at=seq(0,2,by=0.4),labels = seq(0,0.5,by=.1))
mtext("Likelihood", side=4, col="red",line=3)

############
# Collect more data and try again...

moredata <- c(3, 1, 6, 2, 3, 2, 6, 1, 3, 3)

## prior
prior <- dbeta(param.space,shape1=15,shape2=5)

## likelihood
likelihood <- sapply(param.space,function(t) prod(dbinom(moredata,size=10,prob=t)))

## weight the data likelihood by the prior
weighted.likelihood <- likelihood*prior

## compute normalization constant

normalization.constant <- sum(weighted.likelihood)

## Posterior!!

posterior <- weighted.likelihood/normalization.constant

## Plot it out!
par(mai=c(1,1,0,1))
plot(param.space,prior,ylim=c(0,10),type="l",lwd=1,lty=2,col="blue",ylab="Probability Density",xlab="param.space")
points(param.space,posterior*length(param.space),type="l",col="blue",lwd=2,lty=1)
points(param.space,likelihood*1e9,type="l",col="red",lwd=1,lty=2)
axis(4,at=seq(0,6,by=1),labels = seq(0,6e-9,by=1e-9))
mtext("Likelihood", side=4, col="red",line=3)


#######
# Try a very informative prior!

likelihood <- dbinom(data,size=10,prob=param.space)

prior <- dbeta(param.space,shape1=150,shape2=50)
#prior

## weight the data likelihood by the prior

weighted.likelihood <- likelihood*prior

## compute normalization constant

normalization.constant <- sum(weighted.likelihood)

## Posterior!!

posterior <- weighted.likelihood/normalization.constant

## Plot it out!
par(mai=c(1,1,0,1))
plot(param.space,prior,ylim=c(0,15),type="l",lwd=1,lty=2,col="blue",ylab="Probability Density",xlab="param.space")
points(param.space,posterior*length(param.space),type="l",col="blue",lwd=2,lty=1)
points(param.space,likelihood*5,type="l",col="red",lwd=1,lty=2)
axis(4,at=seq(0,2,by=0.4),labels = seq(0,0.5,by=.1))
mtext("Likelihood", side=4, col="red",line=3)

###############
# Do it again- this time with conjugate priors...

### PRIOR

prior_beta <- c(shape1=1,shape2=1)
curve(dbeta(x,prior_beta['shape1'],prior_beta['shape2']),ylim=c(0,5),ylab="Prob Density",col="blue",lwd=1,lty=2,xlab="param.space")

### POSTERIOR

curve(dbeta(x,prior_beta['shape1']+data,prior_beta['shape2']+(10-data)),ylim=c(0,4),ylab="Prob Density",col="blue",lwd=2,lty=1,xlab="param.space",add=T)


#########
# With informative prior...

### PRIOR
prior_beta <- c(shape1=15,shape2=5)
curve(dbeta(x,prior_beta['shape1'],prior_beta['shape2']),ylim=c(0,5),ylab="Prob Density",col="blue",lwd=1,lty=2,xlab="param.space")


### POSTERIOR

curve(dbeta(x,prior_beta['shape1']+data,prior_beta['shape2']+(10-data)),ylim=c(0,4),ylab="Prob Density",col="blue",lwd=2,lty=1,xlab="param.space",add=T)


graphics.off()


########
# And with super informative prior...

### PRIOR
prior_beta <- c(shape1=150,shape2=50)
curve(dbeta(x,prior_beta['shape1'],prior_beta['shape2']),ylim=c(0,15),ylab="Prob Density",col="blue",lwd=1,lty=2,xlab="param.space")


### POSTERIOR
curve(dbeta(x,prior_beta['shape1']+data,prior_beta['shape2']+(10-data)),ylim=c(0,15),ylab="Prob Density",col="blue",lwd=2,lty=1,xlab="param.space",add=T)


##########
# Example: bayesian point estimates can differ markedly from MLE

curve(dlnorm(x,4,1),from=0.001,to=200,ylab="prob density")  # use a lognormal distribution for example of skewed posterior dist...


##########
# Compute and plot the mean and the mode of the distribution

param.space2 <- seq(0.001,200,length=10000)
skewed.posterior <- dlnorm(param.space2,4,1)
mean <- mean(rlnorm(10000,4,1))
mode <- param.space2[which.max(skewed.posterior)]
plot(param.space2,skewed.posterior,type="l",ylab="prob density")
abline(v=c(mean,mode),col=gray(0.5),lwd=3,lty=2)   # add to plot


############
# Do the same for the frog detection example from above

graphics.off()
### POSTERIOR
posterior <- dbeta(param.space,1+data,1+(10-data))
mean <- mean(rbeta(10000,1+data,1+(10-data)))
mode <- param.space[which.max(posterior)]
plot(param.space,posterior,type="l",col="blue",lwd=2)
abline(v=c(mean,mode),col=gray(0.5),lwd=3,lty=2)   # add to plot


#############
# A Bayesian confidence interval (95% credible interval)

### POSTERIOR

curve(dbeta(x,1+data,1+(10-data)),ylim=c(0,4),ylab="Prob Density",col="blue",lwd=2,lty=1,xlab="param.space")

### CREDIBLE INTERVAL

credible.interval <- qbeta(c(0.025,0.975),1+data,1+(10-data))     # get the credible interval using the quantile method

abline(v=credible.interval,col=gray(0.5),lwd=3,lty=2)   # add to plot


###############
# Bayesian analysis without a conjugate prior


###########
# Revisit the Myxomatosis example

library(emdbook)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)
head(Myx)


hist(Myx$titer,freq=FALSE)    # visualize the data, again!


###### Error is modeled as gamma distributed

hist(Myx$titer,freq=FALSE)
curve(dgamma(x,shape=40,scale=0.15),add=T,col="red")


##########
# recall our likelihood function for these data (not on log scale this time!)

GammaLikelihoodFunction <- function(params){
  prod(dgamma(Myx$titer,shape=params['shape'],scale=params['scale']))
}

params <- c(40,0.15) 
names(params) <- c("shape","scale")
params
GammaLikelihoodFunction(params)


##############
# define 2-D parameter space (in real probability scale)!
##############

shapevec <- seq(10,100,by=0.1)   
scalevec <- seq(0.01,0.3,by=0.001)

##############
# define the likelihood surface across this grid within parameter space
##############

likelihood2D <- matrix(nrow=length(shapevec),ncol=length(scalevec))   # initialize storage variable

newparams <- params
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  for(j in 1:length(scalevec)){
    newparams['scale'] <- scalevec[j]
    likelihood2D[i,j] <- GammaLikelihoodFunction(newparams) 
  }
}

############
# Visualize the likelihood surface
############

image(x=shapevec,y=scalevec,z=likelihood2D,zlim=c(1e-70,1e-17),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=likelihood2D,levels=c(1e-18,1e-17),add=T)


#############
# compute the area of each pixel (for probability density computation)

pixelArea <- 0.0001  # for determining probability densities

##############
# define the prior probability surface across this grid within parameter space
##############

prior2D <- matrix(1, nrow=length(shapevec),ncol=length(scalevec))   # initialize prior
prior2D <- prior2D/length(prior2D)

############
# Visualize the 2-D prior distribution
############

image(x=shapevec,y=scalevec,z=prior2D,zlim=c(0,0.001),col=rainbow(10))


###########
# Apply Bayes Rule!

weighted.likelihood <- prior2D * likelihood2D    # numerator of Bayes rule
normalization.constant <- sum(weighted.likelihood)    # denominator of Bayes rule

posterior2D <- weighted.likelihood/normalization.constant

############
# Visualize the 2-D posterior distribution
############

image(x=shapevec,y=scalevec,z=(posterior2D/pixelArea),zlim=c(0,5),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=(posterior2D/pixelArea),levels=c(1:4),add=T,drawlabels=FALSE)


###########
# try to find the contour that contains 95% of our degree of belief!

possible.contours <- data.frame(contour = seq(0.13e-4,1e-4,length=100), quantile = NA)
i=1
for(i in 1:nrow(possible.contours)){
  ndx <- which(posterior2D<possible.contours$contour[i],arr.ind = T)
  possible.contours$quantile[i] <- sum(posterior2D[ndx])
}

head(possible.contours,10)


###########
# Visualize the 2D credible region (HPD)

q95 <- 1.739394e-05
image(x=shapevec,y=scalevec,z=posterior2D,zlim=c(0.5e-11,5e-4),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=posterior2D,levels=q95,add=T,lwd=3,col="red",drawlabels=FALSE)

#############
# Visualize a point estimate

image(x=shapevec,y=scalevec,z=posterior2D,zlim=c(0.5e-11,5e-4),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=posterior2D,levels=q95,add=T,lwd=3,col="red",drawlabels=FALSE)

meanshape <- sum(shapevec*posterior2D) 
meanscale <- sum(scalevec*posterior2D)
points(meanshape,meanscale,pch=20,cex=2,col="red")  # posterior mean

mode <- which(posterior2D==max(posterior2D),arr.ind = T)
points(shapevec[mode[1]],scalevec[mode[2]], pch="X",col="black",cex=1.5)  # posterior mode


##########
# Plot out posterior distributions separately for each parameter

marginal.dist.shape <- apply(posterior2D,1,mean)     # factor out the scale param- focus only on the shape param.
plot(shapevec,(marginal.dist.shape/sum(marginal.dist.shape))/0.1,type="l",lwd=2,col="blue",ylab="probability density",main="Posterior probability")
abline(v=meanshape)

marginal.dist.scale <- apply(posterior2D,2,mean)
plot(scalevec,(marginal.dist.scale/sum(marginal.dist.scale))/0.001,type="l",lwd=2,col="blue",ylab="probability density",main="Posterior probability")
abline(v=meanscale)

meanshape
meanscale


################
# Sample parameters from the joint posterior

SampleFromPosterior <- function(n){
  shape <- rep(shapevec,times=length(scalevec))
  scale <- rep(scalevec,each=length(shapevec))
  jointparams <- data.frame(shape=shape,scale=scale)
  probs <- as.vector(posterior2D)
  samples <- sample(c(1:length(probs)),size=n,replace=TRUE,prob=probs)
  jointparams[samples,]
}

samples<-SampleFromPosterior(n=10000)
par(mfrow=c(3,2))
plot(samples,col=1:10000)
plot(samples,type="l")
plot(ts(samples[,1]))
plot(ts(samples[,2]))
hist(samples[,1],40)
hist(samples[,2],40)
par(mfrow=c(1,1))

