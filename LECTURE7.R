
############################################################
####                                                    ####  
####  NRES 746, Lecture 7                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Bayesian analysis #2: MCMC                        ####
############################################################



#################
# Simple example of MCMC sampling
#################


#########
# first, let's build a function that generates random numbers from a bivariate normal distribution

rbvn<-function (n, rho)   #function for drawing an arbitrary number of independent samples from the bivariate standard normal distribution. 
{
        x <- rnorm(n, 0, 1)
        y <- rnorm(n, rho * x, sqrt(1 - rho^2))
        cbind(x, y)
}

#########
# Now, plot the random draws from this distribution, make sure this makes sense!

bvn<-rbvn(10000,0.98)
par(mfrow=c(3,2))
plot(bvn,col=1:10000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))



###############
# Metropolis-Hastings implementation of bivariate normal sampler... 

library(mvtnorm)    # load a package that allows us to compute probability densities for mv normal distribution 

metropolisHastings <- function (n, rho=0.98){    # a MCMC sampler implementation of a bivariate random number generator
    mat <- matrix(ncol = 2, nrow = n)   # matrix for storing the random samples
    x <- 0   # initial values for all parameters
    y <- 0
    prev <- dmvnorm(c(x,y),mean=c(0,0),sigma = matrix(c(1,rho,rho,1),ncol=2))   # probability density of the distribution at the starting values
    mat[1, ] <- c(x, y)        # initialize the markov chain
    counter <- 1
    while(counter<=n) {
      newx <- rnorm(1,x,0.5)     # make a jump. Note the symmetrical proposal distribution
      newy <- rnorm(1,y,0.5)
      
      newprob <- dmvnorm(c(newx,newy),sigma = matrix(c(1,rho,rho,1),ncol=2))    # assess whether the new jump is good!
      ratio <- newprob/prev   # compute the ratio of probabilities at the old (jump from) and proposed (jump to) locations. 
      
      prob.accept <- min(1,ratio)     # decide the probability of accepting the new jump!
      rand <- runif(1)
      if(rand<=prob.accept){
        x=newx;y=newy    # set x and y to the new location
        mat[counter,] <- c(x,y)    # store this in the storage array 
        counter=counter+1
        prev <- newprob    # get ready for the next iteration
      }
      
    }
    return(mat)
}


###########
# Test the new M-H sampler

bvn<-metropolisHastings(10000,0.98)
par(mfrow=c(3,2))
plot(bvn,col=1:10000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))


############
# MCMC implementation of the Myxomatosis example from the Bolker book
############

library(emdbook)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)
head(Myx)


###########
# Visualize the Myxomatosis data for the 100th time!

hist(Myx$titer,freq=FALSE)


#########
# ... and overlay a proposed data-generating model (gamma distribution)

hist(Myx$titer,freq=FALSE)
curve(dgamma(x,shape=40,scale=0.15),add=T,col="red")


##############
# define 2-D parameter space!
##############

shapevec <- seq(3,100,by=0.1)   
scalevec <- seq(0.01,0.5,by=0.001)

##############
# define the likelihood surface across this grid within parameter space
##############

GammaLogLikelihoodFunction <- function(params){
  sum(dgamma(Myx$titer,shape=params['shape'],scale=params['scale'],log=T))
}
surface2D <- matrix(nrow=length(shapevec),ncol=length(scalevec))   # initialize storage variable

newparams <- c(shape=50,scale=0.2)
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  for(j in 1:length(scalevec)){
    newparams['scale'] <- scalevec[j]
    surface2D[i,j] <- GammaLogLikelihoodFunction(newparams) 
  }
}

############
# Visualize the likelihood surface
############

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)


############
# Write a non-log-transformed likelihood function

GammaLikelihoodFunction <- function(params){
  prod(dgamma(Myx$titer,shape=params['shape'],scale=params['scale'],log=F))   
}

params <- c(shape=40,scale=0.15) 
params
GammaLikelihoodFunction(params)


#############
# Function for returning the prior probability density for any point in parameter space  

GammaPriorFunction <- function(params){
  prior <- c(shape=NA,scale=NA)
  prior['shape'] <- dgamma(params['shape'],shape=0.01,scale=100)
  prior['scale'] <- dgamma(params['scale'],shape=0.001,scale=1000)
  # prior['shape'] <- dunif(params['shape'],3,100)        # alternative: could use uniform prior!
  # prior['scale'] <- dunif(params['scale'],0.01,0.5)
  return(prod(prior))
}

curve(dgamma(x,shape=0.01,scale=1000),3,100)

params <- c(shape=40,scale=0.15) 
params
GammaPriorFunction(params)



############
# Function for computing the ratio of posterior densities between any two points in parameter space

PosteriorRatio <- function(oldguess,newguess){
  oldLik <- max(1e-90,GammaLikelihoodFunction(oldguess))   # compute likelihood and prior density at old guess
  oldPrior <- max(1e-90,GammaPriorFunction(oldguess))
  newLik <- GammaLikelihoodFunction(newguess)             # compute likelihood and prior density at new guess
  newPrior <- GammaPriorFunction(newguess)
  return((newLik*newPrior)/(oldLik*oldPrior))          # compute ratio of weighted likelihoods
}

oldguess <- params
newguess <- c(shape=39,scale=0.15)

PosteriorRatio(oldguess,newguess)


############
# Define proposal distribution for jumps in parameter space (use normal distribution)!

     # function for making new guesses
newGuess <- function(oldguess){
  sdshapejump <- 4
  sdscalejump <- 0.07
  jump <- c(shape=rnorm(1,mean=0,sd=sdshapejump),scale=rnorm(1,0,sdscalejump))
  newguess <- abs(oldguess + jump)
  return(newguess)
}
  # set a new "guess" near to the original guess

newGuess(oldguess=params)   
newGuess(oldguess=params)
newGuess(oldguess=params)


##########
# Set a starting point in parameter spacer

startingvals <- c(shape=75,scale=0.28)    # starting point for the algorithm


###########
# Try our new functions

newguess <- newGuess(startingvals)    # take a jump in parameter space
newguess

PosteriorRatio(startingvals,newguess)   # difference in posterior ratio


###############
# Visualize the Metropolis-Hastings routine:

chain.length <- 10
oldguess <- startingvals
guesses <- matrix(0,nrow=chain.length,ncol=2)
colnames(guesses) <- names(startingvals)

counter <- 1
while(counter <= chain.length){
  newguess <- newGuess(oldguess)
  post.rat <- PosteriorRatio(oldguess,newguess)
  prob.accept <- min(1,post.rat)
  rand <- runif(1)
  if(rand<=prob.accept){
    oldguess <- newguess
    guesses[counter,] <- newguess 
    counter=counter+1
  }
}

# visualize!

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")


##########
# Get more MCMC samples

chain.length <- 100
oldguess <- startingvals
guesses <- matrix(0,nrow=chain.length,ncol=2)
colnames(guesses) <- names(startingvals)

counter <- 1
while(counter <= chain.length){
  newguess <- newGuess(oldguess)
  post.rat <- PosteriorRatio(oldguess,newguess)
  prob.accept <- min(1,post.rat)
  rand <- runif(1)
  if(rand<=prob.accept){
    oldguess <- newguess
    guesses[counter,] <- newguess 
    counter=counter+1
  }
}

# visualize!

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")


############
# And more...

chain.length <- 1000
oldguess <- startingvals
guesses <- matrix(0,nrow=chain.length,ncol=2)
colnames(guesses) <- names(startingvals)

counter <- 1
while(counter <= chain.length){
  newguess <- newGuess(oldguess)
  post.rat <- PosteriorRatio(oldguess,newguess)
  prob.accept <- min(1,post.rat)
  rand <- runif(1)
  if(rand<=prob.accept){
    oldguess <- newguess
    guesses[counter,] <- newguess 
    counter=counter+1
  }
}

# visualize!

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")


#############
# Evaluate "traceplot" for the MCMC samples...

##### Shape parameter

plot(1:chain.length,guesses[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")


###### Scale parameter

plot(1:chain.length,guesses[,'scale'],type="l",main="scale parameter",xlab="iteration",ylab="scale")


############
# Remove "burn-in" (allow MCMC routine some time to get to the posterior)

burn.in <- 100
MCMCsamples <- guesses[-c(1:burn.in),]

chain.length=chain.length-burn.in
plot(1:chain.length,MCMCsamples[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")
plot(1:chain.length,MCMCsamples[,'scale'],type="l",main="scale parameter",xlab="iteration",ylab="scale")


##########
# Try again- run for much longer

chain.length <- 20000
oldguess <- startingvals
guesses <- matrix(0,nrow=chain.length,ncol=2)
colnames(guesses) <- names(startingvals)

counter <- 1
while(counter <= chain.length){
  newguess <- newGuess(oldguess)
  post.rat <- PosteriorRatio(oldguess,newguess)
  prob.accept <- min(1,post.rat)
  rand <- runif(1)
  if(rand<=prob.accept){
    oldguess <- newguess
    guesses[counter,] <- newguess 
    counter=counter+1
  }
}

# visualize!

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")


#############
# Use longer "burn-in"

burn.in <- 5000
MCMCsamples <- guesses[-c(1:burn.in),]
chain.length=chain.length-burn.in


plot(1:chain.length,MCMCsamples[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")
plot(1:chain.length,MCMCsamples[,'scale'],type="l",main="scale parameter",xlab="iteration",ylab="scale")


##########
# "thin" the MCMC samples

thinnedMCMC <- MCMCsamples[seq(1,chain.length,by=10),]
plot(1:nrow(thinnedMCMC),thinnedMCMC[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")
plot(1:nrow(thinnedMCMC),thinnedMCMC[,'scale'],type="l",main="scale parameter",xlab="iteration",ylab="scale")


# Visualize the posterior!

plot(density(thinnedMCMC[,'scale']),main="scale parameter",xlab="scale")
plot(density(thinnedMCMC[,'shape']),main="shape parameter",xlab="shape")


#########
# More visual posterior checks...

par(mfrow=c(3,2))
plot(thinnedMCMC,col=1:10000)
plot(thinnedMCMC,type="l")
plot(ts(thinnedMCMC[,1]))
plot(ts(thinnedMCMC[,2]))
hist(thinnedMCMC[,1],40)
hist(thinnedMCMC[,2],40)
par(mfrow=c(1,1))


#############
# Simple example of a Gibbs sampler
#############

########
# first, recall our simple bivariate normal sampler

rbvn<-function (n, rho){  #function for drawing an arbitrary number of independent samples from the bivariate standard normal distribution. 
        x <- rnorm(n, 0, 1)
        y <- rnorm(n, rho * x, sqrt(1 - rho^2))
        cbind(x, y)
}

bvn<-rbvn(10000,0.98)
par(mfrow=c(3,2))
plot(bvn,col=1:10000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))


#############
# Now construct a Gibbs sampler alternative

gibbs<-function (n, rho){    # a gibbs sampler implementation of a bivariate random number generator
    mat <- matrix(ncol = 2, nrow = n)   # matrix for storing the random samples
    x <- 0
    y <- 0
    mat[1, ] <- c(x, y)        # initialize the markov chain
    for (i in 2:n) {
            x <- rnorm(1, rho * y, sqrt(1 - rho^2))        # sample from x conditional on y
            y <- rnorm(1, rho * x, sqrt(1 - rho^2))        # sample from y conditional on x
            mat[i, ] <- c(x, y)
    }
    mat
}


##########
# Test the Gibbs sampler

bvn<-gibbs(10000,0.98)
par(mfrow=c(3,2))
plot(bvn,col=1:10000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))


###########
# Myxomatosis example in BUGS modeling language

##########
# Write the BUGS model to file

cat("
  model {
    
    #############
    # LIKELIHOOD
    ############
    for(obs in 1:n.observations){
      titer[obs] ~ dgamma(shape,rate)
    }
    
    #############
    # PRIORS
    ############
    shape ~ dgamma(0.001,0.001)
    scale ~ dgamma(0.01,0.01)
    rate <- 1/scale
  }
", file="BUGSmodel.txt")



############
# Encapsulate the data into a single "list" object

myx.data.for.bugs <- list(
  titer = Myx$titer,
  n.observations = length(Myx$titer)
)

myx.data.for.bugs


###########
# Function for generating random initial values for all free parameters

init.vals.for.bugs <- function(){
  init.list <- list(
    shape=runif(1,20,100),
    scale=runif(1,0.05,0.3)
  )
  return(init.list)
}

init.vals.for.bugs()
init.vals.for.bugs()
init.vals.for.bugs()


###########
# Run JAGS!!!!
##########

library(R2jags)

library(coda)

params.to.store <- c("shape","scale")

jags.fit <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=5000,model.file="BUGSmodel.txt",n.chains = 3,n.burnin = 0 )

jagsfit.mcmc <- as.mcmc(jags.fit)   # convert to "MCMC" object (coda package)

summary(jagsfit.mcmc)

plot(jagsfit.mcmc)


################
# Run the chains for longer!

jags.fit <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=50000,model.file="BUGSmodel.txt",n.chains = 3, n.burnin=10000,n.thin = 20)

jagsfit.mcmc <- as.mcmc(jags.fit)   # convert to "MCMC" object (coda package)

summary(jagsfit.mcmc)

plot(jagsfit.mcmc)


##############
# Run convergence diagnostics

gelman.diag(jagsfit.mcmc)

