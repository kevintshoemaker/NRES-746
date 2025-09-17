

#  NRES 746, Lecture 7   
#   Bayesian analysis #2: MCMC  ---------------------------------                   


# Simple example of MCMC sampling -----------------------

# first, let's build a function that generates random numbers from a bivariate standard normal distribution

rbvn<-function (n, rho){   #function for drawing an arbitrary number of independent samples from the bivariate standard normal distribution. 
        x <- rnorm(n, 0, 1)
        y <- rnorm(n, rho * x, sqrt(1 - rho^2))
        cbind(x, y)
}

# Now, plot the random draws from this distribution, make sure this makes sense!

bvn_true<-rbvn(10000,0.9)
par(mfrow=c(2,2))
plot(ts(bvn_true[,1]))
plot(ts(bvn_true[,2]))
hist(bvn_true[,1],40)
hist(bvn_true[,2],40)
par(mfrow=c(1,1))



plot(bvn_true[,1],bvn_true[,2])

vis_mcmc2 <- function(mat){   # assume mat has 1 row per sample, and 2 columns (# params)
  par(mfrow=c(3,ncol(mat)))
  plot(mat,col=1:nrow(mat)) ; plot(mat,type="l")
  plot(ts(mat[,1])) ; plot(ts(mat[,2]))
  hist(mat[,1],40) ; hist(mat[,2],40)
  par(mfrow=c(1,1))
}


## Metropolis-Hastings implementation of bivariate normal sampler...   ----------

proposal <- function(x) rnorm(1,x,0.4)
cond_prob <- function(this,xval,rho) dnorm(this,rho*xval,1-rho^2) 
choose_mh <- function(prev,new,other,rho){
  r = min(cond_prob(new,other,rho)/cond_prob(prev,other,rho),1)
  ifelse(r>runif(1),new,prev)
}

metropolisHastings <- function (n, rho){    # an MCMC bivariate random number generator
    mat <- matrix(NA,n,2)   # matrix for storing the random samples
    prev = c(0,0); mat[1,] <- prev       # initial values for all parameters
    for(i in 1:n){
      for(j in 1:2) prev[j] <- choose_mh(prev[j],proposal(prev[j]),prev[setdiff(1:2,j)] ,rho)
      mat[i,] <- prev
    }
    mat
}


# Test our M-H sampler

bvn_mh<-metropolisHastings(5000,rho=0.9)
vis_mcmc2(bvn_mh)


# Simple example of a Gibbs sampler ----------------

# first, recall our 'true' bivariate normal sampler

vis_mcmc2(bvn_true)


## Now construct a Gibbs sampler  ---------------

cond_prob2 <- function(xval,rho) rnorm(1, rho * xval, sqrt(1 - rho^2))   # sample random value from full conditional

gibbs<-function (n, rho){    # a Gibbs sampler for bivariate normal
    mat <- matrix(ncol = 2, nrow = n)   # matrix for storing the random samples
    prev <- c(0,0); mat[1, ] <- prev     # initialize the markov chain
    for (i in 2:n) {
      prev[1] <- cond_prob2(prev[2],rho)   # sample from full conditional
      prev[2] <- cond_prob2(prev[1],rho)
      mat[i,] <- prev
    }
    mat
}


# Test the Gibbs sampler ------------------

bvn_gbs <-gibbs(10000,0.9)
vis_mcmc2(bvn_gbs)


rho = 0.9  # set correlation as a global constant

  # compute the unnormalized log posterior
log_posterior <- function(x) dnorm(x[1],log=T) + dnorm(x[2],rho*x[1],1-rho^2,log=T)


library(pracma)   # load package capable of computing gradient numerically at any point in parameter space

   # compute the gradient of the log posterior density function at a point x
gradient <- function(x) pracma::grad(log_posterior,x)

# Ancillary code for HMC ----------
  # modified from https://jonnylaw.rocks/posts/2019-07-31-hmc/

leapfrog_step <- function(gradient, step_size, position, momentum, d) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum1 + gradient(position1) * 0.5 * step_size
  matrix(c(position1, momentum2), ncol = d*2)
}


leapfrogs <- function(gradient, step_size, l, position, momentum, d) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum, d)
    position <- pos_mom[seq_len(d)]   # position is the first d elements of the row vector
    momentum <- pos_mom[-seq_len(d)] # momentum is the final d elements of the row vector
  }
  pos_mom
}

log_acceptance <- function(propPosition,
                           propMomentum,
                           position,
                           momentum,
                           log_posterior) {
  log_posterior(propPosition) + sum(dnorm(propMomentum, log = T)) - 
    log_posterior(position) - sum(dnorm(momentum, log = T))
}
hmc_step <- function(log_posterior, gradient, step_size, l, position) {
  d <- length(position)
  momentum <- rnorm(d)    # initial momentum- a kick to get the sampler going
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum, d)   # position and momentum vectors  
  propPosition <- pos_mom[seq_len(d)]    # separate into position and momentum vectors
  propMomentum <- pos_mom[-seq_len(d)]
  a <- log_acceptance(propPosition, propMomentum, position, momentum, log_posterior)
  if (log(runif(1)) < a) {    # this is a metropolis procedure! 
    propPosition
  } else {
    position
  }
}
hmc <- function(log_posterior, gradient, step_size, l, initP, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(initP))   # initialize the output matrix
  out[1, ] <- initP   # get the sampler started
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i-1,]) # one HMC step
  }
  out  # fully filled-in matrix of steps in parameter space
}
# Test the HMC sampler ------------------

bvn_hmc <-hmc(log_posterior, gradient, .2, 10, c(0,0), 1000)  
vis_mcmc2(bvn_hmc)

# Using MCMC to fit the Myxomatosis example from the Bolker book --------------

library(emdbook)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)
head(Myx)


# Visualize the Myxomatosis data for the 100th time!

hist(Myx$titer,freq=FALSE)


# define 2-D parameter space!

shapevec <- seq(1,150,length=200)        # divide parameter space into tiny increments
ratevec <- seq(0.5,30,length=200)

# define the likelihood surface  -------------

loglik <- function(pars){
  sum(dgamma(Myx$titer,shape=pars['shape'],rate=pars['rate'],log = T))
}

parmsurface <- expand.grid(shapevec,ratevec); colnames(parmsurface) <- c("shape","rate")
parmsurface$ll <- sapply(1:nrow(parmsurface), function(t) loglik(unlist(parmsurface[t,1:2]))  )

library(ggplot2)
ggplot(parmsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the log likelihood surface
  geom_raster(aes(fill=ll)) +
  scale_fill_gradient(limits=c(-150,-37)) +
  geom_contour(aes(z=ll),breaks=c(-45,-40) ,lwd=1.2)


# Function for returning the log prior probability density for any 2D parameter vector 
logprior <- function(params){
  dgamma(params['shape'],0.001,0.001,log=T)  +  
  dgamma(params['rate'],0.001,0.001,log=T)
}

# curve(dgamma(x,shape=0.001,rate=0.001),3,100)   # visualize gamma
# params <- c(shape=40,rate=7)    # test function
# logprior(params)

parmsurface$pr <- sapply(1:nrow(parmsurface), function(t) logprior(unlist(parmsurface[t,1:2]))  )
ggplot(parmsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the log likelihood surface
  geom_raster(aes(fill=pr)) +
  scale_fill_gradient(limits=c(-25,-13.14)) 


# Function for computing the ratio of posterior densities -----------------

eval_jump <- function(old,new){
  oldnum <- loglik(old) + logprior(old)   # compute likelihood and prior density at old guess
  newnum <- loglik(new) + logprior(new)              # compute likelihood and prior density at new guess
  c(diff= unname(newnum-oldnum) )         # compute ratio of weighted likelihoods (log scale)
}

params <- c(shape=37,rate=6)
old <- params    # test function
new <- c(shape=39,rate=7)
eval_jump(old,new)


# Define proposal distribution --------------------------
    # use bivariate normal distribution to make a guess

jump_vcv <- matrix(c(1,0.3,0.3,0.4),ncol=2)    # in real life, we would tune these parameters to optimize our sampler

     # function for making new guesses
make_guess <- function(old){
  newguess <- mvtnorm::rmvnorm(1,old,jump_vcv)[1,] 
  ifelse(newguess<0.001,0.001,newguess)   # make sure we don't get any negative guesses. This induces an asymmetry, but we won't worry about that right now as it is unlikely to influence our MCMC samples
}
params
make_guess(params)     # set a new "guess" near to the original guess


# Set a starting point in parameter space -------------------

startingvals <- c(shape=75,rate=4)    # starting point for the algorithm


# Try our new functions  ------------------

newguess <- make_guess(startingvals)    # take a jump in parameter space
newguess

eval_jump(startingvals,newguess)   # difference in posterior ratio


# Visualize the Metropolis-Hastings routine: ---------------

mh <- function(n,st){    # function for doing M-H MCMC
  chain <- matrix(nrow=n,ncol=length(st),dimnames=list(1:n,names(st)))
  chain[1,] <- startingvals
  for(i in 2:n){
    prop <- make_guess(chain[i-1,])    # proposal jump
    prob_accept <- min(1,exp(eval_jump(chain[i-1,],prop)) )
    if(prob_accept >= runif(1)){    # choose the jump in accordance with its relative probability under the posterior
      chain[i,] = prop
    }else{
      chain[i,] <- chain[i-1,]  
    }
  }
  chain
}

chain= mh(100,startingvals)

# visualize!
parmsurface$post <- parmsurface$ll + parmsurface$pr
post_plot = ggplot(parmsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the log likelihood surface
  geom_raster(aes(fill=post)) +
  scale_fill_gradient(limits=c(-150,-57)) +
  geom_contour(aes(z=post),breaks=c(-65,-59) ,lwd=1.2)
post_plot +  geom_path(data=chain,aes(x=shape,y=rate),lwd=1.2,col="white")


# Get more MCMC samples --------------

chain= mh(1000,startingvals)
post_plot +  geom_path(data=chain,aes(x=shape,y=rate),lwd=1.2,col="white")


# And more... -------------------

chain= mh(10000,startingvals)
post_plot +  geom_path(data=chain,aes(x=shape,y=rate),lwd=1.2,col="white")


# Evaluate "traceplot" for the MCMC samples... ---------------------

## Shape parameter

plot(1:nrow(chain),chain[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")


## Rate parameter

plot(1:nrow(chain),chain[,'rate'],type="l",main="rate parameter",xlab="iteration",ylab="rate")


# Remove "burn-in" (allow MCMC routine some time to get to the posterior) --------------

chain <- chain[-c(1:1000),]    # remove first 1000 sample

plot(1:nrow(chain),chain[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")
plot(1:nrow(chain),chain[,'rate'],type="l",main="rate parameter",xlab="iteration",ylab="rate")


# Change the VCV for proposal distribution

jump_vcv <- matrix(c(1.1,0.7,0.7,.8),ncol=2) 

# Try again- run for much longer ---------------------

chain= mh(100000,startingvals)   # takes ~10 seconds to run


# Use longer "burn-in" and thin ------------------

chain <- chain[-c(1:25000),]    # remove first 25000 sample
chain <- chain[seq(1,nrow(chain),5),]   # keep every fifth sample


plot(1:nrow(chain),chain[,'shape'],type="l",main="shape parameter",xlab="iteration",ylab="shape")
plot(1:nrow(chain),chain[,'rate'],type="l",main="rate parameter",xlab="iteration",ylab="rate")


acf(chain[,"shape"],lag.max=500)


# Visualize the posterior!

plot(density(chain[,'rate']),main="rate parameter",xlab="rate")
plot(density(chain[,'shape']),main="shape parameter",xlab="shape")


# More visual posterior checks... -----------------

vis_mcmc2(chain)


data {
  int<lower=0> N;
  vector[N] titer;
}

parameters {
  real<lower=0> shape;
  real<lower=0> rate;
}

model {
  shape ~ gamma(0.001,0.001);   // prior on shape
  rate ~ gamma(0.001,0.001);   // prior on rate
  titer ~ gamma(shape, rate);   // likelihood
}


library(cmdstanr)
mod1 <- cmdstan_model("myx.stan") # Compile stan model


# Encapsulate the data into a single "list" object ------------------

stan_data <- list(
    N = nrow(Myx),
    titer = Myx$titer
)


# Run stan  ------------------

fit1 <- mod1$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 200,
  iter_sampling = 500
)

fit1$summary()

samples <- fit1$draws(format="draws_df")
bayesplot::mcmc_trace(samples,"shape")
bayesplot::mcmc_trace(samples,"rate")

acf(samples$shape,lag.max=100)

# Run convergence diagnostics  ------------------
fit1$summary()

