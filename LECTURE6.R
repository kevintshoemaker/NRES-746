
#  NRES 746, Lecture 6                             
#   Bayesian analysis #1: concepts   -------------------------             


# Bayesian analysis example using Binomial distribution  -----------------------

# first visualize the "prior" for the probability *p* as a uniform distribution

ggplot() +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1.5)) +
  stat_function(fun=dbeta,args = list(shape1 = 1, shape2 = 1),lwd=2,col=gray(.4)) +
  labs(x="parameter \"p\"",y="probability") +
  theme_classic()


# frog call example: --------------

#    imagine we detected the frog in 3 of 10 visits to a known-occupied wetland
#    visualize the data likelihood alongside the prior probability
#    recall that the likelihood surface is not a probability distribution (thus, the 2 y axes)

lik = function(p) dbinom(3,10,p)

ggplot() +
  xlim(0,1) + ylim(0,2) + 
  stat_function(fun=dbeta,args = list(shape1 = 1, shape2 = 1),lwd=1.5,col=gray(.4)) +
  stat_function(fun=lik,lwd=1.5,col="darkorange") +
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


# Brute-force Bayes ----------------------------

# prior across parameter space 

prior <- function(p) dbeta(p,1,1)    # flat prior

## Numerator for Bayes rule: weight the data likelihood by the prior
numer <- function(p) lik(p)*prior(p)      # Numerator for Bayes rule

## Denominator for Bayes rule: compute normalization constant
marg_lik <- function() integrate(numer,0,1)$value   # marginal likelihood is a single constant number: the sum of the weighted likelihoods!

## Posterior (numerator/denominator)
posterior <- function(p) numer(p)/marg_lik()   # this is Bayes' rule!

## Plot it out!

ggplot() +
  xlim(0,1) +
  stat_function(fun=prior,lwd=1.5,col=gray(.4)) +
  stat_function(fun=lik,lwd=1.5,col="darkorange") +
  stat_function(fun=posterior,lwd=1.5,col="darkgreen") +
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


# Try an informative prior!  -------------

prior <- function(p) dbeta(p,15,5)    # more informative prior

ggplot() +
  xlim(0,1) +
  stat_function(fun=prior,lwd=1.5,col=gray(.4)) +
  stat_function(fun=lik,lwd=1.5,col="darkorange") +
  stat_function(fun=posterior,lwd=1.5,col="darkgreen") +
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()

# Collect more data and try again...   -----------

dat <- c(3, 1, 6, 2, 3, 2, 6, 1, 3, 3)
lik = function(p) sapply(p, function(t) prod(dbinom(dat,10,t)) )  # vectorized with sapply

ggplot() +
  xlim(0,1) +
  stat_function(fun=prior,lwd=1.5,col=gray(.4)) +
  # stat_function(fun=lik,lwd=1.5,col="darkorange") +    # likelihood is too small to show up
  stat_function(fun=posterior,lwd=1.5,col="darkgreen") +
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


prior <- function(p) dbeta(p,150,50)    # super informative prior

ggplot() +
  xlim(0,1) +
  stat_function(fun=prior,lwd=1.5,col=gray(.4)) +
  stat_function(fun=posterior,lwd=1.5,col="darkgreen") +
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()

## Conjugate priors ---------------------------

# Do it again- this time with conjugate priors...

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 1, shape2 = 1),lwd=1.5,col=gray(.4)) +  # PRIOR
  stat_function(fun=dbeta,args = list(shape1 = 1+3, shape2 = 1+(10-3)),lwd=1.5,col="darkgreen") +  # POSTERIOR (after observing 3 successes out of 10)
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


# With informative prior...

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 15, shape2 = 5),lwd=1.5,col=gray(.4)) +  # PRIOR
  stat_function(fun=dbeta,args = list(shape1 = 15+3, shape2 = 5+(10-3)),lwd=1.5,col="darkgreen") +  # POSTERIOR (after observing 3 successes out of 10)
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


graphics.off()


# And with super informative prior...

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 150, shape2 = 50),lwd=1.5,col=gray(.4)) +  # PRIOR
  stat_function(fun=dbeta,args = list(shape1 = 150+3, shape2 = 50+(10-3)),lwd=1.5,col="darkgreen") +  # POSTERIOR (after observing 3 successes out of 10)
  labs(x="parameter \"p\"",y="probability/likelihood") +
  theme_classic()


## Bayesian point estimate -----------------------
# Example: bayesian point estimates can differ from MLE

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 0.5+1, shape2 = 0.5+2),lwd=1.5,col="darkgreen") +  # skewed posterior
  labs(x="parameter",y="probability/likelihood") +
  theme_classic()


# Compute and plot the mean and the mode of the distribution

posterior = function(p) dbeta(p,1.5,2.5)
mean <- integrate(function(p){posterior(p)*p},0,1)$value
mode <- optimize(posterior,c(0,1),maximum=T)$maximum

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 0.5+1, shape2 = 0.5+2),lwd=1.5,col="darkgreen") +  # skewed posterior
  geom_vline(xintercept=c(mean,mode),col=c("red","blue"),lwd=2,lty=2) +
  labs(x="parameter",y="probability/likelihood") +
  theme_classic()


# A Bayesian confidence interval (95% credible interval)  -------------------

credible.interval <- qbeta(c(0.025,0.975),1+3,1+(10-3))     # get the credible interval using the quantile method

ggplot() +
  xlim(0,1) +
  stat_function(fun=dbeta,args = list(shape1 = 1+3, shape2 = 1+(10-3)),lwd=1.5,col="darkgreen") +
  geom_vline(xintercept=credible.interval,col="blue",lwd=1,lty=2) +
  labs(x="p",y="probability/likelihood") +
  theme_classic()


## Bayesian analysis without a conjugate prior

# Revisit the Myxomatosis example  --------------------------

library(emdbook)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)
head(Myx)


hist(Myx$titer,freq=FALSE)    # visualize the data, again!


### Error is modeled as gamma distributed

ggplot(Myx, aes(x = titer)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, fill = "lightblue", color = "black") +
      stat_function(fun=dgamma,args = list(shape = 40, rate = 6),lwd=1,lty=2,col="darkgreen") +
      labs(x = "Titer", y = "Density") +
      theme_classic()


# recall our likelihood function for these data (not on log scale this time!)

lik <- function(pars){
  prod(dgamma(Myx$titer,shape=pars['shape'],rate=pars['rate']))
}

pars <- c(shape=40,rate=6)    # test the function
lik(pars)


# define 2-D parameter space (in real probability scale)!

shapevec <- seq(0,150,length=100)        # divide parameter space into tiny increments
ratevec <- seq(0.5,30,length=100)

parmsurface <- expand.grid(shapevec,ratevec)
colnames(parmsurface) <- c("shape","rate")
parmsurface$lik <- sapply(1:nrow(parmsurface), function(t) lik(unlist(parmsurface[t,1:2]))  )

likplot_2D <- ggplot(parmsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the likelihood surface
  geom_raster(aes(fill=lik)) +
  scale_fill_gradient(limits=c(1e-70,1e-17)) 

likplot_2D


# compute the area of each pixel in parameter space (for probability density computation)
pixelArea <- diff(range(shapevec))/100 * diff(range(ratevec))/100  
npixels <- 100*100

# define the prior probability surface across this grid within parameter space

prior <- function(pars) replicate(nrow(pars),pixelArea/(npixels*pixelArea) )    # set as uniform across the support
parmsurface$prior <- prior(as.matrix(parmsurface[,c("shape","rate")]) )


# Apply Bayes Rule!   -------------------

parmsurface$numer <- with(parmsurface, lik*prior  )    # numerator of Bayes rule
denom <- sum(parmsurface$numer)     # denominator of Bayes rule

parmsurface$post <- parmsurface$numer/denom    # apply Bayes rule

# Visualize the 2-D posterior distribution

postplot_2D <- ggplot(parmsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the likelihood surface
  geom_raster(aes(fill=post)) +
  scale_fill_gradient(limits=c(1e-25,0.02562)) 

postplot_2D


# find the contour that encloses approx. 95% of our degree of belief!  ----------------

cntrconf <- function(c) with(parmsurface, sum(post[which(post>=c)]) )
find_cntr <- function(c,conf) (cntrconf(c)-conf)^2  
c = optimize(find_cntr, interval=c(1e-10,0.02562),conf=0.95)$objective 
c


# Visualize the 2D credible region (HPD)

postplot_2D <- postplot_2D +
  geom_contour(aes(z=post),breaks=6.824e-06 ,lwd=1.2) 
postplot_2D

# Visualize a point estimate

meanpars <- with(parmsurface, c(shape=sum(shape*post),rate=sum(rate*post) ) ) 
modepars <- unlist(parmsurface[which.max(parmsurface$post),c("shape","rate")] ) 
pts <- as.data.frame(rbind(meanpars,modepars)); pts$method=c("mean","mode")

postplot_2D + geom_point(data=pts,aes(x=shape,y=rate,col=method),size=3)


library(cowplot)

shape_marginal = parmsurface |> 
  group_by(shape) |> 
  summarise(post = sum(post))

rate_marginal = parmsurface |> 
  group_by(rate) |> 
  summarise(post = sum(post))

# Plot out posterior distributions separately for each parameter  -------------

g1 = ggplot(shape_marginal,aes(shape,post)) + geom_path(lwd=2) + theme_classic() + ylab("density")
g2 = ggplot(rate_marginal,aes(rate,post)) + geom_path(lwd=2) + theme_classic() + ylab("density")
# g1

cowplot::plot_grid(g1,g2)


cdf_shape <- cumsum(shape_marginal$post)
cdf_rate <- cumsum(rate_marginal$post)

meanshape = with(shape_marginal, sum(shape*post) ) ; meanrate = with(rate_marginal, sum(rate*post) )
ci95shape = shapevec[c(tail(which(cdf_shape<0.025),1),tail(which(cdf_shape<0.975),1) )]
ci95rate = ratevec[c(tail(which(cdf_rate<0.025),1),tail(which(cdf_rate<0.975),1) )]
cis = rbind(ci95shape,ci95rate); colnames(cis) <- c("lb","ub")

stats = data.frame(parm=c("shape","rate"), mean=c(meanshape,meanrate))
stats = cbind(stats,cis ) 
stats


# Sample parameters from the joint posterior

SampleFromPosterior <- function(n){
  samples <- sample(c(1:nrow(parmsurface)),size=n,replace=TRUE,prob=parmsurface$post)
  parmsurface[samples,c("shape","rate")]
}

samples<-SampleFromPosterior(n=10000)
par(mfrow=c(2,2))
plot(ts(samples[,1]),xlab="sample",ylab="shape")
plot(ts(samples[,2]),xlab="sample",ylab="rate")
hist(samples[,1],40,xlab="shape",main="histogram of shape param")
hist(samples[,2],40,xlab="scale",main="histogram of rate param")
par(mfrow=c(1,1))

