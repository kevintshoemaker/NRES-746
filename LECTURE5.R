
#  NRES 746, Lecture 5            
#   University of Nevada, Reno                       
#   Optimization     --------------------------                                
#      Searching parameter space to identify the MLE  
#      While thwarting the curse of dimensionality   



# Explore Bolker's myxomatosis example   -------------------------

library(emdbook)    # this is the package provided to support the textbook!
library(ggplot2)
library(ggthemes)

MyxDat <- MyxoTiter_sum         # load Bolker's example data
MyxDat$grade <- as.factor(MyxDat$grade)

ggplot(MyxDat,aes(day,titer)) + 
  geom_point(aes(col=grade))  +
  facet_wrap(vars(grade), scales = "free") +
  theme_clean()

Myx <- subset(MyxDat,grade==1)    # subset: select most virulent
head(Myx)


hist(Myx$titer,freq=FALSE)    # distribution of virus loads


# Overlay a gamma distribution on the histogram -------------------

hist(Myx$titer,freq=FALSE)     # note the "freq=FALSE", which displays densities of observations, and therefore makes histograms comparable with probability density functions
curve(dgamma(x,shape=40,rate=6),add=T,col="red")


# Build gamma likelihood function  ---------------------

GammaLikelihoodFunction <- function(params){           # only one argument (params)- the data are hard-coded here (this is often the case with simple likelihood functions)
  -sum(dgamma(Myx$titer,shape=params['shape'],rate=params['rate'],log=T))     # use params and data to compute likelihood 
}

params <- c(shape=40,rate=6) 
GammaLikelihoodFunction(params)    # test the function!


# Optimize using R's built-in "optim()" function: find the maximum likelihood estimate

MLE <- optim(params,GammaLikelihoodFunction)  

MLE$par
MLE$value


# visualize the maximum likelihood fit

hist(Myx$titer,freq=FALSE)
curve(dgamma(x,shape=MLE$par["shape"],rate=MLE$par["rate"]),add=T,col="red")


# BRUTE FORCE ALTERNATIVE    ------------------------------

# define 2-D parameter space!

shapevec <- seq(10,100,by=0.1)        # divide parameter space into tiny increments
ratevec <- seq(0.5,30,by=0.05)

# define the likelihood surface across this grid within parameter space


surface2D <- matrix(nrow=length(shapevec),ncol=length(ratevec))   # initialize storage variable

newparams <- params
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  for(j in 1:length(ratevec)){
    newparams['rate'] <- ratevec[j]
    surface2D[i,j] <- -1*GammaLikelihoodFunction(newparams)   # compute likelihood for every point in 2-d parameter space
  }
}

# Visualize the likelihood surface

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-250,-35),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-150),add=T)


# Find the MLE (brute force)  ------------------------

ndx <- which(surface2D==max(surface2D),arr.ind=T)  # index of the max likelihood grid cell
shapevec[ndx[,1]]     
ratevec[ndx[,2]]

MLE$par  # compare with the answer from "optim()"


# Derivative-based optimization methods   ------------------

# function for estimating the slope of the likelihood surface at any point in parameter space....

## NOTE: even here I'm using a coarse, brute force method for estimating the first and second derivative of the likelihood function

params <- MLE$par
SlopeFunc <- function(shape_guess,tiny=0.001){      
  params['shape'] <- shape_guess
  high <- GammaLikelihoodFunction(params+c(tiny,0))
  low <- GammaLikelihoodFunction(params-c(tiny,0))
  slope <- (high-low)/(tiny*2)
  return(slope)
}

SlopeFunc(shape_guess=30)    #try it!


# Visualize the slope of the likelihood function at different points in parameter space

shapevec <- seq(10,100,by=0.1)   

# define the likelihood surface

surface1D <- numeric(length(shapevec))   # initialize storage variable

newparams <- params
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  surface1D[i] <- GammaLikelihoodFunction(newparams) 
}

plot(surface1D~shapevec,type="l")
point <- GammaLikelihoodFunction(c(shape=30,MLE$par['rate']))
slope <- SlopeFunc(shape_guess=30)
lines(c(20,40),c(point-slope*10,point+slope*10),col="red")


# function for estimating the curvature of the likelihood function at any point in parameter space

params <- MLE$par
CurvatureFunc <- function(shape_guess,tiny=0.001){
  params['shape'] <- shape_guess
  high <- SlopeFunc(shape_guess+tiny)
  low <- SlopeFunc(shape_guess-tiny)
  curvature <- (high-low)/(tiny*2)   # how much the slope is changing in this region of the function
  return(curvature)
}

CurvatureFunc(shape_guess=30)   # try it!


# First- visualize the gradient of the likelihood function

firstderiv <- numeric(length(shapevec))   # initialize storage variable
for(i in 1:length(shapevec)){
  firstderiv[i] <- SlopeFunc(shapevec[i]) 
}

plot(firstderiv~shapevec,type="l")
abline(h=0,col="red")


# Now we can perform a simple, derivative-based optimization!

### Pick "80" as the starting value

firstderiv <- SlopeFunc(80)           # evaluate the first and second derivatives
secondderiv <- CurvatureFunc(80)
firstderiv
secondderiv


# Use this info to estimate the root

oldguess <- 80
newguess <- oldguess - firstderiv/secondderiv   # estimate the root (where first deriv is zero)
newguess


# Repeat this process

oldguess <- 41.31
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess) 
newguess


# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


# Implement the Newton Method as a function!  ------------------

NewtonMethod <- function(firstguess,tolerance=0.0000001){
  deriv <- SlopeFunc(firstguess)
  oldguess <- firstguess
  counter <- 0
  while(abs(deriv)>tolerance){
    deriv <- SlopeFunc(oldguess)
    newguess <- oldguess - deriv/CurvatureFunc(oldguess)
    oldguess<-newguess
    counter=counter+1
  }
  mle <- list()
  mle$estimate <- newguess
  mle$likelihood <- GammaLikelihoodFunction(c(shape=newguess,MLE$par['rate']))
  mle$iterations <- counter
  return(mle)
}


newMLE <- NewtonMethod(firstguess=80)
newMLE


# SIMPLEX OPTIMIZATION METHOD!   -----------------------

# set up an "initial" simplex

firstguess <- c(shape=70,rate=5)   # "user" first guess 

simplex <- list()
 
           # set up the initial simplex based on the first guess...
simplex[['vertex1']] <- firstguess + c(3,1)
simplex[['vertex2']] <- firstguess + c(-3,-1)
simplex[['vertex3']] <- firstguess + c(3,-1)

simplex


    ## first let's make a function to plot the simplex on a 2-D likelihood surface...

addSimplex <- function(simplex,col="red"){
  temp <- as.data.frame(simplex)    # easier to work with data frame here
  points(x=temp[1,c(1,2,3,1)], y=temp[2,c(1,2,3,1)],type="b",lwd=2,col=col)
}

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-300,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-120),add=T)
addSimplex(simplex)



# Evaluate log-likelihood at each vertex of the simplex

SimplexLik <- function(simplex){
  newvec <- -1*unlist(lapply(simplex,GammaLikelihoodFunction))   # note use of apply instead of for loop...
  return(newvec)
}

SimplexLik(simplex)



# Helper Functions

## this function reflects the worst vertex across the remaining vector

# values <- SimplexLik(simplex)
# oldsimplex=simplex[order(values,decreasing = T)]   # note: must be sorted with worst vertex last
ReflectIt <- function(oldsimplex){
  
  # vertnames <- names(oldsimplex)
  n=length(oldsimplex[[1]])
  centroid <- apply(t(as.data.frame(oldsimplex[1:n])),2,mean)
  
  reflected <- centroid + (centroid - oldsimplex[[n+1]])
  expanded <- centroid + 2*(centroid - oldsimplex[[n+1]])
  contracted <- centroid + 0.5*(centroid - oldsimplex[[n+1]])
  
  alternates <- list()
  alternates$reflected <- oldsimplex
  alternates$expanded <- oldsimplex 
  alternates$contracted <- oldsimplex 
  alternates$reflected[[n+1]] <- reflected
  alternates$expanded[[n+1]] <- expanded
  alternates$contracted[[n+1]] <- contracted
  return(alternates)
}
# ReflectIt(oldsimplex)


ShrinkIt <- function(oldsimplex){
  n <- length(oldsimplex[[1]])
  X.vert <- t(as.data.frame(oldsimplex[(1:(n+1))]))
  temp <- sweep(0.5*sweep(X.vert, 2, oldsimplex[[1]], FUN = "-"), 2, X.vert[1, ], FUN="+")
  temp2 <- as.data.frame(t(temp))
  lapply(temp2,function(t) c(shape=t[1],rate=t[2])  )
}


MoveTheSimplex <- function(oldsimplex){     # (incomplete) nelder-mead algorithm
  newsimplex <- oldsimplex  # 
           # Start by sorting the simplex (worst vertex last)
  VertexLik <- SimplexLik(newsimplex)
  newsimplex <- newsimplex[order(VertexLik,decreasing=T)]
  liks <- VertexLik[order(VertexLik,decreasing=T)]
  worstLik <- liks[3]
  secondworstLik <- liks[2]
  bestLik <- liks[1]
  
  candidates <- ReflectIt(oldsimplex=newsimplex)      # reflect across the remaining edge
  CandidateLik <- sapply(candidates,SimplexLik)                          # re-evaluate likelihood at the vertices...
  CandidateLik <- apply(CandidateLik,c(1,2), function(t) ifelse(is.nan(t),-99999,t))
  bestCandidate <- names(which.max(CandidateLik[3,]))
  bestCandidateLik <- CandidateLik[3,bestCandidate]
  
  if((CandidateLik[3,"reflected"]<=bestLik)&(CandidateLik[3,"reflected"]>secondworstLik)){
    newsimplex <- candidates[["reflected"]]
  }else if (CandidateLik[3,"reflected"]>bestLik){
    if(CandidateLik[3,"expanded"]>CandidateLik[3,"reflected"]){
      newsimplex <- candidates[["expanded"]]
    }else{
      newsimplex <- candidates[["reflected"]]
    }
  }else{
    if(CandidateLik[3,"contracted"]>worstLik){
      newsimplex <- candidates[["contracted"]]
    }else{
      newsimplex <- ShrinkIt(newsimplex)
    }
  }

  return(newsimplex)
}

# image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
# contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
# addSimplex(oldsimplex,col="red")
# addSimplex(candidates$reflected,col="green")
# addSimplex(candidates$half,col="green")

# Visualize the simplex  ---------------------

oldsimplex <- simplex
newsimplex <- MoveTheSimplex(oldsimplex)
image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-125),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex,col="green")



# Make another move  -------------

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-125),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex,col="green")


# Make another move  ----------------------

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-125),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex,col="green")


# Make another move  ----------------

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-125),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex,col="green")


# Make another few moves  ----------------------

par(mfrow=c(2,2))

for(i in 1:4){
  oldsimplex <- newsimplex
  newsimplex <- MoveTheSimplex(oldsimplex)
  
  image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
  contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-125),add=T)
  addSimplex(oldsimplex,col="red")
  addSimplex(newsimplex,col="green")
}



# Build a simplex optimization function!  -----------------

SimplexMethod <- function(firstguess,tolerance=0.00001){
  initsimplex <- list()
  initsimplex[['vertex1']] <- firstguess + c(5,0.5)
  initsimplex[['vertex2']] <- firstguess + c(-5,-0.5)
  initsimplex[['vertex3']] <- firstguess + c(5,-0.5)
  VertexLik <- SimplexLik(initsimplex)
  oldbestlik <- VertexLik[which.max(VertexLik)]
  deltalik <- 100
  counter <- 0
  oldsimplex <- initsimplex
  while((counter<250)&(any(abs(diff(VertexLik))>tolerance))){
    newsimplex <- MoveTheSimplex(oldsimplex)
    VertexLik <- SimplexLik(newsimplex)
    bestlik <- VertexLik[which.max(VertexLik)]
    oldsimplex <- newsimplex
    counter <- counter+1
  }
  mle <- list()
  mle$estimate <- newsimplex[[1]]
  mle$likelihood <- bestlik
  mle$iterations <- counter
  return(mle)
}

SimplexMethod(firstguess = c(shape=39,rate=4))


# Simulated annealing!  -----------------------

startingvals <- c(shape=80,rate=7)
startinglik <- -GammaLikelihoodFunction(startingvals)
startinglik

k = 100   # set the "temperature"
 
     # function for making new guesses
newGuess <- function(oldguess=startingvals){
  maxshapejump <- 5
  maxratejump <- 0.75
  jump <- c(runif(1,-maxshapejump,maxshapejump),runif(1,-maxratejump,maxratejump))
  newguess <- oldguess + jump
  return(newguess)
}
  # set a new "guess" near to the original guess

newGuess(oldguess=startingvals)     # each time is different- this is the first optimization procedure with randomness built in
newGuess(oldguess=startingvals)
newGuess(oldguess=startingvals)


# evaluate the difference in likelihood between the new proposal and the old point

LikDif <- function(oldguess,newguess){
  oldLik <- -GammaLikelihoodFunction(oldguess)
  newLik <- -GammaLikelihoodFunction(newguess)
  return(newLik-oldLik)
}

newguess <- newGuess(oldguess=startingvals)
loglikdif <- LikDif(oldguess=startingvals,newguess)
loglikdif


# run and visualize a Metropolis simulated annealing routine -------------

k <- 100
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=100,ncol=2)
colnames(guesses) <- names(startingvals)
while(counter<100){
  newguess <- newGuess(oldguess)
  loglikdif <- LikDif(oldguess,newguess)
  if(loglikdif>0){ 
    oldguess <- newguess
  }else{
    rand=runif(1)
    if(rand <= exp(loglikdif/k)){
      oldguess <- newguess   # accept even if worse!
    }
  }
  counter <- counter + 1
  guesses[counter,] <- oldguess
}

# visualize!

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-135),add=T)
lines(guesses,col="red")


# Run it for longer!

k <- 10
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=1000,ncol=2)
colnames(guesses) <- names(startingvals)
while(counter<1000){
  newguess <- newGuess(oldguess)
  while(any(newguess<0)) newguess <- newGuess(oldguess)
  loglikdif <- LikDif(oldguess,newguess)
  if(loglikdif>0){ 
    oldguess <- newguess
  }else{
    rand=runif(1)
    if(rand <= exp(loglikdif/k)){
      oldguess <- newguess   # accept even if worse!
    }
  }
  counter <- counter + 1
  guesses[counter,] <- oldguess
}

# visualize!

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-135),add=T)
lines(guesses,col="red")



# cool the "temperature" over time and let the algorithm settle down

k <- 100
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=10000,ncol=2)
colnames(guesses) <- names(startingvals)
MLE <- list(vals=startingvals,lik=-GammaLikelihoodFunction(startingvals),step=0)
while(counter<10000){
  newguess <- newGuess(oldguess)
  while(any(newguess<0)) newguess <- newGuess(oldguess)
  loglikdif <- LikDif(oldguess,newguess)
  if(loglikdif>0){ 
    oldguess <- newguess
  }else{
    rand=runif(1)
    if(rand <= exp(loglikdif/k)){
      oldguess <- newguess   # accept even if worse!
    }
  }
  counter <- counter + 1
  if(counter%%100==0) k <- k*0.8
  guesses[counter,] <- oldguess
  thislik <- -GammaLikelihoodFunction(oldguess)
  if(thislik>MLE$lik) MLE <- list(vals=oldguess,lik=-GammaLikelihoodFunction(oldguess),step=counter)
}

# visualize!

image(x=shapevec,y=ratevec,z=surface2D,zlim=c(-500,-30),col=topo.colors(12))
contour(x=shapevec,y=ratevec,z=surface2D,levels=c(-30,-40,-80,-135),add=T)
lines(guesses,col="red")
points(MLE$vals[1],MLE$vals[2],col="green",pch=20,cex=3)

MLE

optim(params,GammaLikelihoodFunction)$par

