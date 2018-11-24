
############################################################
####                                                    ####  
####  NRES 746, Lecture 5                               ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Optimization
####     Searching parameter space to identify the MLE  ####
####     While thwarting the curse of dimensionality    ####
############################################################



####################
# Explore Bolker's myxomatosis example

library(emdbook)    # this is the package provided to support the textbook!

MyxDat <- MyxoTiter_sum         # load Bolker's example data
Myx <- subset(MyxDat,grade==1)
head(Myx)


hist(Myx$titer,freq=FALSE)    # distribution of virus loads


###########
# Overlay a gamma distribution on the histogram 

hist(Myx$titer,freq=FALSE)     # note the "freq=FALSE", which displays densities of observations, and therefore makes histograms comparable with probability densities
curve(dgamma(x,shape=40,scale=0.15),add=T,col="red")


################
# Build likelihood function

GammaLikelihoodFunction <- function(params){           # only one argument (params)- the data are hard-coded here (this is often the case with simple likelihood functions)
  sum(dgamma(Myx$titer,shape=params['shape'],scale=params['scale'],log=T))     # use params and data to compute likelihood 
}

params <- c(40,0.15) 
names(params) <- c("shape","scale")
params
GammaLikelihoodFunction(params)    # test the function!


############
# USE R's 'OPTIM()' FUNCTION
############

############
# Optimize using R's built-in "optim()" function: find the maximum likelihood estimate

ctrl <- list(fnscale=-1)   # maximize rather than minimize!!
MLE <- suppressWarnings(optim(fn=GammaLikelihoodFunction,par=params,control=ctrl,method="BFGS"))   # stop the warnings!

MLE$par


##############
# visualize the fit

hist(Myx$titer,freq=FALSE)
curve(dgamma(x,shape=MLE$par["shape"],scale=MLE$par["scale"]),add=T,col="red")


######################
# BRUTE FORCE ALTERNATIVE
######################

##############
# define 2-D parameter space!
##############

shapevec <- seq(10,100,by=0.1)   
scalevec <- seq(0.01,0.3,by=0.001)

##############
# define the likelihood surface across this grid within parameter space
##############

surface2D <- matrix(nrow=length(shapevec),ncol=length(scalevec))   # initialize storage variable

newparams <- params
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  for(j in 1:length(scalevec)){
    newparams['scale'] <- scalevec[j]
    surface2D[i,j] <- GammaLikelihoodFunction(newparams)   # compute likelihood for every point in 2-d parameter space
  }
}

############
# Visualize the likelihood surface
############

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)


############
# Find the MLE
############

ndx <- which(surface2D==max(surface2D),arr.ind=T)  # index of the max likelihood grid cell
shapevec[ndx[,1]]     
scalevec[ndx[,2]]

MLE$par  # compare with the answer from "optim()"


###################
# Derivative-based optimization methods
###################

######
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


#########
# Visualize the slope of the likelihood function at different points in parameter space

shapevec <- seq(10,100,by=0.1)   

##############
# define the likelihood surface
##############

surface1D <- numeric(length(shapevec))   # initialize storage variable

newparams <- params
for(i in 1:length(shapevec)){
  newparams['shape'] <- shapevec[i]
  surface1D[i] <- GammaLikelihoodFunction(newparams) 
}

plot(surface1D~shapevec,type="l")
point <- GammaLikelihoodFunction(c(shape=30,MLE$par['scale']))
slope <- SlopeFunc(shape_guess=30)
lines(c(20,40),c(point-slope*10,point+slope*10),col="red")


########
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


######
# First- visualize the gradient of the likelihood function

firstderiv <- numeric(length(shapevec))   # initialize storage variable
for(i in 1:length(shapevec)){
  firstderiv[i] <- SlopeFunc(shapevec[i]) 
}

plot(firstderiv~shapevec,type="l")
abline(h=0,col="red")


#########
# Now we can perform a simple, derivative-based optimization!

### Pick "80" as the starting value

firstderiv <- SlopeFunc(80)           # evaluate the first and second derivatives
secondderiv <- CurvatureFunc(80)
firstderiv
secondderiv


#########
# Use this info to estimate the root

oldguess <- 80
newguess <- oldguess - firstderiv/secondderiv   # estimate the root (where first deriv is zero)
newguess


##########
# Repeat this process

oldguess <- 41.31
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess) 
newguess


#######
# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


#######
# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


#######
# again...

oldguess<-newguess
newguess <- oldguess - SlopeFunc(oldguess)/CurvatureFunc(oldguess)
newguess


##########
# Implement the Newton Method as a function!

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
  mle$likelihood <- GammaLikelihoodFunction(c(shape=newguess,MLE$par['scale']))
  mle$iterations <- counter
  return(mle)
}


newMLE <- NewtonMethod(firstguess=80)
newMLE


#############
# SIMPLEX OPTIMIZATION METHOD!
#############

########
# set up an "initial" simplex

firstguess <- c(shape=40,scale=0.25)   # "user" first guess 

simplex <- list()
 
           # set up the initial simplex based on the first guess...
simplex[['vertex1']] <- firstguess + c(5,0.05)
simplex[['vertex2']] <- firstguess + c(-5,-0.05)
simplex[['vertex3']] <- firstguess + c(5,-0.05)

simplex


    ## first let's make a function to plot the simplex on a 2-D likelihood surface...

addSimplex <- function(simplex,col="green"){
  temp <- as.data.frame(simplex)    # easier to work with data frame here
  points(x=temp[1,c(1,2,3,1)], y=temp[2,c(1,2,3,1)],type="b",lwd=2,col=col)
}

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
addSimplex(simplex)



########
# Evaluate log-likelihood at each vertex of the simplex

SimplexLik <- function(simplex){
  newvec <- unlist(lapply(simplex,GammaLikelihoodFunction))   # note use of apply instead of for loop...
  return(newvec)
}

SimplexLik(simplex)


#####
# Helper Functions
#####

## this function relects the worst vertex across the remaining vector
ReflectIt <- function(oldsimplex,WorstVertex){
  
    ## re-arrange simplex- worst must be first
  worstndx <- which(names(oldsimplex)==WorstVertex)
  otherndx <- c(1:3)[-worstndx]
  newndx <- c(worstndx,otherndx) 
  
    ## translate so that vertex 2 is the origin (0,0)
  translate <- oldsimplex[[newndx[2]]]
  newsimplex <- list(oldsimplex[[1]]-translate,oldsimplex[[2]]-translate,oldsimplex[[3]]-translate)
  
    ## use vector reflection (reflect vertex 1 over a vector containing the origin and vertex 3) to find the reflection across a vector that includes the origin
  vdotl <- sum(newsimplex[[newndx[1]]]*newsimplex[[newndx[3]]])
  ldotl <- sum(newsimplex[[newndx[3]]]*newsimplex[[newndx[3]]])
  
  projection <- (vdotl/ldotl)*newsimplex[[newndx[3]]]
  
  reflected <- 2*projection-newsimplex[[newndx[1]]]
  
    ## translate back to the likelihood surface
  newsimplex[[newndx[1]]] <- reflected
  newsimplex <- list(newsimplex[[1]]+translate,newsimplex[[2]]+translate,newsimplex[[3]]+translate)
    ## return the new simplex
  names(newsimplex) <- names(oldsimplex)
  
    ## generate some alternative jumps (or "oozes"!)...
  oldpoint <- oldsimplex[[worstndx]]
  newpoint <- newsimplex[[worstndx]]
  
  newpoint2 <- newpoint-oldpoint
  double <- newpoint2 * 2
  half <- newpoint2 * 0.25
  
  alternates <- list()
  alternates$reflected <- newsimplex
  alternates$double <- newsimplex 
  alternates$half <- newsimplex 
  alternates$double[[worstndx]] <- double + oldpoint
  alternates$half[[worstndx]] <- half + oldpoint
  return(alternates)
}


ShrinkIt <- function(oldsimplex,BestVertex){
  newsimplex <- oldsimplex
  
      ## indices...
  bestndx <- which(names(oldsimplex)==BestVertex)
  otherndx <- c(1:3)[-bestndx]
  
  translate <- oldsimplex[[bestndx]]
  
  i=2
  for(i in otherndx){
    newvector <- oldsimplex[[i]]-translate
    shrinkvector <- newvector * 0.5
    newsimplex[[i]] <- shrinkvector + translate
  }
  
  return(newsimplex)
}


MoveTheSimplex <- function(oldsimplex){     # (incomplete) nelder-mead
  newsimplex <- oldsimplex  # 
           # Start by identifying the *worst* vertex (the one with the lowest likelihood)
  VertexLik <- SimplexLik(newsimplex)
  WorstLik <- min(VertexLik)
  WorstVertex <- names(VertexLik[which.min(VertexLik)])    # identify vertex with lowest likelihood
  candidates <- ReflectIt(oldsimplex=newsimplex,WorstVertex)      # reflect across the remaining edge
  CandidateLik <- sapply(candidates,SimplexLik)                          # re-evaluate likelihood at the vertices...
  CandidateLik <- apply(CandidateLik,c(1,2), function(t) ifelse(is.nan(t),-99999,t))
  bestCandidate <- names(which.max(CandidateLik[WorstVertex,]))
  bestCandidateLik <- CandidateLik[WorstVertex,bestCandidate]
  if(bestCandidateLik>=WorstLik){
    newsimplex <- candidates[[bestCandidate]]
  } else{
    BestVertex <- names(VertexLik[which.max(VertexLik)])
    newsimplex <- ShrinkIt(oldsimplex,BestVertex)
  }
  return(newsimplex)
}

###########
# Visualize the simplex

oldsimplex <- simplex
newsimplex <- MoveTheSimplex(oldsimplex)
image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex)



############
# Make another move

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex)


############
# Make another move

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex)


############
# Make another move

oldsimplex <- newsimplex
newsimplex <- MoveTheSimplex(oldsimplex)

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
addSimplex(oldsimplex,col="red")
addSimplex(newsimplex)


############
# Make another few moves

par(mfrow=c(2,2))

for(i in 1:4){
  oldsimplex <- newsimplex
  newsimplex <- MoveTheSimplex(oldsimplex)
  
  image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
  contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
  addSimplex(oldsimplex,col="red")
  addSimplex(newsimplex)
}



############
# Build a simplex optimization function!

SimplexMethod <- function(firstguess,tolerance=0.00001){
  initsimplex <- list()
  initsimplex[['vertex1']] <- firstguess + c(5,0.05)
  initsimplex[['vertex2']] <- firstguess + c(-5,-0.05)
  initsimplex[['vertex3']] <- firstguess + c(5,-0.05)
  VertexLik <- SimplexLik(initsimplex)
  oldbestlik <- VertexLik[which.max(VertexLik)]
  deltalik <- 100
  counter <- 0
  while(counter<100){
    newsimplex <- MoveTheSimplex(oldsimplex)
    VertexLik <- SimplexLik(newsimplex)
    bestlik <- VertexLik[which.max(VertexLik)]
    deltalik <- bestlik-oldbestlik
    oldsimplex <- newsimplex
    oldbestlik <- bestlik
    counter <- counter+1
  }
  mle <- list()
  mle$estimate <- newsimplex[[1]]
  mle$likelihood <- bestlik
  mle$iterations <- counter
  return(mle)
}


SimplexMethod(firstguess = c(shape=39,scale=0.28))


################
# Simulated annealing!

startingvals <- c(shape=80,scale=0.15)
startinglik <- GammaLikelihoodFunction(startingvals)
startinglik

k = 100   # set the "temperature"
 
     # function for making new guesses
newGuess <- function(oldguess=startingvals){
  maxshapejump <- 5
  maxscalejump <- 0.05
  jump <- c(runif(1,-maxshapejump,maxshapejump),runif(1,-maxscalejump,maxscalejump))
  newguess <- oldguess + jump
  return(newguess)
}
  # set a new "guess" near to the original guess

newGuess(oldguess=startingvals)     # each time is different- this is the first optimization procedure with randomness built in
newGuess(oldguess=startingvals)
newGuess(oldguess=startingvals)


############
# evaluate the difference in likelihood between the new proposal and the old point

LikDif <- function(oldguess,newguess){
  oldLik <- GammaLikelihoodFunction(oldguess)
  newLik <- GammaLikelihoodFunction(newguess)
  return(newLik-oldLik)
}

newguess <- newGuess(oldguess=startingvals)
loglikdif <- LikDif(oldguess=startingvals,newguess)
loglikdif


############
# run and visualize a Metropolis routine

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

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")


###########
# Run it for longer!

k <- 10
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=1000,ncol=2)
colnames(guesses) <- names(startingvals)
while(counter<1000){
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

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")



#############
# cool the "temperature" over time and let the algorithm settle down

k <- 100
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=10000,ncol=2)
colnames(guesses) <- names(startingvals)
MLE <- list(vals=startingvals,lik=GammaLikelihoodFunction(startingvals),step=0)
while(counter<10000){
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
  if(counter%%100==0) k <- k*0.8
  guesses[counter,] <- oldguess
  thislik <- GammaLikelihoodFunction(oldguess)
  if(thislik>MLE$lik) MLE <- list(vals=oldguess,lik=GammaLikelihoodFunction(oldguess),step=counter)
}

# visualize!

image(x=shapevec,y=scalevec,z=surface2D,zlim=c(-1000,-30),col=topo.colors(12))
contour(x=shapevec,y=scalevec,z=surface2D,levels=c(-30,-40,-80,-500),add=T)
lines(guesses,col="red")
points(MLE$vals[1],MLE$vals[2],col="green",pch=20,cex=3)

MLE

