
#  NRES 746, Lecture 5            
#   University of Nevada, Reno                       
#   Optimization     --------------------------                                
#      Searching parameter space to identify the MLE  
#      While thwarting the curse of dimensionality   



# Explore Bolker's myxomatosis example   -------------------------

library(emdbook)    # this is the package provided to support the textbook!
library(ggplot2)
library(ggthemes)

Myx <- MyxoTiter_sum        

ggplot(Myx,aes(day,titer)) + 
  geom_point(aes(col=grade))  +
  facet_wrap(vars(grade), scales = "free") +
  theme_classic() +
  theme(legend.position = "none") 

Myx <- subset(Myx,grade==1)    # subset: select most virulent
head(Myx)


hist(Myx$titer,freq=FALSE)    # distribution of virus loads


# Overlay a gamma distribution on the histogram -------------------

hist(Myx$titer,freq=FALSE)     # note the "freq=FALSE", which displays densities of observations, and therefore makes histograms comparable with probability density functions
curve(dgamma(x,shape=40,rate=6),add=T,col="red")


# Build gamma LL and NLL function  ---------------------

GammaNLL <- function(params){  
  -sum(dgamma(Myx$titer,shape=params['shape'],rate=params['rate'],log=T))     # use params and data to compute likelihood 
}

params <- c(shape=40,rate=6) 
GammaNLL(params)    # test the function!

GammaLL <- function(params){     # same thing- but not using negative LL
  sum(dgamma(Myx$titer,shape=params['shape'],rate=params['rate'],log=T))     # use params and data to compute likelihood 
}


# Optimize using R's built-in "optim()" function: find the maximum likelihood estimate

opt1 <- optim(params,GammaNLL,hessian = T)  
MLE = opt1$par   # store maximum likelihood estimates for params
maxLL = -opt1$value    # store maximum log likelihood (note minus sign)  


# visualize the maximum likelihood fit

hist(Myx$titer,freq=FALSE,xlab="titer")
curve(dgamma(x,shape=MLE["shape"],rate=MLE["rate"]),add=T,col="darkgreen",lwd=2)


# BRUTE FORCE OPTIMIZATION    ------------------------------

# define 2-D parameter space!

shapevec <- seq(0,150,length=100)        # divide parameter space into tiny increments
ratevec <- seq(0.5,30,length=100)

# define the likelihood surface across this grid within parameter space

LLsurface <- expand.grid(shapevec,ratevec)
colnames(LLsurface) <- c("shape","rate")
LLsurface$LL <- sapply(1:nrow(LLsurface), function(t) -GammaNLL(unlist(LLsurface[t,])))   # note minus sign to turn NLL to LL 

summary(LLsurface)

LL_plot_2D <- ggplot(LLsurface,mapping =aes(x=shape,y=rate)) +  # Visualize the likelihood surface
  geom_raster(aes(fill=LL)) +
  geom_contour(aes(z=LL),breaks=seq(maxLL-25,maxLL,10) ,lwd=1.2) +
  scale_fill_gradient(limits=c(maxLL-100,maxLL)) 
LL_plot_2D


# Find the MLE (brute force)  ------------------------

with(LLsurface,c(shape=shape[which.max(LL)],rate=rate[which.max(LL)]) )

MLE  # compare with the answer from "optim()"


# Derivative-based optimization methods   ------------------

library(pracma)   # load package capable of computing gradient numerically at any point in parameter space

grad(GammaLL,MLE)   # confirm that the gradient at the MLE is around zero
grad(GammaLL,unlist(LLsurface[50,1:2]))   # and it's mu


library(dplyr)
param = slice_sample(LLsurface,n=25)
thisgrad = sapply(1:nrow(param),function(t) grad(GammaLL,unlist(param[t,1:2])) ) 
thisgrad = t(apply(thisgrad, 2, function(t) t / Norm(t,2) ) )*2  # normalize by dividing by the magnitude
colnames(thisgrad) <- c("grad_shape","grad_rate")
param <- cbind(param,thisgrad)

# Visualize the gradient of the likelihood function at different points in parameter space
LL_plot_2D +
  geom_segment(data=param, aes(x=shape,y=rate,
                    xend=shape+grad_shape,yend=rate+grad_rate),
               arrow = arrow(length = unit(0.3, "cm")),col="yellow",lwd=2)


# function for estimating the curvature of the likelihood function at any point in parameter space

hessian(GammaLL,MLE)   # confirm the curvature at the maximum likelihood estimate 

-opt1$hessian

hessian(GammaLL,unlist(LLsurface[50,1:2]))   # and here's the curvature at a different point in space... 


# Now we can perform a simple, derivative-based optimization!

start = c(shape=50,rate=7)

thisgrad <- grad(GammaLL,start)
thiscurv <- hessian(GammaLL,start)
thisgrad
thiscurv


# Use this info to estimate the root
newguess <- start - (solve(thiscurv)%*%thisgrad)[,1]
# grad(GammaNLL,newguess)
# hessian(GammaNLL,newguess)
newguess


# Repeat this process

do_newton <- function(oldguess){
  thisgrad <- grad(GammaLL,oldguess); thiscurv <- hessian(GammaLL,oldguess)
  oldguess - (solve(thiscurv)%*%thisgrad)[,1]
}

newguess = do_newton(newguess)
newguess


# again...

newguess = do_newton(newguess)
newguess


# Implement the Newton Method as a function!  ------------------

NewtonMethod <- function(guess,tolerance=0.0000001){
  counter=0
  while(Norm(grad(GammaLL,guess),2)>tolerance){
    guess = do_newton(guess)
    counter=counter+1
  }
  list(
    estimate = guess,
    likelihood = GammaLL(guess),
    iterations = counter
  )
}


newMLE <- NewtonMethod(start)
newMLE


# SIMPLEX OPTIMIZATION METHOD!   -----------------------

# set up an "initial" simplex

guess <- c(shape=60,rate=8)   # "user" first guess 

make_simplex <- function(guess){
  list(
    vertex1 = guess,
    vertex2 = guess + c(10,0),
    vertex3 = guess + c(-5,-1) 
  )
}

thissimplex = make_simplex(guess)

thissimplex
    ## first let's visualize the simplex on a 2-D likelihood surface...

simplex_dat = as.data.frame(do.call(rbind,thissimplex))

LL_plot_2D  + 
  geom_polygon(data=simplex_dat, aes(x=shape,y=rate),lwd=1.5,fill="pink",alpha=0.5)


# Evaluate log-likelihood at each vertex of the simplex

SimplexLik <- function(simplex){
  newvec <- sapply(simplex,GammaLL)   # note use of apply instead of for loop...
  return(newvec)
}

SimplexLik(thissimplex)



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

# Visualize the simplex  ---------------------

trys = list()
trys[[1]] <- simplex_dat
oldsimplex <- thissimplex
newsimplex <- MoveTheSimplex(oldsimplex)
trys[[2]] <-  as.data.frame(do.call(rbind,newsimplex))
trysdat <- do.call(rbind,trys); trysdat$iter = rep(c(1,2),each=3)

LL_plot_2D  + 
  geom_polygon(data=trysdat, aes(x=shape,y=rate,group=iter),lwd=0.5,fill="pink",linetype=1,alpha=0.7,color="black") 
 

# Make another few moves  -------------

for(i in 3:6){
  newsimplex <- MoveTheSimplex(newsimplex)
  trys[[i]] <-  as.data.frame(do.call(rbind,newsimplex))
}

trysdat <- do.call(rbind,trys); trysdat$iter = rep(c(1:6),each=3)

LL_plot_2D  + 
  geom_polygon(data=trysdat, aes(x=shape,y=rate,group=iter),lwd=0.5,fill="pink",linetype=1,alpha=0.7,color="black")
 


# Build a simplex optimization function!  -----------------

SimplexMethod <- function(firstguess,tolerance=0.00001){
  simplex <- make_simplex(firstguess)
  VertexLik <- SimplexLik(simplex)
  counter <- 0
  while(any(abs(diff(VertexLik))>tolerance)){
    simplex <- MoveTheSimplex(simplex)
    VertexLik <- SimplexLik(simplex)
    bestlik <- VertexLik[which.max(VertexLik)]
    counter <- counter+1
  }
  list(estimate = simplex[[1]],
    likelihood = bestlik,
    iterations = counter)
}

SimplexMethod(c(shape=60,rate=8))

MLE


# Simulated annealing!  -----------------------

startingvals <- c(shape=80,rate=7)
startinglik <- GammaLL(startingvals)
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
  oldLik <- GammaLL(oldguess)
  newLik <- GammaLL(newguess)
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

LL_plot_2D +
  geom_path(data=guesses,aes(x=shape,y=rate),lwd=.4,col="black")


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

LL_plot_2D +
  geom_path(data=guesses,aes(x=shape,y=rate),lwd=.4,col="black")



# cool the "temperature" over time and let the algorithm settle down

k <- 100
oldguess <- startingvals
counter <- 0
guesses <- matrix(0,nrow=10000,ncol=2)
colnames(guesses) <- names(startingvals)
MLE_sa <- list(vals=startingvals,lik=GammaLL(startingvals),step=0)
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
  thislik <- GammaLL(oldguess)
  if(thislik>MLE_sa$lik) MLE_sa <- list(vals=oldguess,lik=GammaLL(oldguess),step=counter)
}

# visualize!

LL_plot_2D +
  geom_path(data=guesses,aes(x=shape,y=rate),lwd=.4,col="black") + 
  annotate("point",x=MLE_sa$vals[1],y=MLE_sa$vals[2],col="green",pch=20,cex=3)

MLE_sa


