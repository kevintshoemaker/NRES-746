
############################################################
####                                                    ####  
####  NRES 746, Lab 3                                   ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  DIY Likelihood Functions                          ####
############################################################



#########
# Challenge 1a

NLL_frogOccupancy <- function(params=0.5,data=c(3,2,6),N=10){
  -sum(dbinom(data,prob=params,size = N,log=T))
}


NLL_frogOccupancy(params=0.5,data=c(3,2,6),N=10)   # test your function


###########
# 1b

xvals <- seq(0.001,0.999,0.001)
nlls <- sapply(1:length(xvals),function(t) NLL_frogOccupancy(xvals[t],data=c(3,2,6),N=10) ) 
plot(nlls~xvals,xlab="parameter \"p\"",ylab="neg log lik",type="l")

MLE <- xvals[which.min(nlls)]
ML <- min(nlls)

CI95 <- range(xvals[nlls<=(ML+2)])

abline(v=MLE,col="green",lwd=2)
abline(h=ML+2,col="blue")
abline(v=CI95,col="green",lty=2)



#####
# 2a

rffr <- read.csv("ReedfrogFuncResp.csv",row.names = 1)
  # alternative: data(ReedfrogFuncresp)     # from Bolker's "emdbook" package
  # ?Reedfrog      # learn more about this dataset
head(rffr)


#########
# define a Holling type II functional response, with an initial guess about parameter values

Holl2<-function(x, a, h){(a*x)/(1+(a*h*x))}
plot(rffr$Killed~rffr$Initial)
curve(Holl2(x, a=0.5, h=1/80), add=TRUE,col="red")


###########
# Write a likelihood function

#    params: vector of params to estimate (a and h from the Holling type II functional response)
#    k: number killed per trial   (data)
#    N: number of tadpoles per trial (data)

binomNLL2<-function(params,N,k){
	a=params[1]
	h=params[2]
	predprob=a/(1+a*h*N)	
	-sum(dbinom(k,prob=predprob,size=N,log=TRUE))
}


quantilefunc <- function(x,pars,q){
  ifelse(x>0.5,
    qbinom(q,prob=Holl2(round(x),pars["a"],pars["h"])/round(x),size=round(x)),
    0
  )
}

inits <- c(a=0.6,h=(1/60))
Rffuncresp <- function(params=inits,data=rffr){
  temp <- suppressWarnings( optim(fn=binomNLL2,  par=inits, N=data[,1], k=data[,2])   )
  MLE <- temp$par
  plot(data[,2]~data[,1],ylab="Killed",xlab="Init Density")
  curve(Holl2(x,inits["a"],inits["h"]),add=T,col="red")
  curve(Holl2(x,MLE["a"],MLE["h"]),add=T,col="green",lwd=2)
  curve(quantilefunc(x,pars=MLE,0.975),add=T,col="green",lty=2)
  curve(quantilefunc(x,pars=MLE,0.025),add=T,col="green",lty=2)
  return(MLE)
}


inits <- c(a=0.6,h=(1/60))    # test the function
Rffuncresp(params=inits,data=rffr)


#######
# 3a


########
# Myxomatosis data

library(emdbook)
data(MyxoTiter_sum)      # load the data
head(MyxoTiter_sum)   


myxdat <- subset(MyxoTiter_sum, grade==1)    # select just the most virulent strain

plot(myxdat$titer~myxdat$day,xlim=c(0,10))    # visualize the relationship



Ricker <- function(x,a,b){
  a*x*exp(-b*x)
}

NLL_myxRicker <- function(params=c(a=1,b=0.2,shape=1),data=myxdat[,-1]){
  exp_titer <- Ricker(data[,1],params["a"],params["b"])
  scale <- exp_titer/params["shape"]
  -sum(dgamma(data[,2],shape=params["shape"],scale=scale,log=TRUE))
}


NLL_myxRicker(params=c(a=4,b=0.2,shape=40),data=myxdat[,-1])   # test the function


##########
# 3b

MyxRicker <- function(params=c(a=2,b=0.2,shape=30),data=myxdat[,-1]){
  temp <- optim(NLL_myxRicker,par = params,data=data)
  MLE <- temp$par
  plot(data[,2]~data[,1],xlab="days",ylab="titer",xlim=c(0,10),ylim=c(0,10))
  curve(Ricker(x,MLE["a"],MLE["b"]),add=T,col="green",lwd=2)
  curve(qgamma(0.975,shape=MLE["shape"],scale=Ricker(x,MLE["a"],MLE["b"])/MLE["shape"]),0.01,10,add=T,col="green",lty=2)
  curve(qgamma(0.025,shape=MLE["shape"],scale=Ricker(x,MLE["a"],MLE["b"])/MLE["shape"]),0.01,10,add=T,col="green",lty=2)
  return(MLE)
}


MyxRicker(params=c(a=2,b=0.2,shape=30),data=myxdat[,-1])   # test the function


MyxRicker_ci <- function(LikFunc = NLL_myxRicker, params=c(a=2,b=0.2,shape=30), params_selected = c("a","b"), data=myxdat[,-1], param1_lims= c(0.1,10),param2_lims=c(0.01,1)){
  temp <- optim(LikFunc,par = params,data=data)
  MLE <- temp$par
  ML <- temp$value
  notselected <- setdiff(names(MLE),params_selected)
  notselected_ndx <- which(names(MLE)==notselected)
  param_ndx <- match(params_selected,names(MLE))
  
  var1 <- seq(param1_lims[1],param1_lims[2],length=100)
  var2 <- seq(param2_lims[1],param2_lims[2],length=100)
  
  pars <- params
  pars[notselected_ndx] <- MLE[notselected]
  likarr <- matrix(NA,nrow=length(var1),ncol=length(var2))
  v1=1;v2=1
  for(v1 in 1:length(var1)){
    pars[param_ndx[1]] <- var1[v1]
    for(v2 in 1:length(var2)){
      pars[param_ndx[2]] <- var2[v2]
      likarr[v1,v2] <- -NLL_myxRicker(params=pars,data=data)
    }
  }
  
  image(x=var1,y=var2,z=likarr,zlim=c(-130,-29),col=topo.colors(12))
  contour(x=var1,y=var2,z=likarr,levels=(-ML-c(2,4,6,8,10)),add=TRUE,lwd=1,col=gray(0.3))
  
  ci1 <- range(var1[apply(likarr,1,max)>(-ML-2)])
  ci2 <- range(var2[apply(likarr,2,max)>(-ML-2)])
  
  out<- data.frame(var1 = ci1,var2=ci2)
  return(out)
}


#  test: Ricker params "a" and "b"
testab <- MyxRicker_ci(LikFunc = NLL_myxRicker, params=c(a=2,b=0.2,shape=30), params_selected = c("a","b"), data=myxdat[,-1], param1_lims= c(0.1,9),param2_lims=c(0.01,0.5))

# test: Ricker param "a" and gamma "shape" 
testas <- MyxRicker_ci(LikFunc = NLL_myxRicker, params=c(a=2,b=0.2,shape=30), params_selected = c("a","shape"), data=myxdat[,-1], param1_lims= c(1,8),param2_lims=c(10,150))


testab
testas

