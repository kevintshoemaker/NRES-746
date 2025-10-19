
#  NRES 746, Lab 3                              
#  University of Nevada, Reno                      
#  DIY Likelihood Functions       ---------------------                


# Exercise 3.1a  ---------------------

NLL_frogOccupancy <- function(p){
  -sum(dbinom(c(3,2,6),prob=p,size = 10,log=T))
}


NLL_frogOccupancy(p=0.5)   # test your function


# 3.1b -----------------

xvals <- seq(0.001,0.999,0.001)
nlls <- sapply(1:length(xvals),function(t) NLL_frogOccupancy(xvals[t]) ) 
plot(nlls~xvals,xlab="parameter \"p\"",ylab="neg log lik",type="l")

MLE <- xvals[which.min(nlls)]
ML <- min(nlls)

CI95 <- range(xvals[nlls<=(ML+2)])

abline(v=MLE,col="green",lwd=2)
abline(h=ML+2,col="blue")
abline(v=CI95,col="green",lty=2)



# 3.2a ------------------

rffr <- emdbook::ReedfrogFuncresp     # from Bolker's "emdbook" package
  # ?Reedfrog      # learn more about this dataset
head(rffr)


## 3.2a ---------------

# define a Holling type II functional response, with an initial guess about parameter values

Holl2<-function(x, a, h){(a*x)/(1+(a*h*x))}
plot(rffr$Killed~rffr$Initial, xlab="Initial density",ylab="# eaten")
curve(Holl2(x, a=0.5, h=1/80), add=TRUE,col="red")


# Write a likelihood function

#    params: vector of params to estimate (a and h from the Holling type II functional response)

binomNLL2<-function(params){
  N = rffr$Initial; k=rffr$Killed   # hard code the data in
	-sum(dbinom(k,Holl2(N,params["a"],params["h"])/N,size=N,log=TRUE))
}

binomNLL2(c(a=0.4,h=1/200))

opt2 <- optim(c(a=0.5,h=(1/80)), binomNLL2)  #use default simplex algorithm
MLE = opt2$par
MaxLik = opt2$value


# 3_2b ----------

Rffuncresp <- function(params,dat=rffr){
  df = data.frame(
    x = (min(dat$Initial)-2):(max(dat$Initial)+2)
  )
  df$y = Holl2(df[,1], a=params["a"], h=params["h"])
  df$lwr = qbinom(0.025,df$x,(df$y/df$x))
  df$upr = qbinom(0.975,df$x,(df$y/df$x))
  p = ggplot(df,aes(x=x,y=y)) +
    geom_ribbon(aes(ymin=lwr,ymax=upr),fill=gray(.6)) +
    geom_path(col="darkgreen",lwd=2) +
    geom_point(data=dat,aes(x=Initial,y=Killed),cex=2,col="darkblue")  +
    theme_classic()
  print(p)
  return(df)
}


Rffuncresp(params=MLE,dat=rffr)


# Exercise 3.3a -----------

library(emdbook)
data(MyxoTiter_sum)      # load the data
# head(MyxoTiter_sum)   
myxdat <- subset(MyxoTiter_sum, grade==1)    # select just the most virulent strain

plot(myxdat$titer~myxdat$day,xlim=c(0,10))    # visualize the relationship


# 3.3a --------

Ricker <- function(x,a,b){
  a*x*exp(-b*x)
}

NLL_myxRicker <- function(params){
  mn <- Ricker(myxdat$day,params["a"],params["b"])
  sh <- mn * params["rate"]
  -sum(dgamma(myxdat$titer,shape=sh,rate=params["rate"],log=TRUE))
}


NLL_myxRicker(params=c(a=4,b=0.2,rate=2))   # test the function


# 3.3b --------------

predict_myxRicker1 <- function(mle){
  df <- data.frame(
    x=seq(min(myxdat$day)-1,max(myxdat$day)+1,length=100)
  )
  df$mean <- Ricker(df$x,mle["a"],mle["b"])
  df$lwr <- qgamma(0.025,df$mean*mle["rate"],mle["rate"] )
  df$upr <- qgamma(0.975,df$mean*mle["rate"],mle["rate"] )
  g = ggplot(df,aes(x,mean)) +
    geom_ribbon(aes(ymin=lwr,ymax=upr),fill="gray") +
    geom_path(col="darkgreen",lwd=2) +
    geom_point(data=myxdat,aes(x=day,y=titer),size=2) +
    labs(x="Days since infection",y="Virus titer") +
    theme_classic()
  print(g)
  return(df)
}


predict_myxRicker1(mle=MLE)   # test the function


# 3.4a -----------

CI_myxRicker1 <- function(mle, H){
  stderr = sqrt(diag(solve(H)))
  moe = stderr*qnorm(0.025,lower=F)
  cbind(mle=mle,se=stderr,lwr=mle-moe,upr=mle+moe)
}


CI_myxRicker1(MLE,H)   # test the function


# 3.4b --------------

ci_fun = function(par, d) par[1]*d * exp(-par[2]*d)
  
pi_fun = function(par, d){
  (par[1]*d * exp(-par[2]*d)) * par[3]
}

predict_myxRicker2 <- function(mle, H, prediction=F){
  VCV = solve(H)
  df <- data.frame(
    x=seq(min(myxdat$day)-1,max(myxdat$day)+1,length=100)
  )
  df$mean <- Ricker(df$x,mle["a"],mle["b"])
  if(prediction){
    df$upr = df$lwr = NA 
    for(i in 1:nrow(df)){
      g <- numDeriv::grad(pi_fun, mle, d=df$x[i])
      se_shap <- sqrt(t(g) %*% VCV %*% g)
      df[i,c("lwr","upr")] <- qgamma(c(0.025,0.975),df$mean[i]*mle["rate"],mle["rate"] )
    }
  }else{
    se <- sapply(df$x,function(t){
      g <- numDeriv::grad(ci_fun, mle[1:2], d=t)
      sqrt(t(g) %*% VCV[1:2,1:2] %*% g)[1,1]  
    })    # fine to use a for loop here!
    df$lwr <- df$mean - se * qnorm(0.025,lower=F)
    df$upr <- df$mean + se * qnorm(0.025,lower=F)
  }
  g = ggplot(df,aes(x,mean)) +
    geom_ribbon(aes(ymin=lwr,ymax=upr),fill="gray") +
    geom_path(col="darkgreen",lwd=2) +
    geom_point(data=myxdat,aes(x=day,y=titer),size=2) +
    labs(x="Days since infection",y="Virus titer") +
    theme_classic()
  print(g)
  return(df)
}


predict_myxRicker2(MLE,H,prediction = T)   # test the function
predict_myxRicker2(MLE,H,prediction = F)


# 3_4c ------------

CI_myxRicker2 <- function(mle, H, nll, profile=T,alpha=0.05){
  if(profile){
    parnames <- names(mle)
    npar = length(mle)
    out = NULL
    for(i in 1:npar){
      thisparname = parnames[i]
      otherparnames = setdiff(parnames,thisparname)
      parvals = seq(mle[thisparname]-0.25*mle[thisparname],mle[thisparname]+0.25*mle[thisparname],length=500)
      thisNLL = function(params,thispar){
        params = c(params,thispar)
        mn <- Ricker(myxdat$day,params["a"],params["b"]); sh <- mn * params["rate"]
       -sum(dgamma(myxdat$titer,shape=sh,rate=params["rate"],log=TRUE))
      }
      thispar= mle[thisparname] 
      thisprofile = suppressWarnings( sapply(parvals,function(t){ thispar[]=t; optim(mle[otherparnames], thisNLL, thispar=thispar )$value } ) )
      plausible = parvals[2*(thisprofile-nll) < qchisq(1-alpha,1)]
      out = rbind(out, c(mle=unname(mle[thisparname]),lb=min(plausible), ub=max(plausible) )  )
    }
    return(out)
  }else{
    stderr = sqrt(diag(solve(H)))
    moe = stderr*qnorm(alpha/2,lower=F)
    return(cbind(mle=mle,lwr=mle-moe,upr=mle+moe))
  }
}


CI_myxRicker2(MLE,H,MinNLL,profile=T,alpha=0.9)

CI_myxRicker2(MLE,H,MinNLL,profile=F,alpha=0.9)

