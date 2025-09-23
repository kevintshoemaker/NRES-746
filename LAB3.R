
#  NRES 746, Lab 3                              
#  University of Nevada, Reno                      
#  DIY Likelihood Functions       ---------------------                


# Reed frog example ----------------

rfp <- emdbook::ReedfrogPred  # load the data using the 'emdbook' package:
head(rfp)


# Take a subset of the data

rfp_sub <- subset(rfp, (rfp$pred=='pred')&(rfp$size=="small")&(rfp$density==10))
rfp_sub


rfp_sub$killed <- with(rfp_sub, density-surv)  # make column for number killed in each replicate trial
with(rfp_sub, sum(dbinom(killed, 10, prob=0.5, log=TRUE)) )    # expression of data likelihood(log scale)


L = dbinom(rfp_sub$killed,size=10,prob=0.5)  # evaluate data likelihood with p=0.5
L

prod(L)    # joint data likelihood


p_seq <- seq(0.01, 1, length=100)     # prepare for visualizing the likelihood across parameter space


Lik <- sapply(p_seq, function(t) prod(dbinom(rfp_sub$killed,10,prob=t)) )

plot(Lik~p_seq,lty="solid",type="l", xlab="Predation Probability", ylab="Likelihood")


# plot out the log-likelihood

LogLik <- sapply(p_seq, function(t) sum(dbinom(rfp_sub$killed,10,prob=t,log=T)) )
plot(LogLik~p_seq,lty="solid",type="l", xlab="Predation Probability", ylab="Log Likelihood")


p_seq[which.max(LogLik)]     # MLE for probability of predation


plot(LogLik~p_seq,lty="solid",type="l", xlab="Predation Probability", ylab="Log Likelihood")
abline(v=0.25,lwd=3)


# Write a likelihood function

binomNLL1 <- function(p) {
  -sum(dbinom(rfp_sub$killed, size=10, prob=p, log=TRUE))
}


# use "optim()" to find the MLE

opt1 <- optim(fn=binomNLL1, par = c(p=0.5), method = "BFGS")   # use "optim()" to estimate the parameter value that maximizes the likelihood function 


opt1    # check out the results of "optim()"

MLE = opt1$par
MaxLik = opt1$value


opt1$convergence


hist(rfp_sub$killed,xlim=c(0,10),freq=F,main="",xlab="Number of tadpoles killed")
curve(dbinom(x,prob=MLE,size=10),add=T,from=0,to=10,n=11,lwd=2, col="darkgreen")


NLL_frogOccupancy(p=0.5)   # test your function


# 3.2a ------------------

rffr <- emdbook::ReedfrogFuncresp     # from Bolker's "emdbook" package
  # ?Reedfrog      # learn more about this dataset
head(rffr)


hist(rffr$Killed/rffr$Initial)

binomNLL2(c(a=0.4,h=1/200),c(5,10,15),c(3,5,6))

opt2 <- optim(c(a=0.5,h=(1/80)), binomNLL2, N=rffr$Initial, k=rffr$Killed)  #use default simplex algorithm
MLE = opt2$par
MaxLik = opt2$value


plot(rffr$Killed~rffr$Initial, xlab="Initial density",ylab="# eaten")
curve(Holl2(x, a=MLE["a"], h=MLE["h"]), add=TRUE,col="red")


Rffuncresp(params=MLE,dat=rffr)


# how to generate 'plug-in' prediction intervals

xvec <- seq(40,80,5)
yvec <- 0.5/(1+0.5*0.015*xvec) * xvec
upper<-qbinom(0.975,prob=yvec/xvec, size=xvec) 
lower<-qbinom(0.025,prob=yvec/xvec, size=xvec)

upper
lower


# Exercise 3.3a -----------

library(emdbook)
data(MyxoTiter_sum)      # load the data
# head(MyxoTiter_sum)   
myxdat <- subset(MyxoTiter_sum, grade==1)    # select just the most virulent strain

plot(myxdat$titer~myxdat$day,xlim=c(0,10))    # visualize the relationship


NLL_myxRicker(params=c(a=4,b=0.2,rate=2))   # test the function

temp <- optim(NLL_myxRicker,par = c(a=2,b=0.2,rate=2))
MLE = temp$par
MinNLL = temp$value 
MLE


# Plug-in prediction intervals!

upper<-qgamma(0.975,shape=?, rate=?)   # remember that the mean of the gamma distribution is shape*scale
lower<-qgamma(0.025,shape=?, rate=?)


predict_myxRicker1(mle=MLE)   # test the function

temp <- optim(NLL_myxRicker,par = c(a=2,b=0.2,rate=2),hessian = T)
MLE = temp$par
MinNLL = temp$value 
H = temp$hessian
H

CI_myxRicker1(MLE,H)   # test the function


# assume X is a random variable with mean of c(2.2, 4.1)

exp_X = c(2.2,4.1)   # expected value of random MVN variable
vcv_X = matrix(c(1,.2,.2,1.5),nrow=2)  # vcv matrix of X

# imagine we want to evaluate a transformation of X: 1/(x1^2 + log(x2))

myfun = function(x) 1/(x[1]^2+log(x[2]))   # write a function for your transformation

# compute the gradient of the transformation with respect to the components of X,
      # evaluated at the mean values of X
f_prime_X = numDeriv::grad(myfun,exp_X)

# apply the delta method to obtain an approximate standard error:
var_fX = (t(f_prime_X) %*% vcv_X %*% f_prime_X)[1,1] 

# compute the standard error:
se_fX = sqrt(var_fX)

# compute the critical value for 95% CI
z_crit = qnorm(0.025,lower=F)

# compute confidence interval
ci_fX = myfun(exp_X) + z_crit*se_fX * c(mean=0,lower=-1,upper=1)
ci_fX


predict_myxRicker2(MLE,H,prediction = T)   # test the function
predict_myxRicker2(MLE,H,prediction = F)


CI_myxRicker2(MLE,H,MinNLL,profile=T,alpha=0.9)

CI_myxRicker2(MLE,H,MinNLL,profile=F,alpha=0.9)

