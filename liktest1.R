
rm(list=ls())

## load data ---------

y = c(2,5,3,0,4,4)   # pika haypiles at 6 talus sites


## write likelihood function --------

    # assume poisson response with constant mean/expectation

lik = function(l) prod(dpois(y,l))     # likelihood function

lik(1)   # test the function

nll = function(l) -sum(dpois(y,l,log=T))  # negative log likelihood function.

nll(5)


## visualize likelihood function --------

lams = seq(0.1,10,length=100)

liks = sapply(lams,lik)

plot(lams,liks,type="l")


## visualize negative log likelihood function  ---------

nlls = sapply(lams,nll)

plot(lams,nlls,type="l")


## find minimum by brute force  ----------

lams[which.max(liks)]   # argmax: value of lambda that maximizes the likelihood!

lams[which.min(nlls)]   # argmin: value of lambda that minimizes the likelihood



## find minimum by numerical optimization  ------------

opt = optim(par=1,fn=nll,method="Brent",lower=0.1,upper=20)

mle = opt$par   # value of lambda that minimizes the objective function (nll)



plot(lams,nlls,type="b")
abline (v=mle,lwd=2,color="darkgreen")



