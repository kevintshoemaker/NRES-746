############


dbeta(0.99999999,1,1)

curve(dbeta(x,10,10),0,1)

dbeta(.6,5,5)

qbeta(1,5,5)

pbeta(.7,5,5)


curve(pbeta(x,5,5),-1,2)

qnorm(.975,0,1)

curve(qnorm(x),0,1,ylim=c(-10,10))


hist(rnorm(1000),freq = F)
curve(dnorm(x),add=T)


hist(rpois(100,1.5),breaks=seq(-0.5,10.5,1),freq=F,ylim=c(0,.5))
points(0:20,dpois(0:20,1.5),col="green")

hist(rnbinom(10000,mu=4,size=1),freq=F,breaks=seq(-0.5,90.5,1),ylim=c(0,.3))
lines(0:50,dnbinom(0:50,mu=4,size=1),col="green")

hist(rpois(1000,4),freq=F,breaks=seq(-0.5,12.5,1),add=T,border ="red",col=NA)
lines(0:13,dpois(0:13,4),col="red")


xvals=round(sample(runif(n,0,10)),2)
simlinreg <- function(xvals=xvals,n=10,b0 = 1, b1 = 0.5, sig = 1){
  yexp <- b0+b1*xvals
  yvals <- rnorm(n,yexp,sig)
  data.frame(x=xvals,y=yvals)
}

simlinreg(xvals)

reps <- lapply(1:100, function(t) simlinreg(xvals) )

mods <- lapply(reps,function(t) lm(y~x,t) )

coefs <- lapply(mods,coef)


