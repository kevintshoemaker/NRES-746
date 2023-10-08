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

reps <- lapply(1:1000, function(t) simlinreg(xvals) )

mods <- lapply(reps,function(t) lm(y~x,t) )

coefs <- lapply(mods,coef)

plot(reps[[1]])
lapply(coefs,abline,col=rgb(0.4,0.4,0.4,0.4))


######

data(mtcars)

df = mtcars[,c("disp","mpg")]

plot(mpg~disp,df)



negexp <- function(a,b,x){
  a*exp(b*x)
}

curve(negexp(35,-0.002,x),add=T)

params <- c(a=35,b=-0.002,sigma=2)
nll_full <- function(params){
  -sum(dnorm(df$mpg,negexp(params["a"],params["b"],df$disp),params["sigma"],log=T))
}
nll_full(params)

mle <- optim(params,nll_full)
bestparams <- mle$par
max_nll <- -mle$value


## develop profile likelihood confidence intervals


nll_fixsig <- function(params,sigma){
  -sum(dnorm(df$mpg,negexp(params["a"],params["b"],df$disp),sigma,log=T))
}

nll_fixa <- function(params,a){
  -sum(dnorm(df$mpg,negexp(a,params["b"],df$disp),params['sigma'],log=T))
}

nll_fixb <- function(params,b){
  -sum(dnorm(df$mpg,negexp(params["a"],b,df$disp),params['sigma'],log=T))
}

as <- seq(25,40,length=1000)
bs <- seq(-0.01,-0.0001,length=1000)
sigs <- seq(1,10,length=1000)

mls_a = sapply(as,function(t) -1*optim(params[setdiff(names(params),"a")],nll_fixa,a=t)$value )
ci_a = range(as[which(abs(mls_a-max_nll)<(qchisq(0.95,1)/2) )])
plot(as,mls_a,type="l")
abline(v=ci_a,lty=2)

mls_b = sapply(bs,function(t) -1*optim(params[setdiff(names(params),"b")],nll_fixb,b=t)$value )
ci_b = range(bs[which(abs(mls_b-max_nll)<(qchisq(0.95,1)/2) )])
plot(bs,mls_b,type="l")
abline(v=ci_b,lty=2)

mls_sig = sapply(sigs,function(t) -1*optim(params[setdiff(names(params),"sigma")],nll_fixsig,sigma=t)$value )
ci_sig = range(sigs[which(abs(mls_sig-max_nll)<(qchisq(0.95,1)/2) )])
plot(sigs,mls_sig,type="l")
abline(v=ci_sig,lty=2)



