# Splines demo for NRES 746

## (examples modified from Simon Wood's GAM book)

library(gamair); data(engine)
?engine
engine = engine[order(engine$size),]
size=engine$size; wear=engine$wear   # safe "attach"

# explore polynomial basis expansion --------

xseq = seq(min(size),max(size),length=100)
X = cbind(rep(1,100), poly(xseq,5) )    # polynomial basis expansion and intercept

## visualize poly basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)

## weighted polynomial basis expansion -----

w = c(3,-1.5,0.6,-0.5,-0.2,-0.1)

lp = X %*% w

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
lines(xseq,X[,1]*w[1],col=cols[1],lwd=2)
lines(xseq,X[,2]*w[2],col=cols[2],lwd=2)
lines(xseq,X[,3]*w[3],col=cols[3],lwd=1)
lines(xseq,X[,4]*w[4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5]*w[5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6]*w[6],col=cols[6],lwd=1,lty=4)

lines(xseq,lp[,1],lwd=2,col="darkgreen")  # 

## use lm to find optimal weights ------

Xo = cbind(rep(1,length(size)), poly(size,5) )    # polynomial basis expansion and intercept
m = lm(wear~0+Xo)
w_fit = m$coefficients

lp = X %*% w_fit

lines(xseq,lp[,1],lwd=2,col="black")  # 

## explore other possible functional forms --------

ws = replicate(10,(X%*%runif(ncol(X),-2,2))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


# explore 'tent' basis ----------
   ## also known as piecewise linear regression

j=3
tf <- function(x,xj,j){  # xj is knots, x is data covariate, j is jth knot
  dj=xj*0;dj[j]=1
  approx(xj,dj,x)$y   # jth tent basis function
}

tf.X <- function(x,xj){
  nk = length(xj); n=length(x)
  X = matrix(NA,n,nk)
  for(j in 1:nk) X[,j] <- tf(x,xj,j)
  X
}

sj=seq(min(size),max(size),length=6)  # knots

X = tf.X(xseq,sj)  # design matrix

## visualize tent basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)

## visualize range of functional forms -----------

ws = replicate(10,(X%*%runif(ncol(X),-3,3))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


## find optimal weights using lm -------

Xo = tf.X(size,sj)
b = lm(wear~0+Xo)  # piecewise linear regression

plot(size,wear)
lines(xseq,X %*% coef(b),lwd=2)


## penalized smooth fitting --------

y=wear;x=size;sp=2  # sp is complexity penalty
prs.fit <- function(y,x,xj,sp){
  X=tf.X(x,xj)  # model matrix
  D = diff(diag(length(xj)),differences=2)
  X=rbind(X,sqrt(sp)*D)
  y=c(y,rep(0,nrow(D)))
  lm(y~X-1)
}

sj=seq(min(size),max(size),length=20)  # note: now we have 20 basis functions
b=prs.fit(wear,size,sj,5)   ## try different smoothing penalties
plot(size,wear)
Xp = tf.X(xseq,sj)
lines(xseq,Xp%*% coef(b))


# cubic regression splines -------

library(splines)

X = splines::bs(xseq,df=6,intercept=T)   # design matrix

## visualize cr basis -------

plot(size,wear,xlim=c(1.3,3.1),ylim=c(-1,5))
cols=rainbow(6)
lines(xseq,X[,1],col=cols[1],lwd=2)
lines(xseq,X[,2],col=cols[2],lwd=2)
lines(xseq,X[,3],col=cols[3],lwd=1)
lines(xseq,X[,4],col=cols[4],lwd=1,lty=2)
lines(xseq,X[,5],col=cols[5],lwd=1,lty=3)
lines(xseq,X[,6],col=cols[6],lwd=1,lty=4)


## visualize range of functional forms -----------

ws = replicate(10,(X%*%runif(ncol(X),-3,3))[,1] )

plot(1,1,xlim=c(1.3,3.1),ylim=c(-5,6),pty="n")
sapply(1:ncol(ws), function(t) lines(xseq,ws[,t],col=sample(cols,1))  )


## find optimal weights using lm -------

Xo = bs(size,df=6,intercept = T)
b = lm(wear~0+Xo)  # piecewise linear regression

plot(size,wear)
lines(xseq,X %*% coef(b),lwd=2)


## penalized smooth fitting --------

mod = gam(wear~s(size,bs="cr",k=8,sp=5),method="REML")
b=coef(mod)   ## try different smoothing penalties
plot(size,wear)
Xp= predict(mod,newdata = data.frame(size=xseq), type = "lpmatrix")
lines(xseq,Xp%*%b, lwd=2)


























