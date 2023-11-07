
# load packages ---------------------

library(glmmTMB)
library(DHARMa)

# load slugs data -------------------------

slugs<-read.table( 'http://www.bio.ic.ac.uk/research/mjcraw/statcomp/data/slugsurvey.txt', header=TRUE)
head(slugs)

# write.csv(slugs,file = "slugs.csv", row.names = F)


# Summarize data ---------------------

out <- table(slugs$slugs,slugs$field)
out


barplot(t(out), beside=TRUE, angle=c(45,135), density=c(20,20), col=c('black','red'), legend.text=TRUE, xlab='# of slugs', ylab='frequency')


coords<-barplot(t(out), beside=TRUE, angle=c(45,135), density=c(20,20), ylim=c(0,27), col=c('black','red'), xlab='# of slugs', ylab='frequency')
box()
legend(coords[1,8], 26, c('nursery','rookery'), density=c(20,20), angle=c(45,135), fill=c('black','red'), cex=c(.8,.8),bty='n')


# poisson likelihood -----------------------------

poi.1<-function(data,p) -sum(log(dpois(data$slugs,lambda=p)))


mean(slugs$slugs)


out1 <- nlm(function(p) poi.1(slugs,p),2)
out1


as.numeric(as.factor(slugs$field))


# poisson, mean varies with field ----------------------------

poi.2<-function(data,p) {
  field.dummy<-as.numeric(as.factor(slugs$field))-1
  mylambda<-p[1]+p[2]*field.dummy
  negloglike<- -sum(dpois(data$slugs,lambda=mylambda,log=T))
  return(negloglike)
}


tapply(slugs$slugs,slugs$field,mean)

out2 <- nlm(function(p) poi.2(slugs,p),c(1.2,1))
out2


out1$minimum
out2$minimum

# compare with canned R functions ------------------

mod1 <- glm(slugs~1, data=slugs, family=poisson())
mod2 <- glm(slugs~field, data=slugs, family=poisson())

logLik(mod1)
logLik(mod2)


# AIC function ----------------

my.aic<-function(output) -2*(-output$minimum) + 2*length(output$estimate)


my.aic(out1)
my.aic(out2)

AIC(mod1,mod2)   # compare with models fitted with 'glm()'



#common lambda and theta
  
zip1<-function(data,p) {
  lambda<-p[1]
  theta<-p[2]
  zero.term<-sum(log(theta+(1-theta)* dpois(data$slugs[data$slugs==0], lambda)))   # stochastic process #1
  nonzero.term<-sum(log((1- theta)* dpois(data$slugs[data$slugs>0], lambda)))      # stochastic process #2
  negloglike<- -(zero.term+nonzero.term)
  negloglike
}


mean(slugs$slugs[slugs$slugs>0])

table(slugs$slugs)[1]/sum(table(slugs$slugs))


out7 <- nlm(function(p) zip1(slugs,p),c(3,.4))
out7

## compare with ZIP model using TMB

mod7 <- glmmTMB(slugs~1,slugs,family=poisson(),ziformula = ~1)
logLik(mod7)



my.aic(out7)

AIC(mod1,mod2,mod7)


# different lambda, same theta 
zip2<-function(data,p) {
  field.dummy<-as.numeric(as.factor(slugs$field))-1
  mylambda<-p[1]+p[3]*field.dummy
  theta<-p[2]
  zero.term<-sum(ifelse(data$slugs==0,log(theta+(1-theta)* dpois(data$slugs,lambda=mylambda)),0))
  nonzero.term<-sum(ifelse(data$slugs>0,log((1-theta)* dpois(data$slugs,lambda=mylambda)),0))
  negloglike<- -(zero.term+nonzero.term)
  negloglike
}


negloglike<- -sum(ifelse(data$slugs==0,log(theta+(1-theta)* dpois(data$slugs,lambda=mylambda)),log((1-theta)* dpois(data$slugs,lambda=mylambda))))


zip2.alt<-function(data,p) {
  field.dummy<-as.numeric(as.factor(slugs$field))-1
  mylambda<-p[1]+p[3]*field.dummy
  theta<-p[2]
  negloglike<- -sum(ifelse(data$slugs==0,log(theta+(1-theta)* dpois(data$slugs,lambda=mylambda)),log((1-theta)*dpois(data$slugs,lambda=mylambda))))
  negloglike
  
}

tapply(slugs$slugs[slugs$slugs>0],slugs$field[slugs$slugs>0],mean)


out8 <- nlm(function(p) zip2(slugs,p),c(3.4,.42,-.4))
out8


#different theta
zip3<-function(data,p){
  field.dummy<-as.numeric(as.factor(data$field))-1
  mylambda<-p[1]
  theta<-p[2]+p[3]*field.dummy
  zero.term<-sum(ifelse(data$slugs==0,log(theta+(1-theta)* dpois(data$slugs,lambda=mylambda)),0))
  nonzero.term<-sum(ifelse(data$slugs>0,log((1-theta)* dpois(data$slugs,lambda=mylambda)),0))
  negloglike<- -(zero.term+nonzero.term)
  negloglike
}


table(slugs$slugs,slugs$field)[1,]/40


out9 <- nlm(function(p) zip3(slugs,p),c(3,.6,-.4))
out9


norm.neglike<-function(data,p) {
  t.y<-log(data$slugs+1)
  mu<-p[1]
  my.sd<-p[2] 
  negloglike<- -sum(log(dnorm(t.y,mean=mu,sd=my.sd)))
  negloglike
}


mean(log(slugs$slugs+1))

sd(log(slugs$slugs+1))


out.norm <- nlm(function(p) norm.neglike(slugs,p),c(.73,.74))
out.norm


#calculate negative loglikelihood for AIC
norm.like<-function(data,out) {
  t.y<-log(data$slugs+1)
  mu<-out$estimate[1]
  my.sd<-out$estimate[2]
  negloglike<- -sum(log(dnorm(t.y,mean=mu, sd=my.sd)*1/(data$slugs+1)))
  out<-list(negloglike,out$estimate)
  names(out)<-c("minimum","estimate")
  out
}


out20 <- norm.like(slugs,out.norm)
out20


tapply(log(slugs$slugs+1),slugs$field,mean)


norm.neglike2<-function(data,p) {
  t.y<-log(data$slugs+1)
  field.dummy<-as.numeric(data$field)-1
  mu<-p[1]+field.dummy*p[3]
  my.sd<-p[2]
  negloglike<- -sum(log(dnorm(t.y,mean=mu,sd=my.sd)))
  negloglike
}


outnorm2 <- nlm(function(p) norm.neglike2(slugs,p),c(.5,.7,.5))
outnorm2


norm.like2<-function(data,out) {
  t.y<-log(data$slugs+1)
  field.dummy<-as.numeric(data$field)-1
  mu<-out$estimate[1]+field.dummy*out$estimate[3]
  my.sd<-out$estimate[2]
  negloglike<- -sum(log(dnorm(t.y,mean=mu, 
    sd=my.sd)*1/(data$slugs+1)))
  out<-list(negloglike,out$estimate)
  names(out)<-c("minimum","estimate")
  out
}


out21 <- norm.like2(slugs,outnorm2)
my.aic(out21)
my.aic(out9)


model.names<-c('Pois.common','Pois.mean','Zip.common', 'Zip.mean','Zip.theta','lognormal','lognormal.mean')


models<-list(out1,out2,out7,out8,out9,out20,out21)


AIC.func<-function(model.list,n,modelnames) {
  output<-NULL
  for (i in 1:length(model.list)) {
    cur.model<-model.list[[i]]
    LL<- -cur.model$minimum
    K<-length(cur.model$estimate)
    AIC<- -2*LL + 2*K
    AICc<-AIC + 2*K*(K+1)/(n-K-1)
    output<-rbind(output,c(LL,K,AIC,AICc))
  }
  colnames(output)<-c('LogL','K','AIC','AICc')
  minAICc<-min(output[,"AICc"])
  deltai<-output[,"AICc"]-minAICc
  rel.like<-exp(-deltai/2)
  wi<-round(rel.like/sum(rel.like),3)
  out<-data.frame(modelnames,output,deltai,wi)[order(deltai,decreasing = F),]
  out
}

AIC.func(models,80,model.names)
