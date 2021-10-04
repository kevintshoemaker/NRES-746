
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



###### Read in the reed frog data set

#rfp <- read.csv("ReedfrogPred.csv")

########
# alternatively, load the data using the 'emdbook' package:

library(emdbook)
rfp <- ReedfrogPred

head(rfp)


##### Take a subset of the data

rfp_sub <- subset(rfp, (rfp$pred=='pred')&(rfp$size=="small")&(rfp$density==10))
rfp_sub


killed <- rfp_sub$density-rfp_sub$surv
N=rfp_sub$density
p=0.5
sum(dbinom(killed, size=N, prob=p, log=TRUE))    # expression of data likelihoodb(log scale)


num_killed <- rfp_sub$density-rfp_sub$surv     # specify vector of "successes" (being eaten!)
num_killed


dbinom(num_killed,size=10,prob=0.5)  # evaluate data likelihood with p=0.5


prod(dbinom(num_killed,size=10,prob=0.5))    # joint data likelihood


p <- seq(0.01, 1, length=100)     # prepare for visualizing the likelihood across parameter space


Lik <- numeric(length=100)


#########
# plot out the likelihood

for(i in 1:100){
  Lik[i] <- prod(dbinom(num_killed,size=10,prob=p[i]))
}
plot(Lik~p,lty="solid",type="l", xlab="Predation Probability", ylab="Likelihood")



########
# plot out the log-likelihood

p <- seq(0.01, 1, by=0.01)
LogLik <- numeric(length=100)
for(i in 1:100){
  LogLik[i] <- sum(dbinom(num_killed, size=10, 
  prob=p[i],log=TRUE))
}
plot(LogLik~p,lty="solid",type="l", xlab="Predation Probability", ylab="Log Likelihood")


p[which(LogLik==max(LogLik))]     # MLE for probability of predation


plot(LogLik~p,lty="solid",type="l", xlab="Predation Probability", ylab="Log Likelihood")
abline(v=0.25,lwd=3)


###########
# Write a likelihood function

#    p: probability of predation per trial (param to estimate)
#    k: number killed per trial   (data)
#    N: number of tadpoles per trial (data)

binomNLL1 <- function(p, k, N) {
  -sum(dbinom(k, size=N, prob=p, log=TRUE))
}


#####
# use "optim()" to find the MLE

opt1 <- optim(fn=binomNLL1, par = c(p=0.5), N = 10, k = num_killed, method = "BFGS")   # use "optim()" to estimate the parameter value that maximizes the likelihood function 


opt1    # check out the results of "optim()"


opt1$convergence


opt1$par  # MLE


opt1$value     # max. likelihood (actually minimum negative-log-likelihood)


exp(-opt1$value)   # convert to likelihood


hist(num_killed,xlim=c(0,10),freq=F)
curve(dbinom(x,prob=opt1$par,size=10),add=T,from=0,to=10,n=11)

