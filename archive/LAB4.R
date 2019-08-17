
############################################################
####                                                    ####  
####  NRES 746, Lab 4                                   ####
####                                                    ####
####  Kevin Shoemaker                                   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Bayesian inference                                ####
############################################################



########
# 1a


library(HDInterval)
library(emdbook)

Ricker <- function(x,a,b){
  a*x*exp(-b*x)
}


filename <- "BUGSmodel_ricker.txt"
cat("
    model {
      
      #############
      # LIKELIHOOD
      ############
      for(obs in 1:n.observations){
        titer[obs] ~ dgamma(shape,rate[obs])
        rate[obs] <- shape/exp.titer[obs]
        exp.titer[obs] <- a*days[obs]*exp(-1*b*days[obs])
      }
      
      #############
      # PRIORS
      ############
      shape ~ dgamma(0.01,0.01)
      a ~ dgamma(0.01,0.01)
      b ~ dgamma(0.01,0.01)
    }
  ",file=filename
)



params.to.store <- c("shape","a","b")    # specify the parameters we want to get the posteriors for

#jags.fit_1a <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=50000,model.file="BUGSmodel_ricker.txt",n.chains = 3,n.burnin = 5000,n.thin = 20 )

#jags.fit_1a   # test the model fit


########
# package into function

Myx_Ricker_JAGS <- function(data=Myx[,-1],JAGScode="BUGSmodel_ricker.txt",CImethod="HPD",CIlevel=0.95){
  myx.data.for.bugs <- list(
    titer = data[,2],
    days =data[,1],
    n.observations = nrow(data)
  )
  
  myx.data.for.bugs
  
  init.vals.for.bugs <- function(){
    list(
      shape=runif(1,20,30),
      a=runif(1,0.05,0.3),
      b=runif(1,1,1.3)
    )
  }
  init.vals.for.bugs()
  
  params.to.store <- c("shape","a","b")    # specify the parameters we want to get the posteriors for
  
  jags.fit_1a <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=50000,model.file=JAGScode,n.chains = 3,n.burnin = 5000,n.thin = 20 )
  jags.fit_1a_coda <- as.mcmc(jags.fit_1a)
  #summary(jagsfit.mcmc.1a)
  plot(jags.fit_1a_coda)
  temp <- densityplot(jags.fit_1a_coda)
  print(temp)
  plot(data[,1],data[,2],xlab="days",ylab="titer",ylim=c(1,13),xlim=c(0,10))
  amean <- mean(jags.fit_1a$BUGSoutput$sims.list$a)
  bmean <- mean(jags.fit_1a$BUGSoutput$sims.list$b)
  shapemean <- mean(jags.fit_1a$BUGSoutput$sims.list$shape)
  curve(Ricker(x,amean,bmean),0,10,add=T,lwd=2,col="green")
  posterior_mean <- c(amean,bmean,shapemean)
  
  cintervals <- matrix(NA,nrow=2,ncol=3)
  colnames(cintervals) <- c("a","b","shape")
  if(CImethod=="HPD"){
    cintervals[,"a"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$a,credMass = CIlevel)
    cintervals[,"b"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$b,credMass = CIlevel)
    cintervals[,"shape"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$shape,credMass = CIlevel)
  }else{
    cintervals[,"a"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$a,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
    cintervals[,"b"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$b,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
    cintervals[,"shape"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$shape,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
  }
  gelmandiag <- gelman.diag(jags.fit_1a_coda)

  output = list(
    post.mean <- posterior_mean,
    conf.ints <- cintervals,
    gelman.rubin <- gelmandiag$mpsrf 
  )
  return(output)  
}



Ricker_results <- Myx_Ricker_JAGS(data=Myx[,-1],JAGScode="BUGSmodel_ricker.txt",CImethod="HPD",CIlevel=0.88)
Ricker_results
  
Ricker_results <- Myx_Ricker_JAGS(data=Myx[,-1],JAGScode="BUGSmodel_ricker.txt",CImethod="quantile",CIlevel=0.88)
Ricker_results


#####
# 2a

MMFunc <- function(x,a,b){
  (a*x)/(b+x)
}


filename <- "BUGSmodel_mm.txt"
cat("
    model {
      
      #############
      # LIKELIHOOD
      ############
      for(obs in 1:n.observations){
        titer[obs] ~ dgamma(shape,rate[obs])
        rate[obs] <- shape/exp.titer[obs]
        exp.titer[obs] <- (a*days[obs]) / (b+days[obs])
      }
      
      #############
      # PRIORS
      ############
      shape ~ dgamma(0.01,0.01)
      a ~ dgamma(0.01,0.01)
      b ~ dgamma(0.01,0.01)
    }
  ",file=filename
)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)  #Data set from grade 1 of myxo data
#head(Myx)

myx.data.for.bugs <- list(
  titer = Myx$titer,
  days = Myx$day,
  n.observations = nrow(Myx)
)

#myx.data.for.bugs

init.vals.for.bugs2 <- function(){
  list(
    shape=runif(1,20,30),
    a=runif(1,8,9),
    b=runif(1,1,1.5)
  )
}
#init.vals.for.bugs()

params.to.store <- c("shape","a","b")    # specify the parameters we want to get the posteriors for

#jags.fit_1a <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=50000,model.file="BUGSmodel_mm.txt",n.chains = 3,n.burnin = 5000,n.thin = 20 )

#jags.fit_1a   # test the model fit


########
# package into function

Myx_MM_JAGS <- function(data=Myx[,-1],JAGScode="BUGSmodel_mm.txt",CImethod="HPD",CIlevel=0.95){
  myx.data.for.bugs <- list(
    titer = data[,2],
    days =data[,1],
    n.observations = nrow(data)
  )
  
  init.vals.for.bugs2 <- function(){
    list(
      shape=runif(1,20,30),
      a=runif(1,8,9),
      b=runif(1,1,1.5)
    )
  }
  
  params.to.store <- c("shape","a","b")    # specify the parameters we want to get the posteriors for
  
  jags.fit_1a <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs2,parameters.to.save=params.to.store,n.iter=50000,model.file=JAGScode,n.chains = 3,n.burnin = 5000,n.thin = 20 )
  jags.fit_1a_coda <- as.mcmc(jags.fit_1a)
  #summary(jagsfit.mcmc.1a)
  plot(jags.fit_1a_coda)
  temp <- densityplot(jags.fit_1a_coda)
  print(temp)
  plot(data[,1],data[,2],xlab="days",ylab="titer",ylim=c(1,13),xlim=c(0,10))
  amean <- mean(jags.fit_1a$BUGSoutput$sims.list$a)
  bmean <- mean(jags.fit_1a$BUGSoutput$sims.list$b)
  shapemean <- mean(jags.fit_1a$BUGSoutput$sims.list$shape)
  curve(MMFunc(x,amean,bmean),0,10,add=T,lwd=2,col="green")
  posterior_mean <- c(amean,bmean,shapemean)
  
  cintervals <- matrix(NA,nrow=2,ncol=3)
  colnames(cintervals) <- c("a","b","shape")
  if(CImethod=="HPD"){
    cintervals[,"a"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$a,credMass = CIlevel)
    cintervals[,"b"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$b,credMass = CIlevel)
    cintervals[,"shape"] <- hdi(jags.fit_1a$BUGSoutput$sims.list$shape,credMass = CIlevel)
  }else{
    cintervals[,"a"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$a,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
    cintervals[,"b"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$b,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
    cintervals[,"shape"] <- quantile(jags.fit_1a$BUGSoutput$sims.list$shape,c((1-CIlevel)/2 ,CIlevel+((1-CIlevel)/2)))
  }
  gelmandiag <- gelman.diag(jags.fit_1a_coda)

  output = list(
    post.mean <- posterior_mean,
    conf.ints <- cintervals,
    gelman.rubin <- gelmandiag$mpsrf 
  )
  return(output)  
}



MM_results <- Myx_MM_JAGS(data=Myx[,-1],JAGScode="BUGSmodel_mm.txt",CImethod="HPD",CIlevel=0.88)
MM_results
  
MM_results <- Myx_MM_JAGS(data=Myx[,-1],JAGScode="BUGSmodel_mm.txt",CImethod="quantile",CIlevel=0.88)
MM_results


######
# 2b

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,10))
curve(MMFunc(x,a=MM_results[[1]][1],b=MM_results[[1]][2]),from=0,to=10,add=T,col="red",lwd=2)
curve(Ricker(x,a=Ricker_results[[1]][1],b=Ricker_results[[1]][2]),from=0,to=10,add=T,col="green",lty=2,lwd=2)
legend("topleft",lwd=c(2,2),lty=c(1,2),col=c("red","green"),legend=c("MM","Ricker"),bty="n")


##########
# 3a

######
# first set up the necessary data objects

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)  #Data set from grade 1 of myxo data
#head(Myx)

myx.data.for.bugs <- list(
  titer = Myx$titer,
  days = Myx$day,
  n.observations = nrow(Myx)
)

#myx.data.for.bugs

init.vals.for.bugs <- function(){
  list(
    shape=runif(1,20,30),
    a=runif(1,0.05,0.3),
    b=runif(1,1,1.3)
  )
}

init.vals.for.bugs2 <- function(){
  list(
    shape=runif(1,20,30),
    a=runif(1,8,9),
    b=runif(1,1,1.5)
  )
}

jags.fit_ricker <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs,parameters.to.save=params.to.store,n.iter=50000,model.file="BUGSmodel_ricker.txt",n.chains = 3,n.burnin = 5000,n.thin = 20 )

jags.fit_mm <- jags(data=myx.data.for.bugs,inits=init.vals.for.bugs2,parameters.to.save=params.to.store,n.iter=50000,model.file="BUGSmodel_mm.txt",n.chains = 3,n.burnin = 5000,n.thin = 20 )



Myx_PostPredCheck <- function(MCMC1=jags.fit_ricker,MCMC2=jags.fit_mm){
  lots <- 1000
  nMCMC <- min(jags.fit_ricker$BUGSoutput$n.keep,jags.fit_mm$BUGSoutput$n.keep)
  
     # hard-code the data
  MyxDat <- MyxoTiter_sum
  Myx <- subset(MyxDat,grade==1)  #Data set from grade 1 of myxo data  
  nobs <- nrow(Myx)
  
  SSEobs_rick <- numeric(lots)
  SSEobs_mm <- numeric(lots)
  SSEsim_rick <- numeric(lots)
  SSEsim_mm <- numeric(lots)
  
  params_rick <- matrix(NA, nrow=lots,ncol=3)
  params_mm <- matrix(NA, nrow=lots,ncol=3)
  
  colnames(params_rick) <- c("a","b","shape")
  colnames(params_mm) <- c("a","b","shape")
  
  alldays <- sort(unique(Myx$day))
  allsims_rick <- matrix(nrow=lots,ncol=length(alldays)) 
  allsims_mm <- matrix(nrow=lots,ncol=length(alldays)) 
  colnames(allsims_rick) <- alldays
  colnames(allsims_mm) <- alldays
  
  i=2
  for(i in 1:lots){
    random_index <- sample(1:nMCMC,1)
    params_rick[i,] <- c(a=MCMC1$BUGSoutput$sims.list$a[random_index],
                     b=MCMC1$BUGSoutput$sims.list$b[random_index],
                     shape=MCMC1$BUGSoutput$sims.list$shape[random_index])
    params_mm[i,] <- c(a=MCMC2$BUGSoutput$sims.list$a[random_index],
                     b=MCMC2$BUGSoutput$sims.list$b[random_index],
                     shape=MCMC2$BUGSoutput$sims.list$shape[random_index])
    
    ### loop through titer observations
    exptiter_rick <- Ricker(Myx$day,params_rick[i,"a"],params_rick[i,"b"]) 
    exptiter_mm <- MMFunc(Myx$day,params_mm[i,"a"],params_mm[i,"b"])
    simobs_rick <- rgamma(nobs,shape=params_rick[i,"shape"],rate=(params_rick[i,"shape"]/exptiter_rick))
    simobs_mm <- rgamma(nobs,shape=params_mm[i,"shape"],rate=(params_mm[i,"shape"]/exptiter_mm))
    SSEobs_rick[i] <- sum((Myx$titer-exptiter_rick)^2)
    SSEobs_mm[i] <- sum((Myx$titer-exptiter_mm)^2)
    SSEsim_rick[i] <- sum((simobs_rick-exptiter_rick)^2)
    SSEsim_mm[i] <- sum((simobs_mm-exptiter_mm)^2)
    
    allsims_rick[i,] <- tapply(simobs_rick,Myx$day,function(t) t[sample.int(length(t))][1] )
    allsims_mm[i,] <- tapply(simobs_mm,Myx$day,function(t) t[sample.int(length(t))][1] )
  }
  
  boxplot(allsims_rick,xlab="Days",ylab="Titer",xlim=c(0,10),main="Ricker",at=as.numeric(colnames(allsims_rick)))
  points(Myx$day,Myx$titer,col="green",pch=19,cex=2)
  
  boxplot(allsims_mm,xlab="Days",ylab="Titer",xlim=c(0,10),main="M-M",at=as.numeric(colnames(allsims_mm)))
  points(Myx$day,Myx$titer,col="green",pch=19,cex=2)
  
  plot(SSEsim_rick~SSEobs_rick,main="Posterior Pred Check, Ricker")
  abline(0,1,col="red")
  
  plot(SSEsim_mm~SSEobs_mm,main="Posterior Pred Check, M-M")
  abline(0,1,col="red")
  
  pval_rick <- length(which(SSEsim_rick>SSEobs_rick))/lots
  
  pval_mm <- length(which(SSEsim_mm>SSEobs_mm))/lots
  
  output <- list()
  output[["p-vals"]] <- c(pval_rick,pval_mm)
  
  output[["PPC_rick"]] <- as.data.frame(cbind(params_rick,SSEsim_rick,SSEobs_rick))
  output[["PPC_mm"]] <- as.data.frame(cbind(params_mm,SSEsim_mm,SSEobs_mm))
  
  colnames <- c("a","b","shape","SSEsim","SSEobs")
  names(output[["PPC_rick"]]) <- colnames
  names(output[["PPC_mm"]]) <- colnames
  
  return(output)
  
}


####
# Test the function

test <- Myx_PostPredCheck(MCMC1=jags.fit_ricker,MCMC2=jags.fit_mm)

test[[1]]   # p-vals

head(test[[2]])    # ppc for ricker

head(test[[3]])    # ppc for M-M



#######
# Question 4

### DIC

DIC_ricker <- jags.fit_ricker$BUGSoutput$DIC 
DIC_MM <- jags.fit_mm$BUGSoutput$DIC 

### WAIC

library(loo)    # for computing WAIC

nmcmc_ricker <- jags.fit_ricker$BUGSoutput$n.keep   # number of mcmc samples
avals <- jags.fit_ricker$BUGSoutput$sims.list$a
bvals <- jags.fit_ricker$BUGSoutput$sims.list$b
svals <- jags.fit_ricker$BUGSoutput$sims.list$shape
liks_ricker <- t(sapply(1:nmcmc_ricker,function(t) dgamma(Myx$titer,svals[t],rate=(svals[t]/Ricker(Myx$day,avals[t],bvals[t])),log = T)  ) )

nmcmc_mm <- jags.fit_mm$BUGSoutput$n.keep   # number of mcmc samples
avals <- jags.fit_mm$BUGSoutput$sims.list$a
bvals <- jags.fit_mm$BUGSoutput$sims.list$b
svals <- jags.fit_mm$BUGSoutput$sims.list$shape
liks_mm <- t(sapply(1:nmcmc_ricker,function(t) dgamma(Myx$titer,svals[t],rate=(svals[t]/MMFunc(Myx$day,avals[t],bvals[t])),log=T)  ) )

WAIC_ricker <- waic(liks_ricker)

WAIC_mm <- waic(liks_mm)



DIC_ricker
DIC_MM


WAIC_ricker$estimates["waic",]
WAIC_mm$estimates["waic",]

