
#  NRES 746, Lecture 9                      
#  University of Nevada, Reno   
#  Model performance evaluation  --------------------------


library(emdbook)

MyxDat <- MyxoTiter_sum
Myx <- subset(MyxDat,grade==1)  #Data set from grade 1 of myxo data
head(Myx)


# Fit the model with ML -----------------------------

Ricker <- function(a,b,predvar) a*predvar*exp(-b*predvar)
  
NLL_myx <- function(params){
  expected <- Ricker(params[1],params[2],Myx$day)
  -sum(dgamma(Myx$titer,shape=expected*params[3],rate=params[3],log = T))
}

params <- c(a=1,b=0.2,rate=1)
NLL_myx(params)

opt <- optim(params, NLL_myx)

MLE = opt$par
maxlik = opt$value

MLE
  

# Plug-in prediction interval -------------------------

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,15))
expected <- Ricker(MLE['a'],MLE['b'],1:10)
points(1:10,expected,type="l",col="green")

upper <- qgamma(0.975,shape=expected*MLE['rate'],rate=MLE['rate'])
lower <- qgamma(0.025,shape=expected*MLE['rate'],rate=MLE['rate'])

points(1:10,upper,type="l",col="red",lty=2)
points(1:10,lower,type="l",col="red",lty=2)


# Parametric bootstrap!  -------------------------------------

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,15),type="n")
expected <- Ricker(MLE['a'],MLE['b'],1:10)
points(1:10,expected,type="l",col="darkgreen")

uniquedays <- sort(unique(Myx$day))
expected <- Ricker(MLE['a'],MLE['b'],uniquedays)

simdata = t(replicate(1000,rgamma(length(uniquedays),shape=expected * MLE['rate'],rate=MLE['rate'])))

upper <- apply(simdata,2,function(t) quantile(t,0.975))
lower <- apply(simdata,2,function(t) quantile(t,0.025))

points(uniquedays,upper,type="l",col="red",lty=2)
points(uniquedays,lower,type="l",col="red",lty=2)

boxplot(x=as.list(as.data.frame(simdata)),at=uniquedays,add=T,boxwex=0.25,xaxt="n",range=0,col="red")
points(Myx$day,Myx$titer,cex=1.5,pch=20)


# Compare observed error statistic with expected range of error statistic as part of parametric bootstrap analysis

expected <- Ricker(MLE['a'],MLE['b'],Myx$day)
simdata = t(replicate(1000,rgamma(length(Myx$day),shape=expected * MLE['rate'],rate=MLE['rate'])))
 
rmse_observed <- sqrt(mean((Myx$titer-expected)^2))
rmse_simulated <- apply(simdata,1,function(t) mean((t-expected)^2))

hist(rmse_simulated,freq=F)
abline(v=rmse_observed,col="darkgreen",lwd=3)


# Bayesian goodness-of-fit  -------------------------

library(cmdstanr)
library(posterior)
library(bayesplot)
options(mc.cores=4)


# test your stan code:

mod1 <- cmdstan_model("myx2.stan") # Compile stan model

stan_data <- list(   # bundle data for stan
    N = nrow(Myx),
    titer = Myx$titer,
    day=Myx$day
)

fit <- mod1$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 0   # don't provide progress updates
)

fit$summary()

samples <- fit$draws(format="draws_df")
bayesplot::mcmc_trace(samples,"a"); bayesplot::mcmc_trace(samples,"b")
bayesplot::mcmc_trace(samples,"rate")



# set range of days
days=1:10

nMCMC = length(samples$a)

# make predictions from posterior
n_samp = 500
do_pred = function(){
  grab = samples[sample(1:nMCMC,1),]
  with(grab,rgamma(length(days), shape=Ricker(a,b,days)*rate,rate=rate) )
}
post_pred = as.data.frame(t(replicate(n_samp,do_pred())))

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,15),type="n")
expected <- Ricker(mean(samples$a),mean(samples$b),1:10)
points(1:10,expected,type="l",col="darkred",lwd=2)

boxplot(x=as.list(post_pred),at=1:10,add=T,boxwex=0.25,xaxt="n",range=0,border="darkred")
points(Myx$day,Myx$titer,cex=1.5,pch=20)



n_samp = 500
do_pred2 = function(){
  grab = samples[sample(1:nMCMC,1),]
  thismean = with(grab, Ricker(a,b,Myx$day))
  thissim = with(grab, rgamma(nrow(Myx),thismean*rate,rate=rate ))
  RMSE_obs = sqrt(mean((Myx$titer-thismean)^2))
  RMSE_sim = sqrt(mean((thissim-thismean)^2))
  c(RMSE_obs = RMSE_obs,RMSE_sim = RMSE_sim)
}
post_pred2 = as.data.frame(t(replicate(n_samp,do_pred2())))

plot(post_pred2$RMSE_sim~post_pred2$RMSE_obs, main="posterior predictive check")
abline(0,1,col="red",lwd=2)
p.value=mean(post_pred2$RMSE_sim>post_pred2$RMSE_obs)
p.value


# Summary statistics of a models "usefulness" (e.g., R-squared)  -------------------

SS_res <- sum((Myx$titer-Ricker(MLE["a"],MLE["b"],Myx$day))^2)
SS_tot <- sum((Myx$titer-mean(Myx$titer))^2)
Rsquared <- 1-SS_res/SS_tot

cat("R-squared = ", Rsquared, "\n")

# Fit the null likelihood model!

NLL_null <- function(params){
  -sum(dgamma(Myx$titer,shape=params[1]*params[2],rate=params[2],log = T))
}

init.params <- c(mean=7,rate=1)

opt_null <- optim(par=init.params, fn=NLL_null)
MLE_null=opt_null$par
maxlik_null = opt_null$value

McFadden <- 1-(maxlik/maxlik_null)
cat("McFadden's R-squared = ", McFadden) 


RMSE = sqrt(mean((Myx$titer-Ricker(MLE["a"],MLE["b"],Myx$day))^2))
cat("RMSE = ", RMSE, "\n")


# Collect new data that were not used in model fitting

newdata <- data.frame(
  grade = 1,
  day = c(2,3,4,5,6,7,8),
  titer = c(4.4,7.2,6.8,5.9,9.1,8.3,8.8)
)
newdata


# Validation #1  -------------------------

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,15),type="n",xlab="days",ylab="titer")
expected <- Ricker(MLE['a'],MLE['b'],1:10)
points(1:10,expected,type="l",col="green")

expected <- Ricker(MLE['a'],MLE['b'],1:10)
simdata <- array(0,dim=c(1000,10))
for(i in 1:1000){
  simdata[i,] <- rgamma(10,shape=expected*MLE['rate'],rate=MLE['rate'])
}

upper <- apply(simdata,2,function(t) quantile(t,0.975))
lower <- apply(simdata,2,function(t) quantile(t,0.025))

points(1:10,upper,type="l",col="green",lty=2)
points(1:10,lower,type="l",col="green",lty=2)

boxplot(x=as.list(as.data.frame(simdata)),at=1:10,add=T,boxwex=0.25,xaxt="n",range=0,border="green")
points(newdata$day,newdata$titer,cex=1.5,pch=20,col="red")
points(Myx$day,Myx$titer,cex=1.5,pch=20,col="black")
legend("topleft",pch=c(20,20),col=c("black","red"),legend=c("original data","validation data"))


SS_res <- sum((newdata$titer-Ricker(MLE["a"],MLE["b"],newdata$day))^2)
SS_tot <- sum((newdata$titer-mean(newdata$titer))^2)
Rsquared_validation <- 1-SS_res/SS_tot

cat("R-squared = ", Rsquared_validation, "\n")

expected <- Ricker(MLE["a"],MLE["b"],newdata$day)
McFadden_validation <- 1-(sum(dgamma(newdata$titer,shape=expected*MLE["rate"],rate=MLE["rate"], log = T))/sum(dgamma(newdata$titer,shape=MLE_null["mean"]*MLE_null['rate'],rate=MLE_null["rate"],log=T)))
cat("pseudo R-squared = ", McFadden_validation, "\n")

RMSE = sqrt(mean((newdata$titer-Ricker(MLE["a"],MLE["b"],newdata$day))^2))
cat("RMSE = ", RMSE, "\n")


# Validation #2 ---------------------------

newdata <- data.frame(        # imagine these are new observations...
  grade = 1,
  day = c(10,11,12,13,14,15,16),
  titer = c(6.8,8.0,4.5,3.1,2.7,1.2,0.04)
)
newdata


plot(Myx$titer~Myx$day,xlim=c(0,20),ylim=c(0,15),type="n",xlab="days",ylab="titer")
expected <- Ricker(MLE['a'],MLE['b'],1:20)
points(1:20,expected,type="l",col="green")

expected <- Ricker(MLE['a'],MLE['b'],1:20)
simdata <- array(0,dim=c(1000,20))
for(i in 1:1000){
  simdata[i,] <- rgamma(20,shape=expected*MLE['rate'],scale=MLE['rate'])
}

upper <- apply(simdata,2,function(t) quantile(t,0.975))
lower <- apply(simdata,2,function(t) quantile(t,0.025))

points(1:20,upper,type="l",col="green",lty=2)
points(1:20,lower,type="l",col="green",lty=2)

boxplot(x=as.list(as.data.frame(simdata)),at=1:20,add=T,boxwex=0.25,xaxt="n",range=0,border="green")
points(newdata$day,newdata$titer,cex=1.5,pch=20,col="red")
points(Myx$day,Myx$titer,cex=1.5,pch=20,col="black")
legend("topleft",pch=c(20,20),col=c("black","red"),legend=c("original data","new data"))


SS_res <- sum((newdata$titer-Ricker(MLE["a"],MLE["b"],newdata$day))^2)
SS_tot <- sum((newdata$titer-mean(newdata$titer))^2)
Rsquared_validation <- 1-SS_res/SS_tot

cat("R-squared = ", Rsquared_validation, "\n")

expected <- Ricker(MLE["a"],MLE["b"],newdata$day)
McFadden_validation <- 1-(sum(dgamma(newdata$titer,shape=expected*MLE["rate"],rate=MLE["rate"], log = T))/sum(dgamma(newdata$titer,shape=MLE_null["mean"]*MLE_null["rate"],rate=MLE_null["rate"],log=T)))
cat("pseudo R-squared = ", McFadden_validation, "\n")

RMSE = sqrt(mean((newdata$titer-Ricker(MLE["a"],MLE["b"],newdata$day))^2))
cat("RMSE = ", RMSE, "\n")


# CROSS-VALIDATION ------------------


# PARTITION THE DATA

n.folds <- nrow(Myx)   # jackknife

Myx$fold <- sample(c(1:n.folds),size=nrow(Myx),replace=FALSE)

init.params <- c(a=1,b=0.2,rate=1)

   # new likelihood function that takes data as an argument
NLL_myx = function(params,data){
  expected <- Ricker(params[1],params[2],data$day)
  -sum(dgamma(data$titer,shape=expected*params[3],rate=params[3],log = T))
}

Myx$pred_CV <- 0
for(i in 1:n.folds){
  Myx2 <- subset(Myx,fold!=i)   # observations to use for fitting 
  newfit <- optim(par=init.params, fn=NLL_myx, data=Myx2)   # fit the model, leaving out this partition
  ndx <- Myx$fold == i
  Myx$pred_CV[ndx] <- Ricker(newfit$par['a'],newfit$par['b'],Myx$day[ndx])
}

Myx$pred_full <- Ricker(MLE['a'],MLE['b'],Myx$day)

Myx


RMSE_full <- sqrt(mean((Myx$titer-Myx$pred_full)^2))
RMSE_CV <- sqrt(mean((Myx$titer-Myx$pred_CV)^2))
RMSE_full
RMSE_CV


VarExplained_full = 1 - mean((Myx$titer-Myx$pred_full)^2)/mean((Myx$titer-mean(Myx$titer))^2)
VarExplained_CV = 1 - mean((Myx$titer-Myx$pred_CV)^2)/mean((Myx$titer-mean(Myx$titer))^2)
VarExplained_full
VarExplained_CV


# Cross-validation: titanic disaster example!  --------------------

# titanic <- read.csv("titanic.csv",header=T)
# head(titanic)

library(titanic)  
titanic <- titanic_train


titanic2 <- na.omit(titanic)
titanic2$Pclass = factor(titanic2$Pclass,ordered=T)
model1 <- glm(Survived ~ Sex + scale(Age) + scale(SibSp) + scale(Parch) + scale(Fare), data=titanic2, family="binomial")    #logistic regression
summary(model1)


params <- c(
  int=1,
  male = -1,
  age = 0,
  sibsp = 0,
  parch = 0,
  fare = 0
)

LikFunc <- function(params){
  linear <- params['int'] + 
    params['male']*as.numeric(titanic2$Sex=="male") +
    params['age']*scale(titanic2$Age) +
    params['sibsp']*scale(titanic2$SibSp) +
    params['parch']*scale(titanic2$Parch) +
    params['fare']*scale(titanic2$Fare)
  meanprob <-  1/(1+exp(-(linear)))
  -sum(dbinom(titanic2$Survived,size=1,prob = meanprob,log=T))
}

LikFunc(params)

MLE <- optim(fn=LikFunc,par = params)

MLE$par


SibSp_range <- range(titanic$SibSp)
Parch_range <- range(titanic$Parch)
Fare_range <- range(titanic$Fare)
Age_range <- range(titanic$Age,na.rm = T)


### 

plot(titanic$Survived~titanic$Fare,pch=16,xlab="FARE ($)",ylab="Survived!")

predict_df <- data.frame(
  Sex = "male",
  Age = mean(titanic$Age,na.rm=T),
  SibSp = mean(titanic$SibSp),
  Parch = mean(titanic$Parch),
  Fare = seq(Fare_range[1],Fare_range[2])
)

probSurv <- predict(model1,predict_df,type="response")

lines(seq(Fare_range[1],Fare_range[2]),probSurv)



### 

predict_df <- data.frame(
  Sex = c("male","female"),
  Age = mean(titanic$Age,na.rm=T),
  SibSp = mean(titanic$SibSp),
  Parch = mean(titanic$Parch),
  Fare = mean(titanic$Fare,na.rm=T)
)

tapply(titanic$Survived,titanic$Sex,mean)[2:1]

probSurv <- predict(model1,predict_df,type="response")
names(probSurv) <- c("male","female")

probSurv



plot(titanic$Survived~titanic$Age,pch=16,xlab="AGE",ylab="Survived!")

predict_df <- data.frame(
  Sex = "male",
  Age = seq(Age_range[1],Age_range[2]),
  SibSp = mean(titanic$SibSp),
  Parch = mean(titanic$Parch),
  Fare = mean(titanic$Fare,na.rm=T)
)

probSurv <- predict(model1,predict_df,type="response")

lines(seq(Age_range[1],Age_range[2]),probSurv)


### 

plot(titanic$Survived~jitter(titanic$SibSp),pch=16,xlab="# of Siblings/spouses",ylab="Survived!")

predict_df <- data.frame(
  Sex = "male",
  Age = mean(titanic$Age,na.rm=T),
  SibSp = seq(SibSp_range[1],SibSp_range[2],0.01),
  Parch = mean(titanic$Parch),
  Fare = mean(titanic$Fare,na.rm=T)
)

probSurv <- predict(model1,predict_df,type="response")

lines(seq(SibSp_range[1],SibSp_range[2],0.01),probSurv)


library(ROCR)
library(rms)


model1 <- glm(Survived ~ Sex + SibSp + Parch + Fare, data=titanic, family="binomial")


# CROSS VALIDATION CODE FOR BINARY RESPONSE  ----------------------------

n.folds = 10       # set the number of "folds"
foldVector = rep(c(1:n.folds),times=floor(length(titanic$Survived)/9))[1:length(titanic$Survived)]


CV_df <- data.frame(
  CVprediction = numeric(nrow(titanic)),      # make a data frame for storage
  realprediction = 0,
  realdata = 0
)

for(i in 1:n.folds){
  fit_ndx <- which(foldVector!=i)
  validate_ndx <- which(foldVector==i)
  model <- glm(formula = Survived ~ Sex + SibSp + Parch + Fare, family = "binomial", data = titanic[fit_ndx,]) 
  CV_df$CVprediction[validate_ndx] <-  plogis(predict(model,newdata=titanic[validate_ndx,])) 
  CV_df$realprediction[validate_ndx]  <-  plogis(predict(model1,newdata=titanic[validate_ndx,]))
  CV_df$realdata[validate_ndx] <- titanic$Survived[validate_ndx]
}

CV_RMSE = sqrt(mean((CV_df$realdata - CV_df$CVprediction)^2))       # root mean squared error for holdout samples in 10-fold cross-validation
real_RMSE = sqrt(mean((CV_df$realdata - CV_df$realprediction)^2))  # root mean squared error for residuals from final model

# print RMSE statistics

cat("The RMSE for the model under cross-validation is: ", CV_RMSE, "\n")

cat("The RMSE for the model using all data for training is: ", real_RMSE, "\n")
   


head(CV_df)


library(ROCR)

par(mfrow=c(2,1))
pred <- prediction(CV_df$CVprediction,CV_df$realdata)     # for holdout samples in cross-validation
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
plot(perf, main="Cross-validation")
text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))

pred <- prediction(CV_df$realprediction,CV_df$realdata)     # for final model
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
plot(perf, main="All data")
text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))


CV_df$CVprediction[which(CV_df$CVprediction==1)] <- 0.9999       # ensure that all predictions are not exactly 0 or 1
CV_df$CVprediction[which(CV_df$CVprediction==0)] <- 0.0001
CV_df$realprediction[which(CV_df$realprediction==1)] <- 0.9999
CV_df$realprediction[which(CV_df$realprediction==0)] <- 0.0001

fit_deviance_CV <- mean(-2*(dbinom(CV_df$realdata,1,CV_df$CVprediction,log=T)-dbinom(CV_df$realdata,1,CV_df$realdata,log=T)))
fit_deviance_real <- mean(-2*(dbinom(CV_df$realdata,1,CV_df$realprediction,log=T)-dbinom(CV_df$realdata,1,CV_df$realdata,log=T)))
null_deviance <- mean(-2*(dbinom(CV_df$realdata,1,mean(CV_df$realdata),log=T)-dbinom(CV_df$realdata,1,CV_df$realdata,log=T)))
deviance_explained_CV <- (null_deviance-fit_deviance_CV)/null_deviance   # based on holdout samples
deviance_explained_real <- (null_deviance-fit_deviance_real)/null_deviance   # based on full model...

# print RMSE statistics

cat("The McFadden R2 for the model under cross-validation is: ", deviance_explained_CV, "\n")

cat("The McFadden R2 for the model using all data for training is: ", deviance_explained_real, "\n")
 



