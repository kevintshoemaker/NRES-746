#######
#Data brought in to test steps:
#######
params=c(a=2,b=0.2,shape=30)

params_selected = c("a","b")

param1_lims= c(0.1,9)

param2_lims=c(0.01,0.5)
library(emdbook)
data(MyxoTiter_sum)      # load the data
head(MyxoTiter_sum)   

myxdat <- subset(MyxoTiter_sum, grade==1)    # select just the least virulent strain

plot(myxdat$titer~myxdat$day,xlim=c(0,10))    # visualize the relationship

myxdat <- subset(MyxoTiter_sum, grade==1)    # select just the least virulent strain

plot(myxdat$titer~myxdat$day,xlim=c(0,10))    # visualize the relationship
params = c(a=4,b=-.2,shape=40)
data=myxdat[,-1]


################################################
#STEPS:

#MLE for parameters
NLL_myxRicker <- function(params,data) {
  meantiter = params[1] * myxdat$day * exp(-params[2] * myxdat$day)
  -sum(dgamma(myxdat$titer,shape=params[3],scale=meantiter/params[3],log=TRUE))
}

opt2 <- suppressWarnings(optim(fn=NLL_myxRicker, params))  #use default simplex algorithm
opt2par <- opt2$par

opt2par #MLE for a, b and shape


param1<- params_selected[1]
param2<- params_selected[2]
select_params<- c(param1,param2)
select_params<- params[select_params]#isolated selected parameters-- a and b
subvars <- names(params) %in% c(param1,param2)
const_param <- opt2par[!subvars] #MLE of constant parameter-- shape
new_params <- c(select_params[1],select_params[2],const_param[1])
new_params <- as.list(new_params)#new list containing selected a, selected b and MLE shape

#limit sequence
params_selected1 <- seq(param1_lims[1],param1_lims[2],length=50)
params_selected2 <- seq(param2_lims[1],param2_lims[2],length=50)

loglikelihood_surface2 <- matrix(0,nrow=50,ncol=50)   # set up storage matrix!

new_params
#Fill in matrix
i=1;j=1
for(i in 1:length(params_selected1)){
  new_params[[1]] <- params_selected1[i]
  for(j in 1:length(params_selected2)){
    new_params[[2]] <- params_selected2[j]
    loglikelihood_surface2[i,j] <- NLL_myxRicker(unlist(new_params),myxdat)
  }
}
loglikelihood_surface2
