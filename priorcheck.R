
# prior predictive check  -----

 ## (KTS fleshed out after starting in class)


## load data ------

library(emdbook)
data(FirDBHFec)
fir <- na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
fir$TOTCONES <- round(fir$TOTCONES)
head(fir)
dim(fir)

library(ggplot2)
ggplot(fir,aes(DBH,TOTCONES)) + geom_point()

## functions -------

powerlaw <- function(x,a,b){     # deterministic relationship between mean reproduction and covariate
  a*x^b
}

powerlaw(10,10,0.5)

r=1000

priors <- function(n){     # generate N samples from prior distribution in parameter space
  a = abs(rnorm(n,5,2))
  b = abs(rnorm(n,1.5,0.25))  #   rgamma(n,10,8)
  k = abs(rnorm(n,4,0.5) )
  data.frame(a=a,b=b,k=k)
}
p = priors(10)
head(p)

sim_mu_from_prior <- function(dat,pars){
  thisa = as.numeric(pars[1])
  thisb =  as.numeric(pars[2])
  powerlaw(dat$DBH,thisa,thisb)
}

sim_dat_from_prior <- function(pars,dat){
  thisk = as.numeric(pars[3])
  mu = sim_mu_from_prior(dat,pars)
  rnbinom(length(mu),mu=mu,size=thisk)
}


## prior check of deterministic relationship  -------------

library(tidyverse)

fir$DBH2 = ceiling(fir$DBH)
df = fir |> 
  group_by(DBH2) |> 
  summarize(MU_CONES = mean(TOTCONES,na.rm=T)) |> 
  rename(DBH = DBH2)

p = priors(1000)

# sim_mu_from_prior(pars=p[1,],dat=df)

mu_pred = as.data.frame(apply(p,1,sim_mu_from_prior,dat=df)) |> 
  set_names(1:nrow(p)) |>
  bind_cols(DBH = df$DBH) |> 
  pivot_longer(-DBH,names_to = "rep", values_to = "MU_CONES")


ggplot(df,aes(DBH,MU_CONES)) + geom_path(lwd=3) +
  geom_path(data=mu_pred,aes(DBH,MU_CONES,group=rep),col="gray")
  


## prior check for stochastic relationship  ---------


d_pred = as.data.frame(apply(p,1,sim_dat_from_prior,dat=fir)) |> 
  set_names(1:nrow(p)) |>
  bind_cols(DBH = fir$DBH) |> 
  pivot_longer(-DBH,names_to = "rep", values_to = "TOTCONES") |> 
  mutate(DBH2 = ceiling(DBH))


ggplot(fir,aes(DBH,TOTCONES)) + #geom_point(cex=2) +
  geom_violin(data=d_pred,aes(DBH2,TOTCONES,group=DBH2),col="gray",alpha=0.5)













