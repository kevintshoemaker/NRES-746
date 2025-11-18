
#  NRES 746, Lab 4: Bayesian inference  ------------                                 
#       University of Nevada, Reno                        
                               

data {
  int<lower=0> N;
  vector[N] titer, day;
}

parameters {
  real<lower=0> a,b,rate;
}

transformed parameters {
  vector[N] m = (a * day) .* exp( -b * day);   // mean expected titer  // note elementwise mult - '.*'
  vector[N] shape = rate * m;     // compute gamma shape parameter
}

model {
  a ~ exponential(.1);   // prior on a
  b ~ exponential(.1);   // prior on b
  rate ~ exponential(.1);   // prior on rate
  titer ~ gamma(shape, rate);   // likelihood
}

generated quantities {
  vector[N] log_lik; // N is the number of data points
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(titer[n] | shape[n], rate); // Replace with your model's likelihood
  }
}
    

## exercise 4.1 ---------

# test your stan code:

mod2 <- cmdstan_model("myx2.stan") # Compile stan model

stan_data <- list(   # bundle data for stan
    N = nrow(Myx),
    titer = Myx$titer,
    day=Myx$day
)

inits <- function(){
  list(
    a=runif(1,1,10),
    b=runif(1,.01,0.1),
    rate=runif(1,1,10)
  )
}
# inits()

fit2 <- mod2$sample(
  data = stan_data,
  chains = 4,
  init = inits,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0   # don't provide progress updates
)

fit2$summary()

samples_rick <- fit2$draws(format="draws_df")
bayesplot::mcmc_trace(samples_rick,"a")
bayesplot::mcmc_trace(samples_rick,"b")
bayesplot::mcmc_trace(samples_rick,"rate")



functions {
  vector mm(vector x, real a, real b){
    return( (a*x) ./ (b+x));
  } 
}

data {
  int<lower=0> N;
  vector[N] titer, day;
}

parameters {
  real<lower=0> a,b,rate;
}

transformed parameters {
  vector[N] m = mm(day,a,b);   // mean expected titer  
  vector[N] shape = rate * m;     // compute gamma shape parameter
}

model {
  a ~ exponential(.1);   // prior on a
  b ~ exponential(.1);   // prior on b
  rate ~ exponential(.1);   // prior on rate
  titer ~ gamma(shape, rate);   // likelihood
}

generated quantities {
  vector[N] log_lik; // N is the number of data points
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(titer[n] | shape[n], rate); // Replace with your model's likelihood
  }
}
    

## exercise 4.1 ---------

# test your stan code:

mod3 <- cmdstan_model("myx3.stan") # Compile stan model

stan_data <- list(   # bundle data for stan
    N = nrow(Myx),
    titer = Myx$titer,
    day=Myx$day
)

inits <- function(){
  list(
    a=runif(1,1,10),
    b=runif(1,.01,0.1),
    rate=runif(1,1,10)
  )
}
# inits()

fit3 <- mod3$sample(
  data = stan_data,
  chains = 4,
  init = inits,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 0   # don't provide progress updates
)

fit3$summary()

samples_mm <- fit3$draws(format="draws_df")
bayesplot::mcmc_trace(samples_mm,"a")
bayesplot::mcmc_trace(samples_mm,"b")
bayesplot::mcmc_trace(samples_mm,"rate")




######
# 2b

ricker = function(x,a,b) (a * x) * exp( -b * x)

plot(Myx$titer~Myx$day,xlim=c(0,10),ylim=c(0,10))
curve(mm(x,a=mean(samples_mm$a),b=mean(samples_mm$b)),add=T,col="red",lwd=2)
curve(ricker(x,a=mean(samples_rick$a),b=mean(samples_rick$b)),add=T,col="green",lty=2,lwd=2)
legend("topleft",lwd=c(2,2),lty=c(1,2),col=c("red","green"),legend=c("MM","Ricker"),bty="n")


## question 4_3  -------------

simfun <- function(par,df,predfun,r){
  df$titer_e = predfun(x=df$day,a=par$a,b=par$b)
  shape = df$titer_e*par$rate
  df$titer_p = rgamma(nrow(df),shape,par$rate)
  df$res1 = df$titer - df$titer_e      # obs - exp
  df$res2 = df$titer_p - df$titer_e     # pred - exp
  df$rep = r
  return(df)
}

# MCMC=samples_rick; predfun=ricker; dat=Myx   # for debug
Myx_PostPredCheck <- function(MCMC,predfun,dat){
  lots <- 1000; nMCMC <- length(MCMC$a); nobs = nrow(dat); nobs <- nrow(Myx)
  ret = list()   # initialize return list
  parnames = c("a","b","rate")   # hard coding the parameter names (not best coding practice!)
  ndx = sample(1:nMCMC,lots,replace = T)
  params <- as.data.frame(MCMC)[ndx,parnames]     
  reps = lapply(1:lots, function(t) simfun(par=params[t,c("a","b","rate")],df=dat,predfun,t)     ) 
  reps = do.call("rbind",reps)   # put everything into a bit data frame
  ppc1 <- reps |> 
    group_by(rep) |> 
    summarise(
      SSE_obs = sum(res1^2),
      SSE_sim = sum(res2^2)
    )
  ppc1 <- cbind(ppc1,params)
  g1 = ggplot(ppc1,aes(SSE_obs,SSE_sim)) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, lty=2, lwd=2, col = "darkgreen") 
  ppc2 <- reps |> 
    group_by(day) |> 
    summarise(
      mean = mean(titer_p),
      upper = quantile(titer_p,0.975),
      lower = quantile(titer_p,0.025)
    )
  g2 = ggplot(ppc2,aes(day,mean)) + geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2) + 
    geom_point(data = dat, aes(day,titer),shape = 1,col="darkgreen",size=5) + labs(y="Titer") +
    geom_point(size=3,shape=20) + geom_path(lty=2)
  
  print(cowplot::plot_grid(g1,g2))
  
  ret$bayesian_p = sum(ppc1$SSE_sim>ppc1$SSE_obs)/lots
  ret$df = ppc1
  return(ret)
}


####
# Test the function

test <- Myx_PostPredCheck(MCMC=samples_mm,mm,Myx)
test$bayesian_p
head(test$df)


## Question 4: WAIC and LOO-IC with PSIS  ------------

library(loo)

test_rick <- fit2$loo()    # extract pointwise log likelihoods and compute 'loo' metrics
test_mm <- fit3$loo()

test_rick;test_mm



## Question 4: WAIC and LOO-IC with PSIS  ------------

library(loo)

test_rick <- fit2$loo()    # extract pointwise log likelihoods and compute 'loo' metrics
test_mm <- fit3$loo()

test_rick;test_mm

