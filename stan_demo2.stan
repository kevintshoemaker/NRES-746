

// multilevel with random intercept


data {
  int<lower=1> N;   
  int<lower=1> G;   // number of groups
  vector[N] x;
  vector[N] y;
  array [N] int<lower=1> gg; // index of group by individual
}

parameters {
  real b0,b1;
  real<lower=0> sigma, tau;    // std dev terms:  note that tau is a 'hyperparameter'
  vector [G] gamma;    // random intercepts
}

model {
  b0 ~ normal(0,1);    
  b1 ~ normal(0,1);
  sigma ~ exponential(1);
  tau ~ exponential(1);
  
  // intercept for each group
  
  gamma ~ normal(0,tau);   // group effect
  
  vector [N] mu = b0 + b1*x + gamma[gg];
  y ~ normal(mu, sigma);
}

