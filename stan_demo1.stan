
// data block: these variables must be fed to stan as a list from R
data {
  int<lower=0> N;    // sample size
  int<lower=0> G;    // number of groups
  vector[N] x1;      // covariate vector
  vector[N] y;       // response vector
  array [N] int<lower=1,upper=G> gg;   // index of group for each individual
}

// transformed data block: any data post processing or setting of constants that stan can use through the simulation: only run once per stan simulation
transformed data {
       // empty for now
}

// parameters block: define the dimensionality of the skate park for HMC to explore
parameters {
  real b1;                  // regression coefficient 
  real<lower=0.0> sigma;    // residual standard deviation
  real alpha0;              // global mean intercept term
  real<lower=0.0> tau;      // among-group variation for random intercept
  vector[G] alpha;          // group level intercept terms
}

transformed parameters {
       // empty for now
}

model {
         // hyperpriors
  tau ~ exponential(1);
  alpha0 ~ normal(0,1);
  
         // priors
  b1 ~ normal(0,1);     
  sigma ~ exponential(1);
  alpha ~ normal(alpha0,tau);            // 'prior' on intercept with partial pooling
 
  vector [N] mu = alpha[gg] + b1*x1;    // linear predictor
  y ~ normal(mu, sigma);                // data likelihood
}

generated quantities {
      // empty for now
}

