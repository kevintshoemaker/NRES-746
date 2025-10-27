
// data block: these variables must be fed to stan as a list from R
data {
  int<lower=0> N;    // sample size
  int<lower=0> G;    // number of groups
  int<lower=0> K;    // number of predictors
  matrix[N,K] X;      // covariate matrix
  vector[N] y;       // response vector
  array [N] int<lower=1,upper=G> gg;   // index of group for each individual
}

// transformed data block: any data post processing or setting of constants that stan can use through the simulation: only run once per stan simulation
transformed data {
       // empty for now
}

// parameters block: define the dimensionality of the skate park for HMC to explore
parameters {
  real<lower=0.0> sigma;    // residual standard deviation
  corr_matrix[K+1] Omega;   // correlation matrix for random intercept and slopes
  real alpha0;              // global mean intercept term
  vector[K] beta0;           // global mean regression coefficients
  real<lower=0.0> tau_alpha;      // among-group variation for random intercept
  vector<lower=0.0>[K] tau_betas;      // among-group variation for slope terms
  array[G] real alpha;      // group level intercept terms
  array[G] vector[K] betas;  // group level slope terms
}

transformed parameters {
  vector[K+1] beta = append_row(alpha0, beta0);
  vector[K+1] tau = append_row(tau_alpha, tau_betas);  // standard deviations
  array[G] vector[K+1] mu_beta;
  for(g in 1:G) mu_beta[g] = append_row(alpha[g],betas[g]);
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

