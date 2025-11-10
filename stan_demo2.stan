
// data block: these variables must be fed to stan as a list from R
data {
  int<lower=0> N;    // sample size
  int<lower=0> G;    // number of groups
  vector[N] x;      // covariate vector
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
  corr_matrix[2] Omega;   // correlation matrix for random intercept and slopes
  real alpha0;              // global mean intercept term
  real beta0;             // global mean regression coefficients
  vector<lower=0.0>[2] tau;      // among-group variation for intercept and slope terms
  array[G] real alpha;      // group level intercept terms
  array[G] real beta;  // group level slope terms
}

transformed parameters {
  vector[2] mu0 = [alpha0, beta0]';
  array[G] vector[2] group_betas;
  for(g in 1:G) group_betas[g] = [alpha[g],beta[g]]';
}

model {
         // hyperpriors
  tau ~ exponential(1);
  alpha0 ~ normal(0,1);
  beta0 ~ normal(0,1);
  Omega ~ lkj_corr(4);   // set LKJ prior on correlation matrix
  
         // priors
  sigma ~ exponential(1);
  
  for(g in 1:G){
    group_betas[g] ~ multi_normal(mu0, quad_form_diag(Omega, tau));
  }
  
  for (n in 1:N) {      // data likelihood
    y[n] ~ normal(group_betas[gg[n]][1] + group_betas[gg[n]][2]*x[n], sigma);
  }

}

generated quantities {
      // empty for now
}

