

// noncentered_model.stan
data {
  int<lower=1> N;               // number of observations
  int<lower=1> K;               // number of groups
  array [N] int<lower=1, upper=K> gg; //index of group (for every group tell what observation in)
  vector[N] X;                  // predictor
  vector[N] y;                  // response
}

parameters {
  real alpha0;                      // population intercept
  real b1;                      // slope
  real<lower=0> tau;            // group-level SD
  real<lower=0> sigma;          // residual SD
  vector[K] alpha_raw;          // unscaled group effects
}

transformed parameters {
  vector[K] alpha;
  alpha = alpha0 + tau * alpha_raw;      // scale group effects by tau
}

model {
  // priors
  alpha_raw ~ std_normal();     // standard normal for raw group effects
  
  // likelihood
  y ~ normal(b1 * X + alpha[gg], sigma);
}
