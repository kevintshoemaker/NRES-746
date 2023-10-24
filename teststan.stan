// saved as teststan.stan

// The input data specification
data {
  int<lower=0> J;               // number of schools 
  array[J] real y;              // estimated treatment effects
  array[J] real<lower=0> sigma; // standard error of effect estimates 
}

// The parameters accepted by the model. 
parameters {
  real mu;                // population treatment effect
  real<lower=0> tau;      // standard deviation in treatment effects
  vector[J] eta;          // unscaled deviation from mu by school
}

transformed parameters {
  vector[J] theta = mu + tau * eta;        // school treatment effects
}

// The model to be estimated (likelihood function)
model {
  target += normal_lpdf(eta | 0, 1);       // prior log-density
  target += normal_lpdf(y | theta, sigma); // log-likelihood
}
