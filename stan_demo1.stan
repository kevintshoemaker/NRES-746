
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real b0,b1;
  real<lower=0> sigma;
}

model {
  vector [N] mu = b0 + b1*x;
  y ~ normal(mu, sigma);
}

