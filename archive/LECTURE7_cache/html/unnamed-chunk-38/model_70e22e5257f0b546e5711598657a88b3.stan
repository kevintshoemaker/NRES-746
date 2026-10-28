
data {
  int<lower=0> N;
  vector[N] titer;
}

parameters {
  real<lower=0> shape;
  real<lower=0> rate;
}

model {
  shape ~ gamma(0.001,0.001);   // prior on shape
  rate ~ gamma(0.001,0.001);   // prior on rate
  titer ~ gamma(shape, rate);   // likelihood
}

