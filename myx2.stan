data {
  int<lower=0> N;
  vector[N] titer, day;
}

parameters {
  real<lower=0> a,b,shape;
}

transformed parameters {
  vector[N] m = (a * day) .* exp( -b * day);   // mean expected titer  // note elementwise mult - '.*'
  vector[N] rate = shape ./ m;     // gamma rate = shape / mean (shape held constant)
}

model {
  a ~ exponential(.1);   // prior on a
  b ~ exponential(.1);   // prior on b
  shape ~ exponential(.01);   // weak prior on shape (prior mean 100)
  titer ~ gamma(shape, rate);   // likelihood
}

generated quantities {
  vector[N] log_lik; // N is the number of data points
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(titer[n] | shape, rate[n]); // Replace with your model's likelihood
  }
}
