data {
  int<lower=0> N;
  vector[N] titer, day;
}

parameters {
  real<lower=0> a,b,rate;
}

model {
  a ~ exponential(.1);   // prior on a
  b ~ exponential(.1);   // prior on b
  rate ~ exponential(.1);   // prior on rate
  vector[N] m = (a * day) .* exp( -b * day);   // mean expected titer
  vector[N] shape = rate * m;     // compute gamma shape parameter
  titer ~ gamma(shape, rate);   // likelihood
}
  
