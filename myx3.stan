functions {
  vector mm(vector x, real a, real b){
    return( (a*x) ./ (b+x));
  }
}

data {
  int<lower=0> N;
  vector[N] titer, day;
}

parameters {
  real<lower=0> a,b,rate;
}

model {
  a ~ exponential(0.1);   // prior on a
  b ~ exponential(0.1);   // prior on b
  rate ~ exponential(0.1);   // prior on rate
  vector[N] m = mm(day,a,b);   // mean expected titer  // note elementwise mult - '.*'
  vector[N] shape = rate * m;     // compute gamma shape parameter
  titer ~ gamma(shape, rate);   // likelihood
}
