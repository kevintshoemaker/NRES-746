

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N; //number of observations
  int<lower=1> K; //number of groups 
  array [N] int<lower=1, upper=K> gg; //index of group (for every group tell what observation in)
  vector[N] X; //covariate
  vector[N] y;
}

parameters{
  real alpha0;
  real b1;
  real<lower=0> tau; //among group variaibility
  real<lower=0> sigma; 
  vector[K] alpha;
}

//model block
model{
  
  alpha~ normal(alpha0, tau);
  
  vector[N] LP = b1*X;
  y ~ normal (LP + alpha[gg], sigma);

}


