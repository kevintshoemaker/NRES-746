
data {
  int<lower = 1> N;
  array [N] int<lower=0> obs_cones;
  vector<lower=0> [N] DBH;
  array [N] int<lower=1,upper=2> wave_ndx; 
}

parameters {
  vector[2] loga, b, logbeta;
}

transformed parameters {
  vector[2] a = exp(loga);
  vector[2] beta = exp(logbeta);
}

model  {
  vector[N] mean_cones = a[wave_ndx] .* DBH .^ b[wave_ndx];   // power function: a*DBH^b
  vector[N] alpha = mean_cones .* beta[wave_ndx];
  obs_cones ~ neg_binomial(alpha,beta[wave_ndx]);
}

generated quantities {   // need log_lik of each data point for model selection
  vector[N] log_lik; // N is the number of data points
  {
    real m2, a2;
    for (n in 1:N) {
       m2 = a[wave_ndx[n]] * DBH[n]^b[wave_ndx[n]];   // power function: a*DBH^b
       a2 = m2 * beta[wave_ndx[n]];
       log_lik[n] = neg_binomial_lpmf(obs_cones[n] | a2, beta[wave_ndx[n]]);
    }
  }
}


