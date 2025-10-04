
data {
  int<lower = 1> N;
  array [N] int<lower=0> obs_cones;
  vector<lower=0> [N] DBH;
}

parameters {
  real loga, b, logbeta;
}

transformed parameters {
  real a = exp(loga);
  real beta = exp(logbeta);
}

model  {
  vector[N] mean_cones = exp(loga + b .* DBH);   // power function: a*DBH^b
  vector[N] alpha = mean_cones .* beta;
  obs_cones ~ neg_binomial(alpha,beta);
}

generated quantities {   // need log_lik of each data point for model selection
  vector[N] log_lik; // N is the number of data points
  {
    real m2, a2, b2;
    for (n in 1:N) {
       m2 = exp(loga + b * DBH[n]);   // power function: a*DBH^b
       a2 = m2 * beta;
       log_lik[n] = neg_binomial_lpmf(obs_cones[n] | a2, beta);
    }
  }
}

