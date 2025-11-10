data {
  int<lower=0> N;                    // number of training observations
  int<lower=0> N_test;               // number of observations in test dataset
  int<lower=0> K;                    // number of predictors
  int<lower=0> C;                    // number of classes
  array[N] int<lower=1,upper=C> cc;   // class index variable
  array[N_test] int<lower=1,upper=C> cc_test;   // class index variable
  matrix[N, K] X;                    // predictor matrix
  matrix[N_test, K] X_test;          // predictor matrix (test dataset)
  array[N] int<lower=0,upper=1> y;   // binary outcome
  array[N_test] int<lower=0,upper=1> y_test;   // binary outcome
  real<lower=0> lambda;              // L2 regularization parameter
}

parameters {
  vector[C] alpha;                   // intercept (one per class)
  vector[K] beta;                    // coefficients
}

model {
  // L2 regularization prior (ridge regression)
  alpha ~ normal(0, 10);             // weakly informative prior on intercept
  beta ~ normal(0, 1.0 / sqrt(lambda));  // L2 penalty on coefficients
  
  // Likelihood
  y ~ bernoulli_logit(alpha[cc] + X * beta);
}

generated quantities {
  vector[N] log_lik;                 // log-likelihood for each observation
  array[N_test] int y_pred;               // predicted classes
  int tp = 0;
  int tn = 0;
  real accuracy;
  
  {
    real eta;
    for (n in 1:N) {
      eta = alpha[cc[n]] + X[n] * beta;
      log_lik[n] = bernoulli_logit_lpmf(y[n] | eta);
    }
 
    for (n in 1:N_test) {
      eta = alpha[cc_test[n]] + X_test[n] * beta;
      y_pred[n] = bernoulli_rng(inv_logit(eta));
      if(y_pred[n]==1 && y_test[n]==1) tp += 1;
      if(y_pred[n]==0 && y_test[n]==0) tn += 1;
    }
    accuracy = (tp*1.0 + tn)/N_test;
  }  // end local scoping
}
