// hurdle model with random effects for nest

// data block 
data {
  int<lower=1> N;  // number of observations
  array[N] int<lower=0> y; //array for count response (Sibling negotiation)
  vector[N] food; // predictor: 1 = Deprived, 0 = Satiated
  int<lower=1> J;  //number of nests (grouping factor)
  array[N] int<lower=1, upper=J> nest_id; //nest ID for each observation array
}



// parameters block 
parameters {
  
  //Fixed effects 
  real alpha0;           //intercept for the zero hurdle (logit scale)
  real alpha_food;       //slope for FoodTreatment 
  real beta0;           //intercept for positive counts
  real beta_food;       //slope for FoodTreatment in count model 
  
  //Random effects
  vector[J] a_nest;     //random intercepts for hurdle component
  vector[J] b_nest;     //random intercepts for count compenent
  real<lower=0> sigma_a; //sd for random effects (hurdle compenent)
  real<lower=0> sigma_b; //sd for random effects (count compenent)
  
}

//model block
model {
  //broad priors to help show where it should be looking
  alpha0 ~ normal(0,2);     //for intercept (zero part)
  beta0 ~ normal(0,2);      //for intercept (count part)
  alpha_food ~ normal(0,1);  //effect of food on zero probability
  beta_food ~ normal(0,1); //effect of food on count amount
  a_nest ~ normal(0, sigma_a); //random intercepts for hurdle part
  b_nest ~ normal(0, sigma_b); //random intercepts for count part
  sigma_a ~ exponential(1);   //prior for random effect sd (this keeps it above 0)
  sigma_b ~ exponential(1);   //same thing ^
  
  real theta,lambda;
  
  //actual hurdle model (based on Stan user's guide)
  for(n in 1:N){
    //compute linear predictors (for both hurdle and count)
    //theta = log odds of not calling (hurdle)
    //lamba = expected count rate (based on mean) when there is a call
    theta = inv_logit( alpha0 + alpha_food * food[n] + a_nest[nest_id[n]]); //hurdle predictor
    lambda = exp(beta0 + beta_food * food[n] + b_nest[nest_id[n]]); //count predictor (log)
  
  if(y[n] == 0){
    //for zeros:likelihood is the probability of no call
    target += bernoulli_lpmf(1 | theta);
  } else{
    //for positive counts: 
    target += bernoulli_lpmf(0 | theta)
              + poisson_lpmf(y[n] | lambda)
              - poisson_lccdf(0 | lambda);   // truncated poisson where zero is not allowed
              
         }
    }
}

generated quantities {
  vector [N] log_lik;
  real theta2, lambda2;
  
  
  for(n in 1:N){
    theta2 = inv_logit(alpha0 + alpha_food * food[n] + a_nest[nest_id[n]]);  //hurdle predictor
    lambda2 = exp(beta0 + beta_food * food[n] + b_nest[nest_id[n]]); //count predictor (log)
    if(y[n] == 0){
        //for zeros:likelihood is the probability of no call
      log_lik[n] = bernoulli_lpmf(1 | theta2);
    } else{
        //for positive counts: 
      log_lik[n] = bernoulli_lpmf(0 | theta2)
              + poisson_lpmf(y[n] | lambda2)
              - poisson_lccdf(0 | lambda2);   // truncated poisson where zero is not allowed
              
    }
  
  }
}



