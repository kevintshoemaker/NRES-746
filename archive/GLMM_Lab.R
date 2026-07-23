# Prepare data for analysis with JAGS, including the number of sites
Nsites <- length(levels(as.factor(snakes3$site)))
jagsdata_s3 <- with(snakes3, list(b_mass = b_mass, b_length = b_length, site = site,
                                  N = length(b_mass), Nsites = Nsites))

# JAGS model specification function

filename = "JAGS.txt"

cat("
model{
  # Likelihood:
  for (i in 1:N){
    b_mass[i] ~ dnorm(mu[i], tau)  # Likelihood of the response variable
    mu[i] <- alpha + a[site[i]] + beta * b_length[i]  # Model for the mean
  }
  # Priors:
  alpha ~ dnorm(0, 0.01)      # Prior for the overall intercept
  sigma_a ~ dunif(0, 100)      # Prior for the standard deviation of random effect
  tau_a <- 1 / (sigma_a * sigma_a)  # Convert standard deviation to precision
  for (j in 1:Nsites){
    a[j] ~ dnorm(0, tau_a)    # Prior for random intercept for each site
  }
  beta ~ dnorm(0, 0.01)        # Prior for the slope
  sigma ~ dunif(0, 100)        # Prior for the standard deviation of fixed effect
  tau <- 1 / (sigma * sigma)   # Convert standard deviation to precision
}
",file=filename)
# Initial values function for JAGS
init_values <- function(){
  list(alpha = rnorm(1), sigma_a = runif(1), beta = rnorm(1), sigma = runif(1))
}

# Parameters to be saved from the JAGS model
params <- c("alpha", "beta", "sigma", "sigma_a")

# Fit the JAGS model to the data
fit_lm3 <- jags(data = jagsdata_s3, inits = init_values, parameters.to.save = params, model.file = filename,
                n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 10, DIC = FALSE)

# Display the fitted model
fit_lm3

