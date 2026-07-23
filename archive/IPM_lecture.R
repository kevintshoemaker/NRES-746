library(IPMbook); library(jagsUI)
data(wryneck)
str(wryneck)

tail(wryneck)
# Identify failed broods
fail <- which(wryneck$x==0)
# Create encounter histories
y <- matrix(NA, nrow=length(wryneck$f), ncol=max(wryneck$k))
for (i in 1:length(wryneck$f)){
y[i,wryneck$f[i]] <- 1
y[i,wryneck$j[i]] <- 1
}
for (i in 1:length(fail)){
y[fail[i],wryneck$k[fail[i]]] <- 0
}
y[176,]
y[177,]
# Bundle data
jags.data <- with(wryneck, list(y=y, f=f, k=k, n.nest=nrow(y), T=21, age=age))
str(jags.data)
# Write JAGS model file
cat(file="model13.txt", "
model {
# Priors and linear models
for (i in 1:n.nest){
for (t in f[i]:(k[i]-1)){
phi[i,t] <- phia[age[i] + t - f[i]]
} #t
} #i
for (a in 1:T){
phia[a] <- ilogit(alpha + beta * a)
}
alpha ~ dnorm(0, 0.001)
beta ~ dnorm(0, 0.001)
# Likelihood
for (i in 1:n.nest){
for (t in (f[i]+1):k[i]){
y[i,t] ~ dbern(phi[i,t-1] * y[i,t-1])
} #t
} #i
# Derived parameter: nest success
nu <- prod(phia[1:T])
}
")
# Initial values
inits <- function(){list(alpha=runif(1, 4, 5), beta=runif(1, 0, 0.1))}
# Parameters monitored
parameters <- c("phia", "nu", "alpha", "beta")
# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000
# Call JAGS from R (ART 1 min) and check convergence
out16 <- jags(jags.data, inits, parameters, "model13.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
n.thin=nt, n.adapt=na, parallel=TRUE)
#traceplot(out16) # Not shown
print(out16, 3)
# Choose constants in simulation
nbrood <- 1000 # Number of broods with young counted
brood.mean <- 1.5 # Average brood size
sd.brood <- 0.3 # log-linear brood random effect
# Draw Poisson random numbers
set.seed(24)
expNyoung <- exp(log(brood.mean) + rnorm(nbrood, 0, sd.brood))
C <- rpois(nbrood, expNyoung)
table(C)
# Data bundle
jags.data <- list(C1=C, sumC1=sum(C), C1copy=C, nbrood=nbrood)
str(jags.data)

library(jagsUI)

# Write JAGS model file
cat(file="model8.txt", "
model {
# Priors and linear models
rho1 ~ dunif(0, 5) # Mean brood size in model 1
rho2 ~ dunif(0, 5) # Mean brood size in model 2
rho3 ~ dunif(0, 5) # Mean brood size in model 3
tau.rho3 <- pow(sd.rho3, -2)
sd.rho3 ~ dunif(0, 3) # Brood-level overdispersion in model 3
# Likelihoods for three separate models
# Model 1: Poisson GLM for disaggregated data
for (i in 1:nbrood){
C1[i] ~ dpois(rho1)
}
# Model 2: Poisson GLM for aggregated data
sumC1 ~ dpois(rho2 * nbrood)
# Model 3: Poisson GLMM for aggregated data with brood-level overdispersion
for (i in 1:nbrood){
C1copy[i] ~ dpois(pois.mean[i])
log(pois.mean[i]) <- logmean[i]
logmean[i] ~ dnorm(log(rho3), tau.rho3)
}
}
")
# Initial values
inits <- function(){list(rho1=runif(1, 0.5, 2.5), rho2=runif(1, 0.5, 2.5),
rho3=runif(1, 0.5, 2.5))}
# Parameters monitored
parameters <- c("rho1", "rho2", "rho3", "sd.rho3")
# MCMC settings
ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000
# Call JAGS from R (ART 3 min), check convergence and summarize posteriors
out11 <- jags(jags.data, inits, parameters, "model8.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
n.thin=nt, n.adapt=na, parallel=TRUE)
#traceplot(out11) # Not shown
print(out11, 3)


# Choose values for data simulation
nmarked <- 10 # Number of marked individuals at each occasion
nyears <- 11 # Number of years
phi <- 0.8 # Constant apparent survival probability
p <- 0.4 # Constant recapture probability
# Determine occasion when an individual first captured and marked
f <- rep(1:(nyears-1), each=nmarked)
nind <- length(f) # Total number of marked individuals
# State or ecological process
z <- array(NA, dim=c(nind, nyears)) # Empty alive/dead matrix
# Initial conditions: all individuals alive at f(i)
for (i in 1:nind){
z[i,f[i]] <- 1
}
set.seed(1) # Initialize the RNGs in R
# Propagate alive/dead process forwards via transition rule:
# Alive individuals survive with probability phi
for (i in 1:nind){
for (t in (f[i]+1):nyears){
z[i,t] <- rbinom(1, 1, z[i,t-1] * phi)
} #t
} #i
head(z); tail(z) # look at start and end of z
# Observation process: simulate observations
y <- array(0, dim=c(nind, nyears))
for (i in 1:nind){
y[i,f[i]] <- 1
for(t in (f[i]+1):nyears){
y[i,t] <- rbinom(1, 1, z[i,t] * p)
} #t
} #i
head(y) # Complete simulated capture-recapture data set 
for (i in 1:10){ # Look at true and observed states of first 10 individuals
print(rbind("True state (z)" = z[i,], "Observed state (y)" = y[i,]))
browser()
}
# Data bundle
jags.data <- list(y=y, f=f, nind=nind, nyears=ncol(y))
str(jags.data)
# Write JAGS model file
cat(file="model14.txt", "
model {
# Priors and linear models
phi.const ~ dunif(0, 1) # Vague prior for constant phi
p.const ~ dunif(0, 1) # Vague prior for constant p
for (i in 1:nind){ # Loop over individuals
for (t in f[i]:(nyears-1)){ # Loop over time intervals/occasions
phi[i,t] <- phi.const # Here we model pattern in phi ...
p[i,t] <- p.const # ... and p
} #t
} #i
# Likelihood
for (i in 1:nind){
# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):nyears){
# State process
z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
# Observation process
y[i,t] ~ dbern(z[i,t] * p[i,t-1])
} #t
} #i
}
")

# Initial values
inits <- function(){list(z=zInit(y))}

# Parameters monitored
parameters <- c("phi.const", "p.const") # Could also add "z"
# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000
# Call JAGS from R (ART < 1 min) and check convergence
out17 <- jags(jags.data, inits, parameters, "model14.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
n.thin=nt, n.adapt=na, parallel=TRUE)
#traceplot(out17) # Not shown
print(out17, 3)


# Chose constants
nyears <- 25            # Number of years
N1 <- 30                # Initial population size, t=1 
mu.lam <- 1.02          # Mean annual population growth rate
sig2.lam <- 0.02        # Process (temporal) variation of the growth rate
sig2.y <- 400           # Variance of the observation error

# Simulate true system state
N <- numeric(nyears)
N[1] <- N1                     # Set initial abundance
set.seed(1)
lambda <- rnorm(nyears-1, mu.lam,sqrt(sig2.lam))   # Draw random lambdas
for(t in 1:(nyears-1)){     
  N[t+1] <- N[t]*lambda[t]     # Propagate population size forward
}

# Simulate observations
eps <- rnorm(nyears,0,sqrt(sig2.y))   # Draw random residuals 
y <- N + eps                   # Add residual error to value of true state


library(IPMbook); library(jagsUI)
# Data bundle
jags.data <- list(y=y, T=length(y))
str(jags.data)



# Write JAGS model file
cat(file="model1.txt", "
model {
# Priors and linear models
mu.lam ~ dunif(0, 10) # Prior for mean growth rate
sig.lam ~ dunif(0, 1) # Prior for sd of growth rate
sig2.lam <- pow(sig.lam, 2)
tau.lam <- pow(sig.lam, -2)
sig.y ~ dunif(0.1, 100) # Prior for sd of observation process
sig2.y <- pow(sig.y, 2)
tau.y <- pow(sig.y, -2)
# Likelihood
# Model for the initial population size: uniform priors
N[1] ~ dunif(0, 500)
# Process model over time: our model of population dynamics
for (t in 1:(T-1)){
lambda[t] ~ dnorm(mu.lam, tau.lam)
N[t+1] <- N[t] * lambda[t]
}
# Observation process
for (t in 1:T){
y[t] ~ dnorm(N[t], tau.y)
}
}
")
# Initial values
inits <- function(){list(sig.lam=runif(1, 0, 1))}
# Parameters monitored
parameters <- c("lambda", "mu.lam", "sig2.y", "sig2.lam", "sig.y", "sig.lam", "N")
# MCMC settings
ni <- 200000; nb <- 10000; nc <- 3; nt <- 100; na <- 5000
# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
n.thin=nt, n.adapt=na, parallel=TRUE)
#traceplot(out1) # Not shown
print(out1, 3)
