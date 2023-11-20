# Exercise 1

df <- read.csv("Exercise_1.csv")
df$id <- as.factor(as.numeric(as.factor(df$Individual_ID)))
df$Season <- as.factor(df$Season)

table(df$id)

df <- df[,-1]

plot(df$Hemoglobin~df$Season)

hist(df$Hemoglobin)

mod1 <- lm(Hemoglobin~Season,df)
summary(mod1)

hist(mod1$residuals)  # looks normally distributed...

library(glmmTMB)


mod2 <- glmmTMB(Hemoglobin~Season+(1|id),family=gaussian(),data=df)

AIC(mod1,mod2)   # simpler model is better


library(DHARMa)

res <- simulateResiduals(mod1)  # violation of homogeneity
plot(res)

res <- simulateResiduals(mod2)
plot(res)    # still a violation



# Exercise 2

df <- read.csv("beesurv.csv")

table(df$site)  # 10 sites, 10 observations

hist(df$visits)   # definitely non-normal! no zeroes... 

names(df)

plot(visits~temp,df)
plot(visits~wind,df)

df$site <- as.factor(df$site)
mod1 <- glmmTMB(visits~temp+wind+((temp+wind)|site),
                family = Gamma(link="log"),data=df)

mod2 <- glmmTMB(visits~temp+wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = Gamma(link="log"),data=df)


summary(mod2)

ranef(mod2)

res <- simulateResiduals(mod2)
plot(res)   # some issues, but not terrible

df$visits


mod3 <- glmmTMB(visits~temp+wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = nbinom1(link="log"),data=df)


summary(mod3)  # not converged

res <- simulateResiduals(mod3)
plot(res)   # some issues, but not terrible


mod4 <- glmmTMB(log(visits)~temp+wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = gaussian(link="identity"),data=df)


summary(mod4)  # not converged

res <- simulateResiduals(mod4)
plot(res)   # some issues, but not terrible  Mod 2 is still the best


mod5 <- glmmTMB(visits~temp*wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = Gamma(link="log"),data=df)


summary(mod5)

ranef(mod5)

res <- simulateResiduals(mod5)
plot(res)   # some issues, but not terrible

df$visits

AIC(mod1,mod2,mod3,mod5) # AIC says mod3 is best!


mod6 <- glmmTMB(visits~temp+wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = genpois(link="log"),data=df)

summary(mod6)
res <- simulateResiduals(mod6)
plot(res)   # some issues, but not terrible


AIC(mod2,mod3,mod6)   # mod 3 and 6 are identical...



mod7 <- glmmTMB(visits~temp*wind+(1|site) + ((0+temp)|site) + ((0+wind)|site),
                family = poisson(link="log"),data=df)

summary(mod7)
res <- simulateResiduals(mod7)
plot(res)   # some issues, but not terrible


AIC(mod2,mod3,mod6,mod7)   # mod 3 and 6 are identical...
    # poisson model fits okay





# Exercise 3 (bayesian)

set.seed(123)

# Define parameters
samplesize <- 200         # Number of observations
nsites <- 10              # Number of sites
b_length <- sort(rnorm(samplesize))  # Explanatory variable (body length)

# Generate a random grouping variable 'sites' with replacement
sites <- sample(1:10, samplesize, replace = TRUE)

# Display the count of observations per site
table(sites)


# True intercept parameters
int_true_mean <- 45       # True mean intercept
int_true_sigma <- 10      # True standard deviation of intercepts
int_true_sites <- rnorm(n = nsites, mean = int_true_mean, sd = int_true_sigma)  # True intercept of each site

# Create a matrix to represent the intercept of each snake individual based on its site
sitemat <- matrix(0, nrow = samplesize, ncol = nsites)
for (i in 1:nrow(sitemat)) sitemat[i, sites[i]] <- 1
int_true <- sitemat %*% int_true_sites

# True slope parameter
slope_true <- 10

# Calculate true means and standard deviation
mu <- int_true + slope_true * b_length
sigma <- 5

# Generate response variable 'b_mass' based on normal distributions
b_mass <- rnorm(samplesize, mean = mu, sd = sigma)

# Create a data frame 'snakes3' with explanatory and response variables
snakes3 <- data.frame(b_length = b_length, b_mass = b_mass, site = sites)

# Display the first few rows of the data frame
head(snakes3)


plot(b_mass ~ b_length, col = site, data = snakes3)


# True intercept parameters
int_true_mean <- 45       # True mean intercept
int_true_sigma <- 10      # True standard deviation of intercepts
int_true_sites <- rnorm(n = nsites, mean = int_true_mean, sd = int_true_sigma)  # True intercept of each site

# Create a matrix to represent the intercept of each snake individual based on its site
sitemat <- matrix(0, nrow = samplesize, ncol = nsites)
for (i in 1:nrow(sitemat)) sitemat[i, sites[i]] <- 1
int_true <- sitemat %*% int_true_sites

# True slope parameter
slope_true <- 10

# Calculate true means and standard deviation
mu <- int_true + slope_true * b_length
sigma <- 5

# Generate response variable 'b_mass' based on normal distributions
b_mass <- rnorm(samplesize, mean = mu, sd = sigma)

# Create a data frame 'snakes3' with explanatory and response variables
snakes3 <- data.frame(b_length = b_length, b_mass = b_mass, site = sites)

# Display the first few rows of the data frame
head(snakes3)

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

library(jagsUI)

# Fit the JAGS model to the data
fit_lm3 <- jags(data = jagsdata_s3, inits = init_values, parameters.to.save = params, model.file = filename,
                n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 10, DIC = FALSE)

fit_lm3

plot(fit_lm3)

coda::autocorr(fit_lm3$samples)  # some strong autocorrelation

traceplot(fit_lm3,"beta")


?autocorr




