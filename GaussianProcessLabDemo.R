# GP data
GP_data <- read.csv("GP_data.csv")
names(GP_data) <- c("x", "year", "births", "UE.rate", "MH.income")
# Visualize the births per year
mean_births <- mean(GP_data$births)
plot(
  GP_data$year,
  GP_data$births,
  type = "b",
  pch = 16,
  col = "blue",
  lwd = 2,
  xlab = "Year",
  ylab = "Total Births",
  main = "Total Births per Year with Mean Line"
)
abline(
  h = mean_births,
  col = "red",
  lwd = 2,
  lty = 2
)
legend(
  "topleft",
  legend = c("Total births", "Mean births"),
  col = c("blue", "red"),
  lwd = c(2, 2),
  lty = c(1, 2),
  pch = c(16, NA),
  bty = "n"
)
# Visualize UN rates
mean_UN <- mean(GP_data$UE.rate)
plot(
  GP_data$year,
  GP_data$UE.rate,
  type = "b",
  pch = 16,
  col = "blue",
  lwd = 2,
  xlab = "Year",
  ylab = "Average UN rate",
  main = "Average Unemployment Rates per Year"
)
abline(
  h = mean_UN,
  col = "red",
  lwd = 2,
  lty = 2
)
legend(
  "topleft",
  legend = c("Average UN rates", "Mean UN rate"),
  col = c("blue", "red"),
  lwd = c(2, 2),
  lty = c(1, 2),
  pch = c(16, NA),
  bty = "n"
)
# Visualize Median Household income
mean_MH <- mean(GP_data$MH.income)
plot(
  GP_data$year,
  GP_data$MH.income,
  type = "b",
  pch = 16,
  col = "blue",
  lwd = 2,
  xlab = "Year",
  ylab = "Median Household Income",
  main = "Median Household Income per Year"
)
abline(
  h = mean_MH,
  col = "red",
  lwd = 2,
  lty = 2
)
legend(
  "topleft",
  legend = c("Median Household Income", "Mean MH"),
  col = c("blue", "red"),
  lwd = c(2, 2),
  lty = c(1, 2),
  pch = c(16, NA),
  bty = "n"
)
#### Prepare data for models ####
# Ensure numeric columns
GP_data$year <- as.numeric(GP_data$year)
GP_data$births <- as.numeric(GP_data$births)
GP_data$UE.rate <- as.numeric(GP_data$UE.rate)
GP_data$MH.income <- as.numeric(GP_data$MH.income)
# Scale births for numeric stability
GP_data$births_scaled <- as.numeric(scale(GP_data$births))
GP_data$year_scaled <- as.numeric(scale(GP_data$year))
GP_data$UE_scaled <- as.numeric(scale(GP_data$UE.rate))
GP_data$MH_scaled <- as.numeric(scale(GP_data$MH.income))
                             
   
#### Check for colinearity ####
corDF <- data.frame(UE_rate = GP_data$UE.rate, MH_income = GP_data$MH.income, Year = GP_data$year)
cor(corDF)


#### Set priors for model ####
# See what brms would give
get_prior(births_scaled ~ gp(year_scaled + MH_scaled),
          data = GP_data, family = gaussian())
# How to set your own priors
priors <- c(
  prior(normal(0, 5), class = "sdgp"),       # marginal SD of GP
  prior(exponential(1), class = "lscale"),   # length-scale
  prior(exponential(1), class = "sigma"),
  prior(normal(0, 10), class = "Intercept")
)


#### Specify the different models
# Gaussian Process Model
GP_model <- brm(
  births_scaled ~ gp(year_scaled, MH_scaled), #cov = exp_quad,
  data = GP_data,
  family = gaussian(),
  #prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000
)
# Bayesian GLM model
BRMS_glm <- brm(
  births_scaled ~ year_scaled + MH_scaled,
  data = GP_data,
  family = gaussian(),
  chains = 4,
  cores = 4,
  iter = 2000
)
# GAM with splines
BRMS_gam <- brm(
  births_scaled ~ s(year_scaled) + s(MH_scaled),
  data = GP_data,
  family = gaussian(),
  chains = 4,
  cores = 4,
  iter = 2000
)


#### Loo Model Comparison ####
# Get loo scores
loo_gp  <- loo(GP_model)
loo_glm <- loo(BRMS_glm)
loo_gam <- loo(BRMS_gam)
# Loo table
print(loo_compare(loo_gp, loo_glm, loo_gam))


#### Model summaries ####
summary(GP_model)
summary(BRMS_glm)
summary(BRMS_gam)
#### Trace plots ####
plot(GP_model)
mcmc_plot(GP_model, type = "trace")


#### Posterior draws ####
pp_check(GP_model, ndraws = 50)


#### Prior predictive check ####
# Fit prior predictive model
prior_mod <- brm(
  births_scaled ~ gp(year_scaled, MH_scaled),
  data = GP_data,
  family = gaussian(),
  prior = priors,
  sample_prior = "only",
  chains = 4,
  iter = 2000
)
# Prior predictive check
pp_check(prior_mod)


#### Coefficient effect plots and model fit plot ####
conditional_effects(GP_model)
conditional_effects(BRMS_glm)
conditional_effects(BRMS_gam)


#### Get and look at stan code ####
GP_stan <- stancode(GP_model)
GP_stan
# Write out file
#write_stan_file(GP_stan, basename = "GP_stan.stan")
