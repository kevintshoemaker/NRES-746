### RIKZ dataset-----------------
#as described in Zuur et al. (2007) and Zuur et al. (2009)

##download data to follow along:
rikz_data <- "https://uoftcoders.github.io/rcourse/data/rikz_data.txt"
download.file(rikz_data, "rikz_data.txt")
rikz_data <- read.table("rikz_data.txt", header = TRUE, sep="\t")

### Data exploration-------------------------

str(rikz_data)
rikz_data$Beach <- as.factor(rikz_data$Beach)
str(rikz_data)
head(rikz_data)

### Run basic linear model using all of the data--------------------

basic.lm <- lm(Richness~ NAP, data = rikz_data)
summary(basic.lm)

library(ggplot2)

# Plot relationship from above model
ggplot(rikz_data, aes(x = NAP, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# Check assumptions.
par(mfrow=c(2,2))
plot(basic.lm)

# Function to find polygons
find_hull <- function(df) df[chull(df$Richness, df$NAP), ]

# Identify polygons in data
library(plyr)
hulls <- ddply(rikz_data, "Beach", find_hull)

# Plot
ggplot(rikz_data, aes(x = NAP, y = Richness, colour = Beach)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  geom_polygon(data=hulls, aes(fill = Beach), alpha = 0.2)


basic.lm <- lm(Richness ~ NAP + Beach, data = rikz_data)
summary(basic.lm)


# Random intercept model with NAP as fixed effect and Beach as random effect
library(lme4)
mixed_model_IntOnly <- lmer(Richness ~ NAP + (1|Beach),
                            data = rikz_data, REML = FALSE)
summary(mixed_model_IntOnly)

# Let's predict values based on our model and add these to our dataframe
# These are the fitted values for each beach, which are modelled separately.
rikz_data$fit_InterceptOnly <- predict(mixed_model_IntOnly)

# Let's plot
ggplot(rikz_data, aes(x = NAP, y = Richness, colour = Beach)) +
  # Add fixed effect regression line (i.e. NAP)
  geom_abline(aes(intercept = `(Intercept)`, slope = NAP),
              linewidth = 2,
              as.data.frame(t(fixef(mixed_model_IntOnly)))) +
  # Add fitted values (i.e. regression) for each beach
  geom_line(aes(y = fit_InterceptOnly), size = 1) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_brewer(palette="Set1")

# Random intercept and slope model
mixed_model_IntSlope <- lmer(Richness ~ NAP + (1 + NAP|Beach),
                             data = rikz_data, REML = FALSE)
summary(mixed_model_IntSlope)

rikz_data$fit_IntSlope <- predict(mixed_model_IntSlope)
ggplot(rikz_data, aes(x = NAP, y = Richness, colour = Beach)) +
  geom_abline(aes(intercept = `(Intercept)`, slope = NAP),
              size = 2,
              as.data.frame(t(fixef(mixed_model_IntSlope)))) +
  geom_line(aes(y = fit_IntSlope), size = 1) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_brewer(palette="Set1")

mixed_model_NoFix <- lmer(Richness ~ 1 + (1|Beach),
                          data = rikz_data, REML = TRUE)
summary(mixed_model_NoFix)

# Define 2 models. Fit both with REML.

mixed_model_IntOnly <- lmer(Richness ~ NAP*Exposure + (1|Beach), REML = TRUE, 
                            data = rikz_data)
mixed_model_IntSlope <- lmer(Richness ~ NAP*Exposure + (1 + NAP|Beach), REML = TRUE, 
                             data = rikz_data)
library(MuMIn)
AICc(mixed_model_IntOnly, mixed_model_IntSlope)

# Full model with both fixed effects and their interaction
mixed_model_IntOnly_Full <- lmer(Richness ~ NAP*Exposure + (1|Beach), REML = FALSE, 
                                 data = rikz_data)

# No interaction
mixed_model_IntOnly_NoInter <- lmer(Richness ~ NAP + Exposure + (1|Beach), 
                                    REML = FALSE, 
                                    data = rikz_data)

# No interaction or main effect of exposure
mixed_model_IntOnly_NAP <- lmer(Richness ~ NAP + (1|Beach), 
                                REML = FALSE, 
                                data = rikz_data)

# No interaction or main effect of NAP
mixed_model_IntOnly_Exp <- lmer(Richness ~ Exposure + (1|Beach), 
                                REML = FALSE, 
                                data = rikz_data)

# No fixed effects
mixed_model_IntOnly_NoFix <- lmer(Richness ~ 1 + (1|Beach), 
                                  REML = FALSE, 
                                  data = rikz_data)


AICc(mixed_model_IntOnly_Full, mixed_model_IntOnly_NoInter,
     mixed_model_IntOnly_NAP, mixed_model_IntOnly_Exp,
     mixed_model_IntOnly_NoFix)

# Summarize best-fit model

summary(update(mixed_model_IntOnly_Full, REML = TRUE))

