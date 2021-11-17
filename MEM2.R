#_____________________________________#
#     NRES 746, Fall 2021             #
#     Diethelm, McKinnon, Mizell     #
#     University of Nevada, Reno      #
#_____________________________________#


#_____________________________________#
####  Mixed Effects Models  (MEMs) ####
#_____________________________________#

##load packages
library(lme4) 
library(nlme)
library(glmmTMB)
library(ggplot2)
##RIKZ dataset

#as described in Zuur et al. (2007) and Zuur et al. (2009)

##download data to follow along:

#rikz_data <- "https://uoftcoders.github.io/rcourse/data/rikz_data.txt"

#download.file(rikz_data, "rikz_data.txt")

rikz_data  <- read.table("rikz_data.txt", header = TRUE, sep="\t")


rikz_data$Beach <- as.factor(rikz_data$Beach)
str(rikz_data)

##This code is adapted from: https://uoftcoders.github.io/rcourse/lec08-linear-mixed-effects-models.html#accounting_for_non-independence

# Start with a linear model:
lm1 <- lm(Richness~ NAP, data = rikz_data)
summary(lm1)


# Check model assumptions.
par(mfrow=c(2,2))
plot(lm1)

#The QQ plot looks okay, but first panel suggests that the model assumption for
#homogeneity is violated.

#Here we see an increasing variance in the residuals with increasing fitted
#values.

#For now, we will ignore these violations of the model assumptions to explore
#mixed-effects modelling strategies on untransformed data.

#Furthermore, we also know the observations in these data are not independent.

lm2 <- lm(Richness ~ NAP + Beach, data = rikz_data)
summary(lm2)

#Q: What is the influence of NAP on species richness, while accounting for
#variation within beaches?

# Although a Poison might work well here given that richness a count of species,
# we will use a Gaussian distribution to keep things simple

#The (1|Beach) is the random effect term, where the 1 denotes this is a
#random-intercept model and the term on the right of the | is a nominal variable
#(or factor) to be used as the random effect.

#This model is fit using maximum likelihood, rather than restricted maximal
#likelihood by specifying REML = FALSE

#If your data are balanced (i.e., similar sample sizes in each factor group) and
#your random effects are not nested, then you can set REML to FALSE to use
#maximum likelihood.


mem.intercept <- lmer(Richness ~ NAP + (1|Beach), data = rikz_data, REML = FALSE)

summary(mem.intercept)
# Let's plot the linear model for each beach, with a varying intercept only: 
rikz_data$fit_mem.intercept <- predict(mem.intercept)

ggplot(rikz_data, aes(x = NAP, y = Richness, colour = Beach)) +
    # Add fixed effect regression line (i.e. NAP)
    geom_abline(aes(intercept = `(Intercept)`, slope = NAP),
                size = 2,
                as.data.frame(t(fixef(mem.intercept)))) +
    # Add fitted values (i.e. regression) for each beach
    geom_line(aes(y = fit_mem.intercept), size = 1) +
    geom_point(size = 3) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_colour_brewer(palette="Set1")

# Random intercept and slope model
mem_intslope <- lmer(Richness ~ NAP + (1 + NAP|Beach), REML = FALSE,
                             data = rikz_data)

summary(mem_intslope)


# Let's plot the linear model for each beach, with a varying intercept AND slope: 
rikz_data$fit_IntSlope <- predict(mem_intslope)

ggplot(rikz_data, aes(x = NAP, y = Richness, colour = Beach)) +
    geom_abline(aes(intercept = `(Intercept)`, slope = NAP),
                size = 2,
                as.data.frame(t(fixef(mem_intslope)))) +
    geom_line(aes(y = fit_IntSlope), size = 1) +
    geom_point(size = 3) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_colour_brewer(palette="Set1")


mod1 <- glmmTMB(Richness ~ NAP, data = rikz_data)
mod2 <- glmmTMB(Richness ~ NAP + (1|Beach), data = rikz_data)

#The relative fit of two nested models can be evaluated using a chi-square difference statistic:
anova(mod1, mod2, test="Chisq")
#First we want to try scaling our numeric predictor to help with convergence
NAP.sc <- scale(rikz_data$NAP, center = TRUE, scale = TRUE)

#this function calculates the mean and standard deviation of the entire vector,
#then scales each element by those values by subtracting the mean and dividing
#by the sd.

#REML=FALSE is particularly important for linear mixed model selection
mod5 <- lmer(Richness ~ NAP.sc  + (1|Beach), REML=FALSE, data = rikz_data)
summary(mod5)
#AIC = 249.8

mod6 <- lmer(Richness ~ NAP.sc + (1 + NAP|Beach), REML=FALSE, data = rikz_data)
summary(mod6)
#AIC = 246.7

anova(mod5, mod6)

# we still need check model fit
plot(residuals(mod6))
qqnorm(resid(mod6))# QQ-plot
qqline(resid(mod6))

#REML=FALSE does not work in generalized linear mixed model selection
#glmer() uses Maximum Likelihood (ML) as default rather than Restricted Maximum Likelihood (REML)

mod7 <- glmer(Richness ~ NAP.sc  + (1|Beach), family = poisson, data = rikz_data)
summary(mod7)
#AIC = 220.8

mod8 <- glmer(Richness ~ NAP.sc + (1 + NAP|Beach), family = poisson, data = rikz_data)
summary(mod8)
#AIC = 218.7

anova(mod7, mod8)

# we still need check model fit
plot(residuals(mod8))
qqnorm(resid(mod8))# QQ-plot
qqline(resid(mod8))

