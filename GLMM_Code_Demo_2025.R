#install.packages("lme4")
#install.packages("glmmTMB")
#install.packages("ggplot2")
#install.packages("broom.mixed")
#install.packages("performance")
#install.packages("DHARMa")
#install.packages("effects")
library(glmmTMB)
data("Owls")
head(Owls)
library(ggplot2)

ggplot(Owls, aes(SiblingNegotiation)) +
       geom_histogram(binwidth = 1) +
       labs(title = "Begging calls per visit", x = "Calls", y = "Frequency")
owldata <- subset(Owls, SiblingNegotiation > 0) # remove zeros
summary(owldata$SiblingNegotiation) # double check distribution

ggplot(owldata, aes(SiblingNegotiation)) +
       geom_histogram(binwidth = 1) +
       labs(title = "Begging calls per visit", x = "Calls", y = "Frequency")
var(owldata$SiblingNegotiation) # Var(Y)
mean(owldata$SiblingNegotiation) # E(Y)
library(lme4)
owls_lme4_nb <- glmer.nb(
SiblingNegotiation ~ FoodTreatment * SexParent # This is shorthand for each predictor and the interaction between both

+ (1|Nest), # This is the syntax for including Nest as a random effect
data = owldata
)
summary(owls_lme4_nb)
library(DHARMa)
par(mfrow = c(2, 2))
res <- simulateResiduals(owls_lme4_nb, plot = FALSE)
plot(res)
library(broom.mixed)
library(performance)
tidy(owls_lme4_nb, effects = "fixed", conf.int = TRUE)
r2_nakagawa(owls_lme4_nb)
library(glmmTMB)
owls_tmb_nb1 <- glmmTMB(
SiblingNegotiation ~ FoodTreatment * SexParent + (1|Nest),
data = owldata,
family = nbinom1()
)
summary(owls_tmb_nb1)

owls_tmb_nb2 <- glmmTMB(
SiblingNegotiation ~ FoodTreatment * SexParent + (1|Nest),
data = owldata,
family = nbinom2()
)
summary(owls_tmb_nb2)
res_tmb1  <- simulateResiduals(owls_tmb_nb1, plot = FALSE)
plot(res_tmb1)

res_tmb2  <- simulateResiduals(owls_tmb_nb2, plot = FALSE)
plot(res_tmb2)
tidy(owls_tmb_nb1, effects = "fixed", conf.int = TRUE)
tidy(owls_tmb_nb2, effects = "fixed", conf.int = TRUE)
r2_nakagawa(owls_tmb_nb1)
r2_nakagawa(owls_tmb_nb2)
AIC(owls_lme4_nb, owls_tmb_nb1, owls_tmb_nb2)
r2_nakagawa(owls_lme4_nb)
r2_nakagawa(owls_tmb_nb1)
r2_nakagawa(owls_tmb_nb2)
logLik(owls_lme4_nb)
logLik(owls_tmb_nb1)
logLik(owls_tmb_nb2)
# Baseline NB (no ZI factor)
m_nb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + (1|Nest),
data = Owls, family = nbinom2())

# Zero-inflated NB: a ZI component (~ 1) adds extra zeros
m_zinb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + (1|Nest),
ziformula = ~ 1,
data = Owls, family = nbinom2())

# Hurdle NB: zero vs positive modeled separately
m_hurd <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + (1|Nest),
ziformula = ~ 1,              # hurdle uses 'truncated_nbinom2' for count
data = Owls,
family = truncated_nbinom2())
AIC(m_nb, m_zinb, m_hurd)
summary(m_zinb)$coefficients
library(effects)
pred_zinb <- allEffects(m_zinb)
plot(pred_zinb, main = "ZI-NB model effects")

