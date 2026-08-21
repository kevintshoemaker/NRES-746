

# can you tell us about how DHARMa works? How general is it? 

# what about random slopes? What about correlations in random slopes and intercepts? 

# talk more about the syntax for specifying random effects in TMB and lme4

# doesn't lme4 do correlated random effects well?

# do you need to take out the zeros?    

# can you model zero inflation as a function of covariates?

# what is template model builder?

# does glmer.nb fit theta or do you need to specify overdispersion directly.

# can we run a poisson regresion and run DHARMa to see if the fit is rejected?

# any idea when to use nbinom2 vs nbinom1?   


hist(Owls$SiblingNegotiation)


# where to go if we need help?




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

# why is it better to use a mixed model here/

Owls$NegPerChick


library(lme4)

names(Owls)

owls_lme4_nb <- glmer.nb(
  SiblingNegotiation ~ FoodTreatment * SexParent + offset(log(BroodSize)) # This is shorthand for each predictor and the interaction between both
  
  + (1|Nest), # This is the syntax for including Nest as a random effect
  data = Owls
)
summary(owls_lme4_nb)

library(DHARMa)

res = simulateResiduals(owls_lme4_nb)
testResiduals(res)

testZeroInflation(res)   # strong rationale for zero inflation!!!


library(broom.mixed)
library(performance)
tidy(owls_lme4_nb, effects = "fixed", conf.int = TRUE)

r2_nakagawa(owls_lme4_nb)


library(glmmTMB)
owls_tmb_nb1 <- glmmTMB(
  SiblingNegotiation ~ FoodTreatment * SexParent + 
              offset(log(BroodSize)) + (1|Nest),
  data = Owls,
  family = nbinom1()
)
summary(owls_tmb_nb1)

owls_tmb_nb2 <- glmmTMB(
  SiblingNegotiation ~ FoodTreatment * SexParent + 
    offset(log(BroodSize)) + (1|Nest),
  data = owldata,
  family = nbinom2()
)
summary(owls_tmb_nb2)


res_tmb1  <- simulateResiduals(owls_tmb_nb1, plot = FALSE)
plot(res_tmb1)

res_tmb2  <- simulateResiduals(owls_tmb_nb2, plot = FALSE)
plot(res_tmb2)

testZeroInflation(res_tmb1)
testZeroInflation(res_tmb2)

tidy(owls_tmb_nb1, effects = "fixed", conf.int = TRUE)
tidy(owls_tmb_nb2, effects = "fixed", conf.int = TRUE)

r2_nakagawa(owls_tmb_nb1)
r2_nakagawa(owls_tmb_nb2)



# Baseline NB (no ZI factor)
m_nb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                  offset(log(BroodSize)) + (1|Nest),
                data = Owls, family = nbinom2())

# Zero-inflated NB: a ZI component (~ 1) adds extra zeros
m_zinb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                    offset(log(BroodSize)) + (1|Nest),
                  ziformula = ~ FoodTreatment,
                  data = Owls, family = nbinom1())

# Hurdle NB: zero vs positive modeled separately
m_hurd <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                    offset(log(BroodSize)) + (1|Nest),
                  ziformula = ~ FoodTreatment,              # hurdle uses 'truncated_nbinom2' for count
                  data = Owls,
                  family = truncated_nbinom2())


fit = simulateResiduals(m_zinb)

testResiduals(fit)   

summary(m_zinb)

library(effects)
pred_zinb <- allEffects(m_zinb)
plot(pred_zinb, main = "ZI-NB model effects")






m_zinb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                    offset(log(BroodSize)) + ArrivalTime + (ArrivalTime|Nest),
                  ziformula = ~ FoodTreatment,
                  data = Owls, family = nbinom1())

summary(m_zinb)

fit = simulateResiduals(m_zinb)
testResiduals(fit)   


pred_zinb <- allEffects(m_zinb)
plot(pred_zinb, main = "ZI-NB model effects")


### cami and emily - noncentered  -----------

#  start with hierarchical model

   #  nice work on the board-  
   #  generate data from standard normal.

y_global -> y_groupmean -> y_ind

z

# useful to know about.   















