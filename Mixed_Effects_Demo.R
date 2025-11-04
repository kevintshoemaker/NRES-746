library(lme4)
library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom.mixed)
set.seed(123)
n_schools <- 5
n_students <- 50

df_nested <- data.frame(
  school = factor(rep(1:n_schools, each = n_students)),
  treatment = factor(rep(rep(c("control","treated"), each = n_students/2), times = n_schools))
)

# Random intercepts and slopes per school
school_intercepts <- rnorm(n_schools, 0, 1)
school_slopes <- rnorm(n_schools, 2, 0.5)

df_nested$y <- 10 +
  school_intercepts[as.numeric(df_nested$school)] +
  school_slopes[as.numeric(df_nested$school)] * (df_nested$treatment == "treated") +
  rnorm(nrow(df_nested), 0, 1)

head(df_nested)
str(df_nested)
m_nested <- lmer(y ~ treatment + (1 + treatment | school), data = df_nested)
summary(m_nested)
fix <- fixef(m_nested)
ranefs <- ranef(m_nested)$school

school_effects <- data.frame(
  school = rownames(ranefs),
  intercept = fix["(Intercept)"] + ranefs[,"(Intercept)"],
  treatment_effect = fix["treatmenttreated"] + ranefs[,"treatmenttreated"]
)

school_effects
ggplot(school_effects, aes(x = school, y = treatment_effect)) +
  geom_point(size = 3) +
  geom_hline(yintercept = fix["treatmenttreated"], linetype = "dashed", color = "red") +
  labs(y = "Estimated treatment effect", title = "School-specific treatment effects") +
  theme_minimal()
df_nested$pred <- predict(m_nested)

ggplot(df_nested, aes(x = treatment, y = y, group = school, color = school)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.3)) +
  stat_summary(fun = mean, geom = "line", aes(group = school), size = 1) +
  labs(y = "Outcome", title = "Varying treatment effects by school") +
  theme_minimal()
vc <- VarCorr(m_nested)
vc
attr(vc$school, "correlation")[1, 2]
set.seed(123)
n_students <- 40
n_raters <- 6

df_crossed <- expand.grid(
  student = factor(1:n_students),
  rater = factor(1:n_raters),
  treatment = factor(c("control","treated"))
)

rater_intercepts <- rnorm(n_raters, 0, 2)
rater_slopes <- rnorm(n_raters, 3, 0.7)

df_crossed$score <- 50 +
  rater_intercepts[as.numeric(df_crossed$rater)] +
  rater_slopes[as.numeric(df_crossed$rater)] * (df_crossed$treatment == "treated") +
  rnorm(nrow(df_crossed), 0, 2)

head(df_crossed)
str(df_crossed)
m_crossed <- lmer(score ~ treatment + (1 + treatment | rater) + (1 | student), data = df_crossed)
summary(m_crossed)
fix <- fixef(m_crossed)
ranefs <- ranef(m_crossed)$rater

rater_effects <- data.frame(
  rater = rownames(ranefs),
  intercept = fix["(Intercept)"] + ranefs[,"(Intercept)"],
  treatment_effect = fix["treatmenttreated"] + ranefs[,"treatmenttreated"]
)

rater_effects
ggplot(rater_effects, aes(x = rater, y = treatment_effect)) +
  geom_point(size = 3) +
  geom_hline(yintercept = fix["treatmenttreated"], linetype = "dashed", color = "red") +
  labs(y = "Estimated treatment effect", title = "Rater-specific treatment effects") +
  theme_minimal()
df_crossed$pred <- predict(m_crossed)

ggplot(df_crossed, aes(x = treatment, y = pred, group = rater, color = rater)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.3)) +
  stat_summary(fun = mean, geom = "line", aes(group = rater), size = 1) +
  labs(y = "Predicted score", title = "Varying treatment effects by rater") +
  theme_minimal()
VarCorr(m_crossed)
data("sleepstudy")
head(sleepstudy)
