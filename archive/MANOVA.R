
############################################################
####                                                    ####  
####  NRES 746, Student Led Topic                       ####
####                                                    ####
####  Maria Sole Bonarota, Mona Farnisa, Ranae Zauner   #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  ANOVA/MANOVA                                      ####
############################################################


library(heplots)
data(RootStock) # Dataset from an experiment comparing apple trees from from six different root stocks from 1918-1934. 

str(RootStock)
library(tidyverse)
library(tibble)
library(ggpubr)
library(rstatix)
library(car)
library(broom)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(purrr)
library(psych)
library(magrittr)
# Load cleaned data

phenotype2.5 <- read.csv("phenotype_data.csv")
phenotype2.5$X <- NULL
# Visualize dataset with boxplots
ggboxplot(phenotype2.5, x = "Variety", y = c("Height_cm", "True_Leaves"), 
  merge = TRUE, palette = "ggplot2")
# Get Summary Statistics
phenotype2.5 %>%
  group_by(Variety) %>%
  get_summary_stats(Height_cm, True_Leaves, type = "mean_sd")

# n in each group should be greater than number of outcome variables 
phenotype2.5 %>%
  group_by(Variety) %>%
  summarise(N = n())
phenotype2.5 %>% 
  group_by(Variety) %>%
  identify_outliers(Height_cm)
# no outliers!

phenotype2.5 %>% 
  group_by(Variety) %>%
  identify_outliers(True_Leaves)
# no outliers!
phenotype2.5 %>%
  group_by(Variety) %>%
  mahalanobis_distance(-id) %>%
  filter(is.outlier == TRUE) %>%
  as.data.frame()
# no outliers!
phenotype2.5 %>%
  group_by(Variety) %>%
  shapiro_test(Height_cm, True_Leaves) %>%
  arrange(variable)
# p-value > 0.05, so data is normally distributed, except for Auto CBD-True Leaves. 
# However, MANOVA is fairly robust to deviations from normality so we will continue
# phenotype2.5 %>%
#   select(Height_cm, True_Leaves) %>%
#   mshapiro_test()


mshapiro_test(phenotype2.5[,2:3])
# p-value greater than 0.05, assumption met!

phenotype2.5 %>% cor_test(Height_cm, True_Leaves)
# Pearson correlation = 0.37
# p-value < 0.05, so no multicolinearity! 
check.linearity <- ggplot(phenotype2.5, aes(x=Height_cm, y=True_Leaves, color=Variety))+
  geom_point(size = 4)+
  theme_classic()+
  geom_smooth(method=lm, se=FALSE) + 
  ylab('True Leaves' ) +
  xlab('Height (cm)') + 
  labs(color = 'Hemp Variety')  #rename legend title
check.linearity
box_m(phenotype2.5[, c("Height_cm", "True_Leaves")], phenotype2.5$Variety)
# p-value < 0.001, so assumption not met
phenotype2.5 %>% 
  gather(key = "variable", value = "value", Height_cm, True_Leaves) %>%
  group_by(variable) %>%
  levene_test(value ~ Variety)
# p-value < 0.05, so no homogeneity of variances!
# True Leaves p-value is 0.05, we will accept, but use Pillais test because it is robust against violations of assumptions
model <- lm(cbind(Height_cm, True_Leaves) ~ Variety, phenotype2.5)
Manova(model, test.statistic = "Pillai")
# Statistically significant difference between Variety on the combined 
# dependent variables (Height and True Leaves), F(4, 54) = 10.546, p < 0.0001.
grouped.data <- phenotype2.5 %>%
  gather(key = "variable", value = "value", Height_cm, True_Leaves) %>%
  group_by(variable)
grouped.data %>% welch_anova_test(value ~ Variety)
# Statistically significant univariate welch_anova difference in Height (F(2, 16.2) = 27.2, p < 0.0001)
# and True Leaves (F(2, 16.7) = 69.9, p < 0.0001 ) between hemp varieties. 
pwc <- phenotype2.5 %>%
  gather(key = "variables", value = "value", Height_cm, True_Leaves) %>%
  group_by(variables) %>%
  games_howell_test(value ~ Variety) %>%
  select(-estimate, -conf.low, -conf.high) # Remove details
pwc
pwc_alt_code <- phenotype2.5 %>%
  gather(key = "variables", value = "value", Height_cm, True_Leaves) %>%
  group_by(variables) %>%
  games_howell_test(value ~ Variety) # %>%
  # select(-estimate, -conf.low, -conf.high) # Remove details
pwc_alt_code$estimate <- NULL
pwc_alt_code$conf.low <- NULL
pwc_alt_code$conf.high <- NULL
pwc_alt_code

# Visualize Results: box plots with p-values
pwc_alt_code <- pwc_alt_code %>% add_xy_position(x = "Variety")
test.label <- create_test_label(
  description = "MANOVA", statistic.text = quote(italic("F")),
  statistic = 10.546, p= "<0.0001", parameter = "4,54",
  type = "expression", detailed = TRUE
)
ggboxplot(
  phenotype2.5, x = "Variety", y = c("Height_cm", "True_Leaves"), 
  merge = TRUE, palette = "ggplot2"
) + 
  stat_pvalue_manual(
    pwc_alt_code, hide.ns = FALSE, tip.length = 0, 
    step.increase = 0.3, step.group.by = "variables",
    color = "variables"
  ) +
  labs(
    subtitle = test.label,
    caption = get_pwc_label(pwc_alt_code, type = "expression")
  )

# load packages
library("car")
library("MASS")

#Load the data
soils <- Soils

# Make variables factors 
Contour <- factor(soils$Contour)
Depth <- factor(soils$Depth)

head(soils)
# Create a model, have to bind the 3 outcome variables using cbind().

soils.mod <- lm(cbind(pH, Dens, Conduc) ~ Contour + Depth + Contour*Depth -1, data = soils) # Note: -1 removes the intercept effect
summary(soils.mod)

# Note: Manova() is brought to you by the car package

Manova(soils.mod, multivariate=T, type= c("II"), test=("Wilks"))
# Manova(soils.mod, multivariate=T, type= c("II"), test=("Pillai"))
# Manova(soils.mod, multivariate=T, type= c("II"), test=("Hotelling-Lawley"))
# Manova(soils.mod, multivariate=T, type= c("II"), test=("Roy"))

### If we did have an interaction term we would use type III 

manova(soils.mod, multivariate=T, type= c("III"), test=("wilks"))
# manova(soils.mod, multivariate=T, type= c("III"), test=("Pillai"))
# manova(soils.mod, multivariate=T, type= c("III"), test=("Hotelling-Lawley"))
# manova(soils.mod, multivariate=T, type= c("III"), test=("Roy"))
library("heplots")
etasq(Anova(soils.mod),anova=TRUE)

etasq(soils.mod,test="Wilks") # Using Wilks to be consistent with above
# etasq(soils.mod,test="Hotelling")
# etasq(soils.mod,test="Roy")
# etasq(soils.mod,test="Pillai")
