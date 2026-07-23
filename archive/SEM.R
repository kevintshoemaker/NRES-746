
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #2                    ####
####                                                    ####
############################################################


############################################################
####  Structural Equation Modeling                      ####
############################################################


## Install packages
#install.packages("OpenMx")
#install.packages("lavaan")
#install.packages("semPlot")
#install.packages("tidyverse")
#install.packages("GGally")

## Load libraries
suppressMessages(library(lavaan))
suppressMessages(library(OpenMx))
suppressMessages(library(semPlot))
suppressMessages(library(tidyverse))
suppressMessages(library(GGally))

## View your data
summary(mtcars)


ggcorr(mtcars[-c(5,7,9)], nbreaks=NULL, label=T, low="red3", high="green3", label_round=2, name="Correlation Scale", label_alpha=T, hjust=0.75) + 
     ggtitle(label="Correlation Plot") + 
     theme(plot.title=element_text(hjust=0.6))


## Simple Linear Model:
linear <- lm(mpg ~ hp+ cyl + disp + carb + am + wt, data=mtcars)

linear


## Set up your model in R
model = '
     # Blue Relationship
     mpg ~ hp + cyl + disp + carb + am + wt
     
     # Green Relationship
     hp ~ cyl + disp + carb
'


## run sem()
path <- sem(model, data=mtcars)

path


summary(path, standardized = T, fit.measures=T, rsquare=T)


## Generated path
semPaths(path, "std", layout='tree', edge.label.cex=.9, curvePivot=T)


#?HolzingerSwineford1939  # some information about the data

head(HolzingerSwineford1939)  # let's look at the first few lines


## Define your model:
ability_model <- '
     # regression equations
     visual ~ sex + ageyr + school + grade
     textual ~ sex + ageyr + school + grade
     speed ~ sex + ageyr + school + grade

     # latent variables
     visual =~ x1 + x2 + x3
     textual =~ x4 + x5 + x6
     speed =~ x7 + x8 + x9

     # covariances
     visual ~~ textual
     visual ~~ speed
     textual ~~ speed
'


## Run the model
fit_ability <- sem(ability_model, data=HolzingerSwineford1939)

## Summarize including fit measures
summary(fit_ability, standardized=T, rsquare=T, fit.measures=T)


## Visualize your paths (may have to try a few!)
semPaths(fit_ability, 'std', layout='circle')

semPaths(fit_ability,'std', layout='tree', edge.label.cex =.9, curvePivot=T)

