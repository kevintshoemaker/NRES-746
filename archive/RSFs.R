
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #8                    ####
####                                                    ####
############################################################


############################################################
####  RSFs!                                             ####
############################################################



mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
## view the first few rows of the data
head(mydata)

mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")

summary(mylogit)


aod::wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 4:6)


exp(cbind(OR = coef(mylogit), confint(mylogit)))


with(mylogit, null.deviance - deviance) #test statistic 
with(mylogit, df.null - df.residual) #degrees of freedom
with(mylogit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)) #p-value
logLik(mylogit)


library(ResourceSelection)
data(goats)
head(goats)


goats$exp.HLI <- exp(goats$HLI)
goats$sin.SLOPE <- sin(pi * goats$SLOPE / 180)
goats$ELEVATION <- scale(goats$ELEVATION)
goats$ET <- scale(goats$ET)
goats$TASP <- scale(goats$TASP)
head(goats)


m1 <- rspf(STATUS ~ TASP + sin.SLOPE + ELEVATION, goats, m=0, B = 99)
m2 <- rspf(STATUS ~ TASP + ELEVATION, goats, m=0, B = 99)
summary(m1)
summary(m2)
CAIC(m1, m2)


#plot(m1)
mep(m1) # marginal effects similar to plot but with CIs
kdepairs(m1) # 2D kernel density estimates
plot(m2)
kdepairs(m2)
mep(m2)


library(lubridate, quietly = T)
library(amt, quietly = T)
library(raster, quietly = T)
data(deer)
head(deer)


plot(deer$x_, deer$y_)


data("sh_forest")
(sh_forest)
plot(sh_forest)
points(deer$x_, deer$y_)


set.seed(3)
hab <- raster(nrows = sh_forest@nrows, ncol = sh_forest@ncols, ext = sh_forest@extent)
hab <- setValues(hab, rnorm(ncell(hab), mean=15, sd=3))
plot(hab)
points(deer$x_, deer$y_)


avail <- random_points(deer)
plot(avail$x_, avail$y_)
points(deer$x_, deer$y_, col="red", pch=19)


deer <- deer %>% extract_covariates(sh_forest) 
avail <- avail %>% extract_covariates(sh_forest)
deer <- deer %>% extract_covariates(hab) 
avail <- avail %>% extract_covariates(hab)

deer$use <- 1
avail$use <- 0
deer.c <- plyr::rbind.fill(deer, avail)

deer.c <- deer.c %>% 
  mutate(forest = factor(sh.forest, levels = 1:2, labels = c("forest", "non-forest"))) 

head(deer.c)
tail(deer.c)


deer.rsf <- rspf(use ~ forest + layer, data=deer.c, m=0, B=99)
summary(deer.rsf)
mep(deer.rsf) # marginal effects similar to plot but with CIs
kdepairs(deer.rsf) # 2D kernel density estimates
plot(deer.rsf)

