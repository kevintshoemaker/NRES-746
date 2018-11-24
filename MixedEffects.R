
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #5                    ####
####                                                    ####
############################################################


############################################################
####  Mixed Effects Models                              ####
############################################################



library(Matrix)
library(lme4)
library(MASS)
library(arm)
library(sjstats)
library(ResourceSelection)


# load in
df.rikz <- read.csv("RIKZ.csv")

# test for correlation
abs(cor(df.rikz[,c(2,5)]))

# define the model
lme.rikz <- lmer(Richness ~ NAP + (1 | Beach), data = df.rikz)

# residual vs fitted plot
plot(lme.rikz, xlab = "Fitted", ylab = "Residual")


summary(lme.rikz)


rm(list=ls())

df.sloth <- read.csv("sloth_point5only_withPatchGroup.csv", header = TRUE, sep = ",")
str(df.sloth) # 25 variables with 25 observations
# drop metadata, response, random effect, categorical data
df.sloth2 <- subset.data.frame(df.sloth, select = -c(SITE_ID, LONGITUDE, LATITUDE, LAND_USE, patch_grou, OCCURRENCE))
df.sloth2$PR <- as.numeric(df.sloth2$PR)
str(df.sloth2)


IsItNormal <- function(df, save = F) { ### should be able to accept any dataframe of arbitrary size & names, if all columns contain data that Shapiro-Wilk test likes
  slothtest <- list()
  for(i in seq(dim(df)[2])) { 
    slothtest <- append(slothtest, shapiro.test(df[[i]]))
  }
  mx.slothtest <- matrix(slothtest, byrow = T, ncol = 4)
  df.slothtest <- data.frame(var = as.character(colnames(df)),
                             W = as.numeric(mx.slothtest[,1]),
                             p = as.numeric(mx.slothtest[,2]),
                             row.names = NULL)
  if(save == T) {write.table(df.slothtest, file = paste0("IsItNormal", round(as.numeric(Sys.time())), ".csv"), row.names = F, sep = ",")}
  return(df.slothtest)
}

IsItNormal(df.sloth2)

df.sloth.log <- as.data.frame(apply(df.sloth2, 2, log))
df.sloth.sqrt <- as.data.frame(apply(df.sloth2, 2, sqrt))
df.sloth.squared <- as.data.frame(apply(df.sloth2, 2, function(x) x^2))

df.slothtest.transformed <- cbind(IsItNormal(df.sloth2),
                                  IsItNormal(df.sloth.log)[,2:3],
                                  IsItNormal(df.sloth.sqrt)[,2:3],
                                  IsItNormal(df.sloth.squared)[,2:3])

colnames(df.slothtest.transformed)[4:5] <- paste0("log.", c("W", "p"))
colnames(df.slothtest.transformed)[6:7] <- paste0("sqrt.", c("W", "p"))
colnames(df.slothtest.transformed)[8:9] <- paste0("squared.", c("W", "p"))

# Table defining the normality at different transformations saved to your drive

write.table(df.slothtest.transformed, file = paste0("AreTransformedSlothsNormal_", round(as.numeric(Sys.time())), ".csv"), row.names = F, sep = ",")


# data subsetted removing "never normal columns" based on normality tests 
df.sloth3 <- subset.data.frame(df.sloth2, select = -c(CONTAG, DIST_SEC_F, PATCH_SHAP, PLAND_SF, SIEI))

# Log transform "PATCH_AREA", "DIST_ROAD","PD", "LPI","AREA_WM"
df.sloth4 <- subset.data.frame(df.sloth2, select = c(PATCH_AREA, DIST_ROAD,PD, LPI,AREA_WM))
df.sloth4.log <- log(df.sloth4)
colnames(df.sloth4.log) <- paste0("log", colnames(df.sloth4.log))

# data subset with only normal non-transformed variables
df.sloth.normal <- subset.data.frame(df.sloth3, select = c(DIST_RIP_F,ECON_AM,EDGE,GYRATE_WM,PATCH_GYRA,PR,SHAPE_WM))

# sqrt CWED
df.sloth3.sqrt <- sqrt(df.sloth3["CWED"])
colnames(df.sloth3.sqrt) <- paste0("sqrt", colnames(df.sloth3.sqrt))

# squared SIDI
sq <- function(x){
  x^2
}

df.sloth3sq <- sq(df.sloth3["SIDI"])
colnames(df.sloth3sq) <- paste0("square", colnames(df.sloth3sq))


# Working Dataset
allsloth <- cbind(df.sloth[,1:6], df.sloth4.log, df.sloth3.sqrt, df.sloth.normal)

str(allsloth)


par(mfrow=c(3,3))
boxplot(logPATCH_AREA ~ OCCURRENCE, data = allsloth, xlab='Absence/Presence', ylab = 'Patch Area')
boxplot(logDIST_ROAD ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Distance to road") 
boxplot(logPD ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Patch Density") #per 100 ha
boxplot(logLPI ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Prop of landscape in patch") #Proportion of landscape in largest patch (%)
boxplot(logAREA_WM ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Mean Area") #Mean area of all patches, weighted by area
boxplot(sqrtCWED ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Density of edge (contrast)") #Density of edge weighted by edge contrast; approaches 0 when contrast is low (m/ha)
boxplot(DIST_RIP_F ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Dist to nearest riparian") #Distance to nearest patch of riparian forest
boxplot(ECON_AM ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "relative prop of edge contrast")
boxplot(EDGE ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Length of edge") #Total length of edge in landscape (m)


par(mfrow=c(3,3))
boxplot(GYRATE_WM ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Mean extent all patches") # Mean extent of all patches, weighted by area
boxplot(PATCH_GYRA ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Exent of patch") #Extent of patch (m)
boxplot(PR ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "# of patch types") #Number of different patch types
boxplot(SHAPE_WM ~ OCCURRENCE, data = allsloth, xlab= "Absence/Presence", ylab = "Mean shape complexity") #Mean shape complexity of all patches, weighted by area


# GLM
glm.sloth <- glm(OCCURRENCE ~ logPATCH_AREA + LAND_USE + logDIST_ROAD + logPD + logLPI + logAREA_WM + sqrtCWED + DIST_RIP_F + ECON_AM + EDGE + GYRATE_WM + PATCH_GYRA + PR + SHAPE_WM, data=allsloth, family=binomial)

summary(glm.sloth)
par(mfrow=c(2,2))
plot(glm.sloth)


# Test Fit of Selected GLM
hoslem.test(allsloth$OCCURRENCE, fitted(glm.sloth))


# Sloth Data is scaled to make sure models run

allslothnotscaled <- allsloth
str(allsloth)

sc <- scale(allsloth[7:19])
allsloth <- cbind(allsloth[1:6], sc)


headers <- colnames(allsloth)[7:19]
cor.m <- abs(cor(allsloth[,headers]))
cor.m.above75 <- cor.m > 0.75
cor.m.above75


# bootstrap an additional 100 samples

set.seed(537)
vec.boots <- sample(1:25, 100, replace = T)
allsloth <- allsloth[vec.boots,]


# create formulae to define models

ls.formulae.aov <- paste("OCCURRENCE ~", c( # prefix for all model formulae
                        "logDIST_ROAD * DIST_RIP_F",
                        "logDIST_ROAD + DIST_RIP_F",
                        "logDIST_ROAD * DIST_RIP_F + logAREA_WM",
                        "logDIST_ROAD * DIST_RIP_F * logAREA_WM",
                        "logDIST_ROAD * DIST_RIP_F + EDGE",
                        "logDIST_ROAD * DIST_RIP_F * GYRATE_WM",
                        "logDIST_ROAD * DIST_RIP_F + GYRATE_WM",
                        "logDIST_ROAD * DIST_RIP_F * PATCH_GYRA",
                        "logDIST_ROAD * DIST_RIP_F + PATCH_GYRA",
                        "logDIST_ROAD * DIST_RIP_F * PR",
                        "logDIST_ROAD * DIST_RIP_F + PR",
                        "logDIST_ROAD * DIST_RIP_F * SHAPE_WM",
                        "logDIST_ROAD * DIST_RIP_F + SHAPE_WM"
                      ))

ls.formulae <- paste(ls.formulae.aov, "+ (1 | LAND_USE)")


# make the glmer

ls.glmer <- lapply(ls.formulae, function(x) {print(x); glmer(x, data = allsloth, family = binomial)})


# extract AIC
ls.aic <- sapply(ls.glmer, function(x) summary(x)$AICtab[1])


# convert to df, pull out the good models

df.modeleval <- data.frame(model = ls.formulae.aov, aic = ls.aic, index = 1:length(ls.formulae.aov), random = "PatchGroup")
df.modeleval$model <- as.character(df.modeleval$model)

df.modeleval$daic <- df.modeleval$aic - min(df.modeleval$aic)
df.modeleval <- df.modeleval[order(df.modeleval$daic),]
head(df.modeleval)


topmodel <- ls.glmer[[10]]
nexttopmodel <- ls.glmer[[6]]

summary(topmodel)
summary(nexttopmodel)

plot(topmodel)
plot(nexttopmodel)

