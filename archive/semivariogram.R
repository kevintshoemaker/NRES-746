# Identifying Spatial Autocorrelation with a Semivariogram ----

## Load Prevalence Data ----
# load gambiadf.csv, or data of your choice
filename <- "StudentLedTopics/RegModelling_SpatiallyCorrelatedRandomEffects/gambiadf.csv" ## change to reflect the true file location
df <- read.csv(file = filename, header = TRUE)

head(df)
names(df)

responsevar <- df$prev # change this to fit your data
ResponseVarString <- "Malaria Prevalence" # change this to fit your data
xvar <- df$x # column that represents x coordinates, in UTM units (meters)
yvar <- df$y # column that represents y coordinates, in UTM units (meters)

hist(responsevar, main=ResponseVarString)

## Calculate squared differences and distances between all possible pairs of observations ----
# initialize vectors for squared difference between all pairs:
# (prev2-prev1)^2
obs_sqdiff <- 0
# and distance between villages:
# sqrt((x2-x1)^2+(y2-y1)^2)
pair_dist_km <- 0

# initialize counter variable for "for" loop
counter <- 0
# first for loop
for (i in 1:(length(responsevar)-1)){
  # second for loop
  for (j in (i+1):length(responsevar)){
    counter <- counter + 1 # increment counter
    # difference in responsevar, squared
    obs_sqdiff[counter] <- (responsevar[i]-responsevar[j])^2
    # 
    pair_dist_km[counter] <- sqrt((xvar[i] - xvar[j])^2 + (yvar[i] - yvar[j])^2)/1000
  } # end first for loop
} # end second for loop
plot(pair_dist_km, obs_sqdiff, xlab="Distance between obs(km)", ylab="Squared diff in prevalence")

## Bin the Distances ----
bins <- cut(pair_dist_km, breaks=c(seq(0,280,by=10)), right = F)
  # I chose 10km bins, what happens if you vary this?

## Calculate Semivariance ----
semivar <- 0.5*tapply(obs_sqdiff, bins, FUN = mean)
avgdist <- tapply(pair_dist_km, bins, FUN = mean)

## Plot Semivariance ----
plot(avgdist, semivar, ylim = c(0,0.11), main="Semivariogram", xlab="Average Distance (km)", ylab="prevalence semivariance")

# Check with geostats::semivariogram() function ----
library(geostats)
sv1 <- geostats::semivariogram(
  x=xvar,
  y=yvar,
  z=responsevar,
  bw = 10000,
  nb = 28,
  plot = TRUE,
  fit = TRUE,
  model = c("spherical", "exponential", "gaussian")
)

# Semivariogram of Residuals ----
# Model Prevalence ~ Altitude
mod <- lm(df$prev ~ df$alt)
plot(mod)

sv2 <- geostats::semivariogram(
  x=df$x,
  y=df$y,
  z=mod$residuals,
  bw = 10000,
  nb = 28,
  plot = TRUE,
  fit = TRUE,
  model = c("spherical", "exponential", "gaussian")
)