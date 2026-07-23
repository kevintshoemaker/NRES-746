
############################################################
####                                                    ####  
####  NRES 746, Student-led topic #3                    ####
####                                                    ####
############################################################


############################################################
####  Spatial autocorrelation                           ####
############################################################



#first things first: here are the packages you'll need to follow along:

list = c('gstat', 'raster', 'geosphere', 'ape', 'foreach', 'doParallel', 'rgdal')
#install.packages(list)


library(gstat)
library(raster)

xy <- expand.grid(1:100, 1:100) # create a coordinate grid to represent a real landscape.
                                # the larger the grid, the longer it will take to generate data
names(xy) <- c('x','y') # name the axes of the grid


  # Specify a model to create a spatially-autocorrelated z variable.
  # range parameter controls degree of spatial autocorrelation.
saLandscapeModel <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=c(100), 
                 model=vgm(psill=100, range=100, model='Exp'), nmax=20)

saLandscape <- predict(saLandscapeModel, newdata = xy, nsim = 1)
image(saLandscape, axes = FALSE, col = terrain.colors(10))


sample <- xy[sample(nrow(xy), 30),] # Randomly sample pixels from the grid
samplePoints <- SpatialPoints(sample) # Create a SpatialPoints object
colors <- topo.colors(20)


# Extract explanatory values.
names(saLandscape) <- c('x','y','z') # Step 1: rename columns in envImage df for clarity
gridded(saLandscape)=~x+y # Step 2: Data frame must be converted to gridded object to create raster layer
saRaster <- raster(saLandscape) # Step 3: convert gridded object to raster
elev <- raster::extract(x=saRaster, y=samplePoints) # Step 4: extract values.
sample <- cbind(sample, elev)
row.names(sample) <- (1:30) 


print(head(sample))


image(saLandscape, axes = FALSE, col = terrain.colors(10))
points(sample)


semivariogram <- function(value,x,y){
  
  # building empty and null vectors
  dist <- vector(mode="numeric", length=length(value))
  semivar <- vector(mode="numeric", length=length(value))
  distance <- c()
  semivariance <- c()
  
  # calculating all possible
  
  for (i in 1:length(value)) {                    # these loops compare all the values with each other
    for( j in 1:length(value)) {
      dist[j] <-  sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)           # measuring the distance
      semivar[j] <- (value[i]-value[j])^2                       # calcualting the semivariance
    }
    distance <- c(distance,dist)
    semivariance <- c(semivariance,semivar)
  }
  
  plot(distance, semivariance, xlab="Distance Between Points", ylab="Squared Difference")
}


semivariogram(sample$elev,sample$x,sample$y)


semivariogram.ll <- function(value,lat,lon){
  library(geosphere)
  
  # building empty and null vectors
  dist <- vector(mode="numeric", length=length(value))
  semivar <- vector(mode="numeric", length=length(value))
  distance <- c()
  semivariance <- c()
  
  # calculating all possible
  
  for (i in 1:length(value)) {
    for( j in 1:length(value)) {
      dist[j] <- distm(c(lon[i], lat[i]), c(lon[j], lat[j]), fun = distHaversine) # measuring the euclidean distance [meters]
      semivar[j] <- (value[i]-value[j])^2          # calcualting the semivariance
    }
    distance <- c(distance,dist)
    semivariance <- c(semivariance,semivar)
  }
  
  plot(distance, semivariance, xlab="Distance [meters]", ylab="Squaired Difference")
  #data <- cbind(distance,semivariance)
  #return(data)
}

week1_transects <- read.csv("SnowEx17_GPR_Week1_transects.txt")
tail(week1_transects)

########### seperating the data by transect ###########
transect46 = week1_transects[which(week1_transects$TRANSECT == 46),]
transect1 = week1_transects[which(week1_transects$TRANSECT == 1),]

# building a small data set for the semivarance
df.new <- transect46[seq(1, nrow(transect1), 200),]
#df.new <- week1_transects[seq(1, nrow(week1_transects), 1000),] #this will take all the transects

# seprating data for the function
lat <- df.new$LAT
lon <- df.new$LONG
x <- df.new$SNOW.DEPTH..m..assuming.velocity.of.0.234.m.ns.

# run the function
#semivariogram.ll(x,lat,lon)


semivariogram.ll(x,lat,lon)


# building a small data set for the semivarance
df.new <- transect46[seq(1, nrow(transect1), 500),]
#df.new <- week1_transects[seq(1, nrow(week1_transects), 1000),] #this will take all the transects

# seprating data for the function
lat <- df.new$LAT
lon <- df.new$LONG
x <- df.new$SNOW.DEPTH..m..assuming.velocity.of.0.234.m.ns.

#install.packages("ape")
library("ape")
#look at package "ape"
#??ape
#function we will use is Moran.I()
#?Moran.I
#look at the help and see the Moran's forumla
#requires "x" - a numeric vector - which is our variable.
#requires "weight" - a matrix of weights - calculated using dist().


#generating IDW matrix
#use dist() to compute and return the distance matrix between rows of data
#?dist
sample.distances<-as.matrix(dist(cbind(sample$x,sample$y))) #distance weight matrix
sample.distances.inverse<-1/sample.distances #inverse distance weight matrix
sample.distances.inverse[1:5,1:5] #however "infin" problem"


#Remove infinity and replace with 0's - occurs because 1 divide by 0 is infinity
diag(sample.distances.inverse)<-0 
sample.distances.inverse[1:5,1:5]


#Morans I
#?Moran.I
#Weights are obtained, now place desired variable to test as "x"
Moran.I(x=sample$elev,weight=sample.distances.inverse)
#Results!


## Parallel Spatial Autocorrelation Simulation
library(parallel)
library(foreach)
library(doParallel)

numCores <- detectCores()-1
registerDoParallel(numCores)

resultsList<-foreach(i=1:100, .combine = c) %do% {     # %dopar%
  library(gstat)
  library(raster)
  library(sp)
  
  xy <- expand.grid(1:100, 1:100) # create a coordinate grid to represent a real landscape.
  # the larger the grid, the longer it will take to generate data
  names(xy) <- c('x','y') # name the axes
  
  # If desired, model deterministic component of explanatory variable on landscape surface
  # by altering the beta parameters in the gstat() call.
  # Currently set to 0 (no deterministic relationship between space and explanatory variable)    
  envDeterm <- gstat::gstat(formula=z~1+x+y, locations=~x+y, dummy=T, beta=c(1,0,0), 
                     model=gstat::vgm(psill=0, range=1, model='Exp'), nmax=1)
  
  # Model spatially-autocorrelated component of explanatory variable on landscape surface.
  # range parameter controls degree of spatial autocorrelation. Set to 1 to eliminate spatial autocorrelation.
  envSA <- gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=c(0), 
                 model=gstat::vgm(psill=.025, range=20, model='Gau'), nmax=2)
  
  # Model stochastic (error) component of explanatory variable on landscape surface.
  envErr <- gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=c(0), 
                  model=gstat::vgm(psill=.025, range=1, model='Gau'), nmax=1)
  
  # Simulate data from each model across the landscape surface.
  # Can visualize each surface using image() if desired
  eDe <- predict(envDeterm, newdata=xy, nsim=1) 
  eSa <- predict(envSA, newdata=xy, nsim=1) 
  eEr <- predict(envErr, newdata=xy, nsim=1)
  # image(eDe)
  # image(eSa)
  # image(eEr)
  
  # Add each component into a single landscape surface.
  envImage <- eDe+eSa+eEr
  # image(envImage)
  
  # Model spatially-autocorrelated component of response variable on landscape surface.
  # range parameter controls degree of spatial autocorrelation. Set to 1 to eliminate spatial autocorrelation.
  resSA <- gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=c(0), 
                 model=gstat::vgm(psill=.025, range=20, model='Gau'), nmax=2)
  
  # Model stochastic (error) component of response variable on landscape surface.
  resErr <- gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=c(0), 
                  model=gstat::vgm(psill=.025, range=1, model='Gau'), nmax=1)
  
  # Simulate data from response variable SA and error models across the landscape surface.
  # Again, can visualize each surface using image() if desired
  rSa <- predict(resSA, newdata=xy, nsim=1)
  rEr <- predict(resErr, newdata=xy, nsim=1)
  # image(rSA)
  # image(rEr)
  
  # Combine the explanatory variable surface with the response SA and error
  # surfaces to generate predicted response values. If there is no true
  # relationship between the explanatory and response variables, the beta1
  # parameter is set to zero.
  beta1 <- 0
  resImage <- (beta1*envImage$sim1) + rSa + rEr 
  
  # Create random sample points
  sample <- xy[sample(nrow(xy), 100),] # Randomly sample pixels from the grid
  samplePoints <- sp::SpatialPoints(sample) # Create a SpatialPoints object
  
  # Extract explanatory and response values using the sampling locations. There
  # is probably a simpler way to do this, but in this case I am creating rasters
  # from the envImage and resImage data frames, then using raster::extract to
  # extract the explanatory and response values and combine them in a new data
  # frame of simulated data.
  
  # Extract explanatory values.
  names(envImage) <- c('x','y','explanatory') # Step 1: rename columns in envImage df for clarity
  sp::coordinates(envImage) <- ~x+y  # Step 2: Data frame must be converted to gridded object to create raster layer
  sp::gridded(envImage) <- TRUE
  envRaster <- raster::raster(envImage) # Step 3: convert gridded object to raster
  explanatory <- raster::extract(x=envRaster, y=samplePoints) # Step 4: extract values.
  
  # Extract response values. Repeat steps 1-4.
  names(resImage) <- c('x','y','response')
  sp::coordinates(resImage) <- ~x+y
  sp::gridded(resImage) <- TRUE
  resRaster <- raster::raster(resImage)
  response <- raster::extract(x=resRaster, y=samplePoints) 
  
  # Bind locations, explanatory values, and response values into a singe data
  # frame of simulated data
  sample <- cbind(sample,explanatory,response) 
  
  # Create simple linear model and extract the p-value.
  model <- stats::lm(response ~ explanatory)
  pvals <- summary.lm(model)$coefficients[,4]
  # pvals
  i <- as.numeric(pvals[2]) # Store the p-value for the beta1 parameter in the resultsList
}

# Count number of positive results (p value > 0.05)
for (i in 1:length(resultsList)){
  if (resultsList[i] >= 0.05){
    resultsList[i] <- FALSE
  }else(resultsList[i] <- TRUE)
}

typeIerror <- sum(resultsList)/100
pvalue <- 0.05
simulationResults <- c(pvalue, typeIerror)
names(simulationResults) <- c('P-value', 'Type-I error rate')
print(simulationResults)


library(rgdal)
library(raster)
age.data1 <- read.csv('age.data1.csv')
print(head(age.data1))

cropped <- raster('AutoCovRaster.tif')
plot(cropped, axes = FALSE, col = terrain.colors(10))
points(age.data1[,1], age.data1[,2])



semivariogram1 <- function(value,x,y){                     
  dist <- vector(mode="numeric", length=length(value))    
  semivar <- vector(mode="numeric", length=length(value))
  distance <- c()
  semivariance <- c()
  for (i in 1:length(value)) {
    for( j in 1:length(value)) {
      dist[j] <-  sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)           
      semivar[j] <- (value[i]-value[j])^2                       
    }
    distance <- c(distance,dist)
    semivariance <- c(semivariance,semivar)
  }
  plot(distance, semivariance, xlab="Distance Between Points", ylab="Squaired Difference between residuals", ylim=c(0,250))
}


basic.model <- lm(heights~age, data=age.data1)
basic.model$coefficients


semivariogram1(basic.model$residuals,age.data1$x,age.data1$y)


linear.cor <- function(dist.matrix=distance){  # helper function for autocovariate.model()
  weight.part1 <- (dist.matrix*(-1/40) + 1)    # defines how strongly neighboring points are correlated with each other as a function of distance between them
}                                              # ie. determines what the values of the weight matrix will be



# Autocovariate.model() is a function which appends an autocovariate value to each datapoint.
# The autocovariate values are then used as a second predictor variable in you linear model
# The inputs are your data and your correlation function (uses the above function, linear.cor() by default)

autocovariate.model <- function(data=data1, cor.func=linear.cor){
  dist.temp <- dist(data[,1:2])      # create a distance matrix from the data
  distance <- as.matrix(x=dist.temp) # create a distance matrix from the data
  weight.part1 <- cor.func(distance) # next few lines create a weighted distance matrix based on your correlation function
  no.negs <- function(x){
    if(x>=0){x<-x}else{x<-0}
    return(x)
  }
  weight <- apply(weight.part1, c(1,2), no.negs)   # neibors beyond a given threshold should have no influance rather than a negative correlation, so any negative values are converted to 0
  weighted.heights <- matrix(nrow=nrow(data), ncol=nrow(data)) # set up matrix which will be the weight matrix times the response variable
  for(i in 1:nrow(data)){
    weighted.heights[i,] <- (weight[i,]*data$heights)
  }
  autocov <- numeric(nrow(data))                   # calculates the autocovariate value at each site
  for(i in 1:nrow(data)){                          # autocov is the sum a given site's neighbors' (weight*response variable)
    autocov[i] <- sum(weighted.heights[i,])-weighted.heights[i,i]
  }
  data.autocov <- cbind(data,autocov)
  return(data.autocov)
}


age.data2 <- autocovariate.model(age.data1)            
head(age.data2)

basic.model <- lm(heights~age, data=age.data2)
corrected.model <- lm(heights~age+autocov, data=age.data2)


semivariogram1(basic.model$residuals,age.data2$x,age.data2$y)
title(main='basic model')
semivariogram1(corrected.model$residuals,age.data2$x,age.data2$y) # clearly still autocorrelation in the residuals, but it is significantly reduced
title(main='corrected model')
 
require ("ape")                # see how much correlation in the residuals is reduced acording to moran's I with a generic 1/distance weight matrix
sample.distances<-as.matrix(dist(cbind(age.data2$x,age.data2$y)))   # first calculate inverse distance weights (IDWs)
sample.distances.inverse<-1/sample.distances
diag(sample.distances.inverse)<-0             # get rid of infinity for 0's because 1 divide by 0 is infinity

Moran.I(basic.model$residuals,sample.distances.inverse, scaled=TRUE)
Moran.I(corrected.model$residuals,sample.distances.inverse, scaled=TRUE) # Moran's I also shows significant reduction


summary(basic.model)$r.squared
summary(corrected.model)$r.squared


basic.model$coefficients                  
corrected.model$coefficients
