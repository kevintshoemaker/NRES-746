# load packages
library(spatstat)
library(tidyverse)
library(raster)
library(rgdal)
library(RColorBrewer) 
library(sp)
library(sf)
library(maptools)
#load demo dataset from sppm package
sppm <- finpines
plot(unmark(sppm), cols = '#756bb1', main = 'spatial point pattern', border = 'red', pch=19)
plot(sppm, cols = '#756bb1', main = 'spatial point pattern with marks', border = 'red')
#load covariate data
  
  #rainfall
  rainfall.raster <- raster("monthly_rainfall.tif")
  col.rain <- colorRampPalette(rev(brewer.pal(5, 'Spectral')))
 
  #dem
  dem <- raster("SJER2013_DSM.tif")
  
  # plot example covariates
  par(mfrow = c(1,2))
  plot(rainfall.raster, main = "Average monthly rainfall (mm)", col = col.rain(255), labels = F)
  plot(dem, main = "Elevation", col = terrain.colors(100), labels = F)
  
fp <- read.table("finpines.txt", header=TRUE)
head(fp)
fin_pines <- as.ppp(fp, owin(c(-5,5), c(-8,2)))
fin_pines
data(bei)
summary(bei)

lamb <- summary(bei)$intensity
lamb
data("ants")
summary(ants)
#Create object to hold density function
ants_dens <- density(ants)
summary(ants_dens)


# Density plot with contour lines and points showing where the nests of the two species of ants are located.
plot.new()
plot(ants_dens, main ="Density of Ant Nests")
contour(ants_dens,add=T)
points(ants, pch=20)
#Ripley's K for ants dataset

n <- 100
ants_ripleyK <- envelope(ants,fun=Kest,nsim=n, verbose=FALSE) 

plot(ants_ripleyK)

#Create subsets of data for each ant species
cataglyphis <- subset(ants, marks=="Cataglyphis", drop=T)
messor <- subset(ants, marks=="Messor", drop=T)

##Density plot for Cataglyphis species with points added
plot(density(cataglyphis), main = "Density of C. bicolor ant nests")
points(cataglyphis, pch=20)
##Density plot for Messor species with points added
plot(density(messor), main = "Density of M. wasmanni ant nests")
points(messor, pch=20)
#K function plots for each species of ants
#cataglyphis
plot(envelope(cataglyphis,fun=Kest,nsim=n,verbose=F))

#Messor
plot(envelope(messor,fun=Kest,nsim=n,verbose=F))
#Array of K function plots for ants dataset
ants_all_k_plots <- alltypes(ants,"K",envelope=T,verbose=F)
plot(ants_all_k_plots)
# paired distances
pairdist(finpines)[1:3, 1:5]

#nearest neighbor distances 
nndist(finpines)[1:5]

# empty space distanecs 
Z <- distmap(finpines)
plot(Z)
# clark evans test
clarkevans.test(finpines, correction = "donnelly")
plot(Fest(finpines))
summary(bei)
class(bei)
plot(bei)
plot(quadratcount(bei)) # shows a count of all the data in quadrats
fit <- ppm(bei ~ 1)
fit
class(fit)
exp(-4.933)
bei.extra
fit2 <- ppm(bei ~ grad, data = bei.extra)
fit2

plot(effectfun(fit2, "grad", se.fit = T))
fit3 <- ppm(bei ~ elev + grad, data=bei.extra)
fit3
coefs <- fit3$coef

grad = bei.extra$grad
max(grad)
imp.grad <- exp(coefs[3]*max(grad))
imp.grad

elev = bei.extra$elev
dif.elev <- max(elev) - min(elev)
min.elev <- min(elev)
imp.elev <- exp(coefs[2]*dif.elev)
imp.elev
fit4 <- ppm(bei ~ polynom(grad, elev, 2), data=bei.extra)
lamhat <- predict(fit4)

lamB <- predict(fit4)
plot(lamB)

plot(predict(fit4, se=TRUE)$se)
M <- persp(bei.extra$elev, olin=lamhat, colmap=topo.colors,
  shade=0.4, theta=-55, phi=25, expand=6,
  box=FALSE, apron=TRUE, visible=TRUE, main = "Terrain with Intensity of Fitted Model")
perspPoints(bei, Z=bei.extra$elev, M=M, pch=20, cex=0.1)
#this can be computed by using the predict.ppm() and by setting type = "count" and by specifying B as the window. 
#in this way we can compute the expected number of trees at elevations below 130 meters. 
B <- levelset(bei.extra$elev, 130)
predict(fit4, type="count", window=B)
predict(fit4, B, type="count", se=TRUE)

predict(fit, B, type="count", interval="confidence")


#In this case it is also meaningful to compute a prediction interval for the random number of trees in the specified region:
predict(fit, B, type="count", interval="prediction")

#Locations of 62 seedlings and saplings of California redwood trees.
plot(redwood)
#A point pattern giving the locations of 3605 trees in a tropical rain forest. Accompanied by covariate data giving the elevation (altitude) and slope of elevation in the study region.
plot(bei)
