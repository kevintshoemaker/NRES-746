
############################################################
####                                                    ####  
####  NRES 746, Student Lecture                         ####
####                                                    ####
####  Heather Reich, Josh Vasquez, Megan Osterhout      #### 
####  University of Nevada, Reno                        ####
####                                                    #### 
############################################################


############################################################
####  Spatial capture-recapture                         ####
############################################################



library(secr)


capture.history<-read.table("hareCH6capt.txt",col.names = c("Session","ID","Occasion","Detector"))
head(capture.history)
trap.data<-read.table("hareCH6trap.txt",col.names = c("Detector","x","y"))
head(trap.data) ##note that the x,y information is given in meter distances, do NOT use latitude/longitude

hare.secr<-read.capthist("hareCH6capt.txt","hareCH6trap.txt",detector = "single")

summary(hare.secr)

plot(hare.secr, tracks = TRUE) 

movements <- unlist(moves(hare.secr))

hist(movements, breaks = seq(-61/2, 500,61), xlab = "Movement", main = "") #another illustration that majority of movements are within 100m


initialsigma <- RPSV(hare.secr, CC = TRUE)
cat("Quick and biased estimate of sigma =", initialsigma, "m\n") #our spatial scale parameter that determines how rapidly capture probability declines with distance AKA our state model 


fit <- secr.fit (hare.secr, buffer = 4 * initialsigma, trace = FALSE) 

esa.plot(fit)
abline(v = 4 * initialsigma, lty = 2, col = 'red')    
#The theory of SECR tells us that buffer width is not critical as long as it is wide enough that animals with activity centers close to the edge of the buffer have effectively zero chance of appearing in our sample. 

#We can check this with the function esa.plot and see that the estimated density has easily reached a plateau at the chosen buffer width (dashed red line)


print(fit)

#The printed report gives us the following useful information:
#      function call and time stamp
#      summary of the data
#      model description (with log likelihood and AIC)
#      estimates of the coefficients
#      estimates of covariance matrix of coefficients
#      estimates of the real parameters


plot(fit, limits = TRUE) 




fit.halfnormal <- secr.fit (hare.secr, buffer = 4 * initialsigma, detectfn = 'HN', trace = FALSE)
fit.negative.exponential <- secr.fit (hare.secr, buffer = 4 * initialsigma, detectfn = 'EX', trace = FALSE)
fits <- secrlist(HN = fit.halfnormal, EX = fit.negative.exponential)#compare the models
predict(fits)
AIC(fits)

library(gstat)
lizard<-hornedlizardCH
head(lizard)
par(mar=c(1,1,2,1))
plot(hornedlizardCH, tracks = TRUE, varycol = FALSE, lab1cap = TRUE, laboffset = 8,
     border = 10, title ='')


FTHL.fit <- secr.fit(hornedlizardCH, buffer = 80, trace = FALSE)
predict(FTHL.fit)
plot(FTHL.fit, xv = 0:80, xlab = "Distance (m)", ylab = 'p', main="Probability of detection")

datadir <- system.file("extdata", package = "secr")
examplefile1 <- paste0(datadir, '/polygonexample1.txt')
polyexample1 <- read.traps(file = examplefile1, detector = 'polygon')
polygonCH <- sim.capthist(polyexample1, popn = list(D = 1, buffer = 200),
                          detectfn = 'HHN', detectpar = list(lambda0 = 5, sigma = 50),
                          noccasions = 1, seed = 123)
par(mar = c(1,2,3,2))
plot(polygonCH, tracks = TRUE, varycol = FALSE, lab1cap = TRUE, laboffset = 15,
     title = paste("Simulated 'polygon' data", "D = 1, lambda0 = 5, sigma = 50"))

#cuesim.fit <- secr.fit(polygonCH, buffer = 200, trace = FALSE) #Takes a long time to run! 
#predict(cuesim.fit)

discreteCH <- discretize (polygonCH, spacing = 20)
par(mar = c(1,2,3,2))
plot(discreteCH, varycol = FALSE, tracks = TRUE)

discrete.fit <- secr.fit(discreteCH, buffer = 200, detectfn = 'HHN', trace = FALSE)
predict(discrete.fit)

examplefile2 <- paste0(datadir, '/polygonexample2.txt') #this txt file is within secr package
polyexample2 <- read.traps(file = examplefile2, detector = 'polygon')
par(mfrow = c(1,2), mar = c(2,1,1,1))
plot(polyexample2)
text(polyexample2$x, polyexample2$y, 1:nrow(polyexample2), cex = 0.7)
newpoly <- make.poly (list(p1 = polyexample2[11:23,],
                           p2 = polyexample2[c(1:11, 23:27),]))

plot(newpoly, label = TRUE)
