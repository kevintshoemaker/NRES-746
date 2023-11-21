
# load packages --------------------

library(ResourceSelection)
library(lubridate)
library(tidyr)
library(terra)
library(spdep)
library(sf)
library(dplyr)
library(sp)
library(adehabitatHR)
library(scales)
library(raster)
#library(rgdal) # Removed from CRAN in October! Replaced by terra.
library(ggplot2)
library(amt, quietly = T)
library(raster, quietly = T)
library(glmmTMB)


# read in data ----------------

sheeps <- read.csv("rsf_intro.csv")

head(sheeps) # What does the first few rows look like?


# Visualize data -------------

ggplot(sheeps, aes(x = JulianDay)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, max(table(sheeps$JulianDay)), by = 100000)) +
  labs(title = "Spring '20 Data Distribution", x = "April 1-June 30", y = "Number of Records") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

colnames(sheeps)[1] <- 'ID' # Bug fix, don't mind
n_indivs <- sheeps %>% 
  distinct(ID) %>% 
  n_distinct() 
n_indivs ## 20 individuals 
ggplot(sheeps, aes(x=Sex, fill=Sex))+
  geom_bar() + 
  labs(title= "Number of Records for Males vs. Females", 
       x="Sex", 
       y="Number of Records") +
  theme_minimal()


mf_counts <- table(sheeps$Sex)
mf_counts
count_table <- table(sheeps$ID)
barplot(count_table, main = "Number of Records for Each Individual",
        xlab = "Individual", ylab = "Number of Records",
        col = "skyblue", border = "black")
plot(Northing~Easting, sheeps, col=ifelse(Sex=="M","blue","red"), pch=16)



# create study region polygon and make random points ----------------

# (In this example, we will create an MCP using a helper package called SpatialPoints)
# SpatialPoints data frame objects don't like missing values
# Let's remove the first two rows that contain NA values and...
# Only include the x and y coordinates for making our MCP:
sheeps.sp <- sheeps[, c("Easting", "Northing")] 
spatial <- SpatialPoints(sheeps.sp)

# Calculate MCPs for sheep:
sheep.mcp <- mcp(spatial, percent = 100)

# Let's examine the output!
sheep.mcp # Calling the MCP object will give us an estimate of area
plot(Northing~Easting, sheeps, col=ifelse(Sex=="M","blue","red"), pch=16)
plot(sheep.mcp, col = alpha(1:5, 0.5), add = TRUE)

sheep.mcp <- st_as_sf(sheep.mcp) # Let's use our MCP that we created and convert it to a spatial data frame
class(sheep.mcp) # You can use the "class" function to verify that your "sheep.mcp" object is indeed a data frame!

     # generate random points
rand <- st_sample(sheep.mcp, size = nrow(sheeps)*5) # Let's generate a number of random data points equal to the length of our data set (1:1 ratio).

rand_final <- st_coordinates(rand) # Let us also make sure that we save their coordinates 
rand_final <- as.data.frame(rand_final) # and convert this to a data frame




plot(Northing~Easting, sheeps, col=ifelse(Sex=="M","blue","red"), pch=16)
plot(sheep.mcp, col = alpha(1:5, 0.5), add = TRUE)
plot(rand, add=T, pch=16, col="green", cex=0.5) # Let's plot the points!
colnames(rand_final)[1] <- "Easting" 
colnames(rand_final)[2] <- "Northing"
rand_final$Used <- 0 # Random Points
sheeps$Used <- 1 # Used Points
sheeps$Num <- 1:1115
rand_final$Num <- rep(1115:2229,each=5)
rand_merge <- merge(sheeps,rand_final, by = "Num", all=TRUE)
rand_merge[is.na(rand_merge)] <- ""
rand_merge$Easting <- paste(rand_merge$Easting.x, rand_merge$Easting.y)
rand_merge$Northing <- paste(rand_merge$Northing.x, rand_merge$Northing.y)
rand_merge$Used <- paste(rand_merge$Used.x, rand_merge$Used.y)
rand_merge <- subset(rand_merge, select = -c(Easting.x, Easting.y, Northing.x, Northing.y, Used.x, Used.y))
copydata <- c("ID", "Sex", "JulianDay")

## STOP HERE AND JUST START WITH EXERCISE
rand_merge[1116:2229, copydata] <- rand_merge[1:1115, copydata]
all_sheep <- rand_merge
all_sheep$Easting <- as.numeric(all_sheep$Easting)
all_sheep$Northing <- as.numeric(all_sheep$Northing)
all_sheep <- na.omit(all_sheep)
missing <- apply(all_sheep, 2, function(x) any(is.na(x)))
all_sheep <- all_sheep[complete.cases(all_sheep), ]
plot(Northing~Easting, rand_merge, col=ifelse(Sex=="M","blue","red"), pch=16, cex=0.5)
raster <- "samplestudy.tif"
pfg <- rast(raster)
print(pfg) # This provides us some basic metadata for our beautiful raster file
plot(pfg) # What does it look like?
# Let's convert the sheep and raster image to use same datum:  
# epsg <- 26913
crs.dat <- CRS("+init=epsg:26913")
sheep_spat <- SpatialPointsDataFrame(coords = all_sheep[,c("Easting","Northing")], 
                                     data=all_sheep, 
                                     proj4string = crs.dat)

#plot(sheep_spat, pch=16, col="blue") # Where are our data points in relation to the raster image?

crs.pfg <- crs(pfg)
pfg_proj <- terra::project(pfg, crs.dat) # Projection using the "terra" package

plot(pfg_proj)
plot(geometry(sheep_spat), pch=16, cex=0.5, col="blue", add=T) # Plot both the raster and our data together!
# Extract point values
point_values <- terra::extract(pfg_proj, sheep_spat@coords)

# Add point values back to main data frame ##

all_sheep <- cbind(all_sheep,point_values)
colnames(all_sheep)[8] <- "pfg"
# hist(all_sheep$samplestudy) # Distribution of PFG values for you to explore!

# View final data frame #
head(all_sheep)

all_sheep$Used <- as.numeric(all_sheep$Used)
all_sheep$JulianDay <- as.numeric(all_sheep$JulianDay)

# START EXERCISE ------------------

rm(list=ls())

df <- read.csv("rsf_lab.csv")



mean.pfg <- mean(all_sheep$pfg)
sd.pfg <-  sd(all_sheep$pfg)
all_sheep$pfg.s <- (all_sheep$pfg-mean.pfg)/sd.pfg

mean.jul <- mean(all_sheep$JulianDay)
sd.jul <-  sd(all_sheep$JulianDay)
all_sheep$Julian.s <- (all_sheep$JulianDay-mean.jul)/sd.jul

all_sheep$ID2 <- as.numeric(as.factor(all_sheep$ID))

## Build your GLMM:: 

spring_mod0 <- glmmTMB(Used ~ Sex + JulianDay + pfg + (1|ID), 
                      data = all_sheep, family=binomial(link="logit"),
                      na.action = "na.fail", REML=FALSE)
# summary(spring_mod0)


# In this example, we introduce a random effect between individual and perennial forbs and grasses. 
# Additionally, this model will allow us to determine the difference between females and males (isfem[o]). 

library(jagsUI)
cat("
    model{
    
        ## Specify likelihood:: 
        for(o in 1:n){
          logit(exp_use[o]) <- b0.l+b1*isfem[o]+b2*JulianDay[o]+b3[ID[o]]*pfg[o]+reff[ID[o]]
          used[o] ~ dbern(exp_use[o])
        }
        for(i in 1:nn){
          reff[i] ~ dnorm(0,reff.prec)
          b3[i] ~ dnorm(b3.0,b3.prec) 
        }
        ## Specify Priors::
        
        b0 ~ dunif(0,1)
        b0.l <- log(b0/(1-b0))
        b1 ~ dnorm(0,0.01)
        b2 ~ dnorm(0,0.01)
        b3.0 ~ dnorm(0,0.01)
        b3.prec <- pow(b3.sd,-2)
        b3.sd ~ dunif(0,10)
        reff.sd ~ dunif(0,10)
        reff.prec <- pow(reff.sd, -2)
    }
    
    ", file="jagscode.txt" 
    
)
## define isfem, JulianDay, pfg, ID, n, nn

## Package data for JAGS

jagslist <- list(
  isfem = ifelse(all_sheep$Sex=="F", 1, 0), 
  JulianDay = all_sheep$Julian.s,
  pfg = all_sheep$pfg.s, 
  ID = all_sheep$ID2, 
  n = nrow(all_sheep), 
  nn = max(all_sheep$ID2),
  used = all_sheep$Used
  
)
#jagslist

### Define inits for Jags:: 

initials <- function(){
  list(
    b0 = runif(1,0.4,0.6), 
    b1 = rnorm(1,0,0.1),
    b2 = rnorm(1,0,0.1),
    b3.0 = rnorm(1,0,0.1),
    b3.sd = runif(1,0.1,0.5),
    reff.sd = runif(1,0.1,0.5)
  )
  
}
#initials()

params <- c("b0","b1","b2","b3", "b3.0", "b3.prec","reff.sd","reff")
#?jags

model <- jags(data=jagslist, inits=initials, parameters.to.save=params, model.file="jagscode.txt",
              n.chains=3, n.adapt=1000, n.iter=10000, n.burnin=5000, n.thin=2,
              parallel=TRUE)

sims = model$sims.list

# Some histograms for you to look at!
# hist(sims$b0)
# hist(sims$b1)
# hist(sims$b2)
# hist(sims$b3[,5]) # effect for individual 1
# hist(sims$b3.0)
# hist(sims$reff.sd)
# hist(sims$reff[,1])
