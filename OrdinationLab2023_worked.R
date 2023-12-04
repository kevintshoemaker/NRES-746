#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$       NRES 746 Lab Script      $#
#$       Intro to Ordination      $#
#$ Martin Genova & Kierstin Acuna $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Load libraries ----
library(codep)
library(vegan)

# Load Doubs fish Data ----
data(Doubs)

# Process data ---

species <- as.data.frame(Doubs.fish[-8,])  # remove observation with no data
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))

nrow(species)  # 29 sites (30 but one removed)
ncol(species)  # 27 species


## Explore the data ----

# Count the number of species frequencies in each abundance class
length(unlist(species))   # 783 observations in the species data (species*sites)
ab <- table(unlist(species))
ab      # as many as 5 individuals of a single species was observed at a single site...

# Plot distribution of species frequencies 
barplot(ab, xlab = "Abundance class", ylab = "Frequency", col = grey(5:0/5))   # not sure why this is needed


sort(colSums(species))  # most abundant species is LOC (73 observations). Least abundant is CHA and OMB (15 individuals)
sort(rowSums(species))  # at the site level, some sites have as low as 3 captures and others as much as 89


# Look at the spread of the predictor variables
# ?Doubs
summary(vars)    # lots of variability in the scales of the predictor variables... 


# Exercise 1 ----

# Build a function to prepare the data for RDA

# what is x and what is y?
RDA_prep <- function(x, y){
  # Step 1: hellinger-transform the y data
  as.data.frame(vegan::decostand(species, method = "hellinger"))
  # Step 2: center and scale the x data
  
}

# Use the function to prepare our data for RDA
data <- RDA_prep(x = ,
                 y =
                   )

# Exercise 2 ----

# Step 1: Create a global RDA using all the predictor variables, and a null RDA 
#   using none of the predictor variables.
rda_global <- vegan::rda(          )

rda_null <- vegan::rda(          )  

# Step 2: Use the ordiR2step function for variable selection. Define the object as the 
#   null rda, and the scope as the global rda.
selection <- vegan::ordiR2step(         )

# Step 3: View the results. The formula in the "Call" contains the variables that 
#   ordiR2step selected.
selection

# Exercise 3 ----
# Step 1: Use the rda function to run an rda on the variables selected in exercise 2
rda_final <- vegan::rda(     )

# Step 2: Visualize the rda using the ordiplot function
vegan::ordiplot(     )

# Step 3: Calculate the adjusted R-squared and p-values for the rda
r2 <- vegan::RsquareAdj(     )

p <- vegan::anova.cca(     )






# Answers ----
## Exercise 1 ----
RDA_prep <- function(x, y){
 
  # Step 1: hellinger-transform the y data
  yhell <- matrix(NA, nrow = nrow(y), ncol = ncol(y)) #create a storage matrix
  y <- as.matrix(y) #convert y to a matrix
  for (i in 1:nrow(y)) {
    yhell[i,] <- sqrt(y[i,]/sum(y[i,])) #use the hellinger-transformation function
  }
  colnames(yhell) <- colnames(y) #assign the species names as column names to the new
                                 #  hellinger-transformed matrix

  # Step 2: center and scale the x data
  xscale <- scale(x)
  
  # Output
  out <- list()
  out[["x"]] <- as.data.frame(xscale)
  out[["y"]] <- as.data.frame(yhell) 
  return(out)
}

# Compare our answer to the vegan::decostand answer
our_hell <- RDA_prep(vars, species)$y
vegan_hell <- vegan::decostand(species, "hellinger")
which(our_hell != vegan_hell) #they're the same!

# Use our function to prepare our data for RDA
data <- RDA_prep(vars, species)

## Exercise 2 ----
# Step 1: Use the RDA_prep function to prepare the data. Then, create a global RDA 
#   using all the predictor variables, and a null RDA using none of the predictor variables.

rda_global <- vegan::rda(data$y ~ ., data=data$x)

rda_null <- vegan::rda(data$y ~ 1, data=data$x)

# Step 2: Use the ordiR2step function for variable selection
selection <- vegan::ordiR2step(object = rda_null, scope = rda_global)

# Step 3: View the results. The formula in the "Call" contains the variables that 
#   ordiR2step selected.
selection

## Exercise 3 ----

# Step 1: Use the rda function to run an rda on the variables selected in exercise 2
rda_final <- vegan::rda(data$y ~ DFS + oxy + bdo + slo + Lat, data=data$x)
summary(rda_final)

# Step 2: Visualize the rda using the ordiplot function
vegan::ordiplot(rda_final, type = "text")

# Step 3: Calculate the adjusted R-squared and p-values for the rda
r2 <- vegan::RsquareAdj(rda_final)

p <- vegan::anova.cca(rda_final, by="term")
















