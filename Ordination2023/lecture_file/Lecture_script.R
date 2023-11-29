#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$     NRES 746 Lecture Script    $#
#$       Intro to Ordination      $#
#$ Martin Genova & Kierstin Acuna $#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


# Load libraries ----
library(codep)
library(vegan)
library(datasets)
library(tidyverse)


# Load Doubs fish Data ----

data(Doubs)
?Doubs #Information about the Doubs Fish data set.
species <- as.data.frame(Doubs.fish[-8,])
vars <- as.data.frame(cbind(Doubs.env[-8,],Doubs.geo[-8,]))


# GLM example ----
  # Cottus gobio gamma distributed model

CHA.alt <- species$CHA + 1 #get rid of zeros to prep for log link
COGO_mod <- glm(CHA.alt ~ pH + flo + oxy, data = vars, family = Gamma(link = "log"))
summary(COGO_mod) #oxygen is significantly positive!


# Dissimilarity Measures ----


## Types of Distances Coefficients ----


### Euclidean Distance ----

(Y.hmm <- data.frame(hydrophillic_1 = c(1, 0, 0), hydrophillic_2 = c(1, 1, 0),
                    mesic_1 = c(0, 1, 0), mesic_2 = c(0,4,0),
                    xeric_1 = c(0, 1, 3),xeric_2 = c(0, 0, 2), 
                    row.names = c("sample_1_wet", "sample_2_intermediate",
                                                        "sample_3_dry")))

# Calculate Euclidean distance using the dist() function
(Y.hmm.DistEu <- as.matrix(dist(x = Y.hmm, method = "euclidean")))

# Calculate Euclidean Distance by Hand
calc_eu_dist <- function(spe_abun_df) {
  # Create output matrix
  output <- as.data.frame(matrix(NA, nrow = nrow(spe_abun_df), ncol = nrow(spe_abun_df)))
  # Index through the rows of the data frame
  for (i in 1:nrow(spe_abun_df)) {
    x1 <- spe_abun_df[i, ]
    for (t in 1:nrow(spe_abun_df)) {
      x2 <- spe_abun_df[t,]
      # Calculate euclidean distance and place distance into output data frame
      output[i,t] <- sqrt(sum((x1 - x2)^2))
    }
  }
  # Return output
  return(output)
}

# Run Euclidean distance by hand function
(Y.hmm_eu_dist <- calc_eu_dist(Y.hmm))


### Bray-Curtis Coefficient ----

(Y.hmm.BCdist <- vegan::vegdist(Y.hmm, method = "bray", binary = FALSE))

(Y.hmm.BCdist.matrix <- as.matrix(Y.hmm.BCdist))

# Calculate Bray-Curtis coefficient by hand
calc_bc_dist <- function(spe_abun_df) {
  # Create output matrix
  output <- as.data.frame(matrix(NA, nrow = nrow(spe_abun_df), ncol = nrow(spe_abun_df)),
                          row.names = rownames(spe_abun_df))
  colnames(output) <- rownames(spe_abun_df)
  # Index through the rows of the data frame
  for (i in 1:nrow(spe_abun_df)) {
    x1 <- spe_abun_df[i, ]
    for (t in 1:nrow(spe_abun_df)) {
      x2 <- spe_abun_df[t,]
      # Create empty data frame to find the minimum values of each species between two sites
      comp_df <- as.data.frame(matrix(nrow = 2, ncol = ncol(spe_abun_df)))
      # Place the site values into the data frame
      comp_df[1,] = x1
      comp_df[2,] = x2
      # Find the minimum abundance values of each species and sum them.
      min_abundances <- apply(comp_df, 2, min)
      W <- sum(min_abundances)
      # Sum the abundances of site 1
      A = sum(x1)
      # Sum the abundances of site 2
      B = sum(x2)
      # Calculate the Bray-Curtis coefficient
      bc_dist <- (1 - ((2 * W) / (A + B)))
      # Place the BC coefficient into the output data frame
      output[i,t] <- bc_dist
    }
  }
  # Return output
  return(output)
}

# Run Bray-Curtis coefficient by hand function
calc_bc_dist(Y.hmm)


# Unconstrained Ordination ----


## Principal Component Analysis ----


### Breaking Down a PCA ----

# 1. Load the data set

data(Doubs)
species <- Doubs.fish[-8,]

# 2. Apply a Hellinger transformation on the Species abundance data to standardize the data

spe.hel <- as.data.frame(vegan::decostand(species, method = "hellinger"))
head(spe.hel)

# 3. Compute the Covariance Matrix

  # Compute covarience matrix by hand
comp_cov <- function(data) {
  output_df <- as.data.frame(matrix(NA, nrow = ncol(spe.hel), ncol = ncol(spe.hel)), row.names = colnames(data))
  colnames(output_df) <- colnames(data)
  for (i in 1:ncol(data)) {
    mean.dif.x1 <- data[,i] - mean(data[,i])
    for (p in 1:ncol(data)) {
      mean.dif.x2 <- data[,p] - mean(data[,p])
      output <- mean.dif.x1 * mean.dif.x2
      output_df[i, p] <- sum(output)/(nrow(data)-1)
    }
  }
  return(output_df)
}

cov_matrix <- comp_cov(spe.hel)

  # Compare with built-in function
cov_base_func <- cov(spe.hel)

# 4. Perform the Eigen-decomposition of the covariance matrix

eigen_decomp <- eigen(cov_matrix)

  # Extract Eigenvalues
(eig_values <- eigen_decomp$values)
  # Eigenvalues are equal to the sum of squares of the distances of each projected data point in the corresponding principal component
  # The sum of squares is maximized in the first principal component.

  # Extract Eigenvectors
(eig_vectors <- -eigen_decomp$vectors)
rownames(eig_vectors) = colnames(spe.hel)
  
  # Extract the first two eigenvectors
eig_vec_1 <- eig_vectors[,1]
eig_vec_2 <- eig_vectors[,2]

# 5. Show amount of variance contributed by each Principal Component by making a Scree Plot.

  # Calculate the estimated variance for each eigenvalue
(e_var <- eig_values / (nrow(spe.hel) - 1))

  # Data frame with variance percentages
var_per <- data.frame(
  PC  = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06","PC07", "PC08", "PC09",
          "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18",
          "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27"),
  PER = c(e_var) * 100 / sum(e_var)) # Calculate the percentage

  # Scree plot to show amount of variance accounted by each principal component
barplot(PER ~ PC, data = var_per,
        xlab = "Principal Components",
        ylab = "Percent of Variation %")

# 6. Selecting Principal Components 

  # Kaiser-Guttman Criterion
eig_val_PC <- data.frame(
  PC = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06","PC07", "PC08", "PC09",
         "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18",
         "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27"),
  EV = eig_values)
  
  # Plot principal components and overlay average eigenvalue line
barplot(EV ~ PC, data = eig_val_PC,
        xlab = "Principal Components",
        ylab = "Eigenvalues")
abline(h = mean(eig_values), col = "red")

  # Broken stick Model
broken_stick <- function(eig_values) {
  # Calculate Broken Stick Model
  n = length(eig_values)
  bsm = data.frame(j=seq(1:n), prop_var=0)
  bsm$prop_var[1] = 1/n
  for (i in 2:n) {
    bsm$prop_var[i] = bsm$prop_var[i-1] + (1/(n + 1 - i))
    }
  bsm$prop_var = 100*bsm$prop_var/n
  
  # Plot Broken Stick Modol Over 
  barplot(t(cbind(100*eig_values/sum(eig_values), bsm$p[n:1])),
        beside=TRUE, 
        main="Broken Stick Model",
        col=c("red","blue"),
        las=2,
        xlab = "Principal Components", ylab = "Percent of Variation (%)")
  legend("topright", c("Eigenvalues", "Broken stick model"), 
         pch=15,
         col=c("red","blue"), 
         bty="n")

}

broken_stick(eig_values)

# If we use the Kaiser-Guttman Criterion, we can how much how much variation is 
# captured by the 5 selected principal components.

eig_vec_kgc <- eig_values[eig_values > mean(eig_values)]

sum(var_per$PER[1:length(eig_vec_kgc)])


# 7. Plot the Principal Components over the data

  # Plot only the first principal component
plot(BAR ~ BLA, col = as.factor(rownames(spe.hel)), pch = 19,
     xlim = c(-0.25,0.5), ylim = c(-0.5,0.5),
     data = (spe.hel), xlab = "BLA (Standardized)", ylab = "BAR (Standardized)")
abline(v=0 , h=0, col = "dark gray")

  # Overlap first eigenvector/principal component
abline(0, eig_vec_1[11]/eig_vec_1[6], col='purple')

  # Plot lines from the first eigenvector to points
line1 <- c(0, eig_vec_1[11]/eig_vec_1[6])
perp.segment.coord <- function(x0, y0, line1){
  a <- line1[1]  #intercept 
  b <- line1[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}
ss <- perp.segment.coord(spe.hel[,6], spe.hel[,11], line1)
segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1, col = 'purple')
with(spe.hel, text(BAR ~ BLA, labels = as.factor(rownames(spe.hel)), pos = 1, cex=1))
title(main = "First Principal Component over the Standardized Data",
      sub = "Purple Lines Horizontal to the First Principal Components is the Variance", cex.sub = 0.75)

  # Plot both the first and second principal component
plot(BAR ~ BLA, col = as.factor(rownames(spe.hel)), pch = 19,
     xlim = c(-0.25,0.5), ylim = c(-0.5,0.5),
     data = (spe.hel), xlab = "BLA (Standardized)", ylab = "BAR (Standardized)")
abline(v=0 , h=0, col = "dark gray")

  # Overlap pertinent eigenvectors
abline(0, eig_vec_1[11]/eig_vec_1[6], col='purple')
abline(0, eig_vec_2[11]/eig_vec_2[6], col='orange')

  # Plot the lines from second eigenvector to points
line2 <- c(0, eig_vec_2[11]/eig_vec_2[6])
ss <- perp.segment.coord(spe.hel[,6], spe.hel[,11], line2)
segments(x0 = ss$x0, x1 = ss$x1, y0 = ss$y0, y1 = ss$y1,col = 'orange')
with(spe.hel, text(BAR ~ BLA, labels = as.factor(rownames(spe.hel)),pos = 1, cex=1))
title(main = "First (Purple) and Second (Orange) Principal Component over the Standardized Data",
      cex.main = 0.8, sub = "Lines Horizontal to the Principal Components are the Variance", cex.sub = 0.75)
# 8.Loading Scores

  # Elements of each eigenvector are called loadings and can be interpreted as the contribution of each variable in the data set to 
  # the corresponding principal component.

  # You can make a table with these values and see the contributions of each variable to each principal component:
  
  # Get variable loading scores
variable.loads <- data.frame(
  PC01 = eig_vec_1, # First eigenvector
  PC02 = eig_vec_2  # Second eigenvector
)
head(variable.loads)

  # You can also calculate the loading score for each site, which shows how they are placed in relation to the principal
  # components.

  # Get site loading scores
loading.scores <- as.data.frame(as.matrix(spe.hel) %*% eig_vectors)
colnames(loading.scores) = c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06","PC07", "PC08",
                             "PC09", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16",
                             "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24",
                             "PC25", "PC26", "PC27")
head(loading.scores)

# 9. Make a biplot of the Principal Components and Loading Scores

  # Set plot parameters
par(mar = c(5, 5, 10, 5),
     mgp = c(2, 1, 0))

  # Plot site loading scores
plot(loading.scores[,2] ~ loading.scores[,1], 
     xlab = 'PC1', ylab = "PC2",
     xlim = c(-1,1), ylim  = c(-1,1),col = as.factor(rownames(loading.scores)), pch = 19)
abline(v = 0, col = "orange")
abline(h = 0, col = "purple")
with(loading.scores, text(PC02 ~ PC01, labels = as.factor(rownames(loading.scores)),pos = 1, cex=1))

par(new=TRUE)

  # Overlay the variable loading scores
plot(PC02 ~ PC01, 
     xlim = c(-1, 1), ylim = c(-1,1),
     col = "red", pch = 8, axes = F, xlab = "", ylab = "",
     data = variable.loads)
axis(4, ylim = c(-1,1), col = "red")
axis(3, xlim = c(-0.9,1), col = "red")
mtext("Loading Scores",side=3,col="red",line=2.5)  
mtext("Loading Scores",side=4,col="red",line=2.5)  
for (i in 1:nrow(variable.loads)) {
  arrows(x0 = 0, y0 = 0, x1 = variable.loads[i,1],y1 = variable.loads[i,2],
         col = "red", lwd = 1, length = 0.1)
}
with(variable.loads, text(PC02 ~ PC01, labels = as.factor(rownames(variable.loads)),pos = 1, cex=1,
                 col = "red"))
title(main = "Biplot of Site and Variable Loading Scores against the First and Second Principal Components",
      cex.main = 0.9)

### PCA analysis using built-in functions ----

par(mar = c(5, 4, 4, 2) + 0.1,mgp = c(3, 1, 0))

  # Using stats: prcomp()
PCA_prcomp <- stats::prcomp(spe.hel)
biplot(PCA_prcomp, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5))
abline(v= 0, h = 0)


  # Using stats: princomp()
PCA_princomp <- stats::princomp(spe.hel)
biplot(PCA_princomp, xlim = c(-0.5,0.5), ylim = c(-0.5,0.5))
abline(v= 0, h = 0)

  # Using vegan::rda()
PCA_rda <- vegan::rda(spe.hel)
biplot(PCA_rda, scaling = 2)
biplot(PCA_rda, scaling = 1)

### Scaling

  # Type 2 scaling
biplot(PCA_rda, scaling = 2)

  # Type 1 scaling
biplot(PCA_rda, scaling = 1)


## Correspondance Analysis ----

  # Run the CA using the cca() function in the vegan package

    # Load data
data(Doubs)
species <- Doubs.fish[-8,]

    # Run CA using the vegan package
spe.ca <- vegan::cca(species)

    # Identify the eigenvectors using the Kaiser-Guttman Criterion
eig_vec_ca <- spe.ca$CA$eig
(ev_kgc <- eig_vec_ca[eig_vec_ca > mean(eig_vec_ca)])

    # Scree Plot
barplot(eig_vec_ca)
abline(h = mean(eig_vec_ca), col = "red")

  # Plot CA
plot(spe.ca, scaling = 2, type = "none", main = "CA",
     xlab = c("CA1"), ylab = c("CA2"))

points(vegan::scores(spe.ca, display = "sites", choices = c(1, 2), scaling = 2),
       pch = 21, col = "black", bg = "steelblue", cex = 1.2)

text(vegan::scores(spe.ca, display = "species", choices = c(1), scaling = 2),
     vegan::scores(spe.ca, display = "species", choices = c(2), scaling = 2),
     labels = rownames(scores(spe.ca, display = "species", scaling = 2)),
     col = "red", cex = 0.8)

## Principal Coordinate Analysis ----

  # Use the pcoa() function in the ape package to run a PCoA
spe.h.pcoa <- ape::pcoa(dist(spe.hel))
    # Plot PCoA biplot
ape::biplot.pcoa(spe.h.pcoa, spe.hel)

  # Run a PcoA with Bray-Curtis distance
    # Calculate Bray-Curtis coefficients
spe.bc <- vegan::vegdist(species, method = "bray")
    # Run PCoA
spe.bc.pcoa <- ape::pcoa(spe.bc)
    # Plot PCoA biplot
ape::biplot.pcoa(spe.bc.pcoa, spe.hel, dir.axis2 = -1)

## Nonmetric MultiDimensional Scaling ----

  # Load Data
data(Doubs)
species <- Doubs.fish[-8,]

  # Apply Helligner transformation
spe.h <- decostand(species, method = "hellinger")

  # Run the NMDS
spe.nmds <- vegan::metaMDS(spe.h,
                           distance = "bray", # specify the distance coefficient to be calculated in the function
                           k = 2, # specify the number of dimensions
                           autotransform = F # indicates that transformation has already been applied to the data
)

  # Stress score
spe.nmds$stress # A stress score <0.1. Great!

  # Shepard Plot
vegan::stressplot(spe.nmds, main = "Shepard plot")

  # NMDS Biplot
vegan::ordiplot(spe.nmds, type = "t", xlim = c(-2,2),
                ylim = c(-1.5, 1.5), main = "NMDS with Doubs Fish Species Community Abundance Data", sub = "Hellinger Transformed data with Bray-Curtis Coefficients")
abline(v = 0, h = 0)