##For PC users
install.packages("INLA",
                 repos=c(getOption("repos"),
                 INLA="https://inla.r-inla-download.org/R/stable"), 
                 dep=TRUE)
# Core packages required for this tutorial

library(dplyr)
library(geoR)
library(sf)
library(leaflet)
library(viridis)
library(terra)
library(geodata)
library(INLA)
library(geostats)
library(fmesher)
library(inlabru)
# Load the Gambia dataset (included with INLA)
data(gambia)

# Explore the data structure
str(gambia)
head(gambia)

# Quick summary statistics:
cat("\nNumber of observations:", nrow(gambia), "\n")
cat("Number of villages:", length(unique(gambia$x)), "\n")
cat("Malaria prevalence:", mean(gambia$pos) * 100, "%\n")

# Map of malaria prevelance by village
village_prev <- gambia %>%
  group_by(x, y) %>%
  summarise(total = n(),
    positive = sum(pos),
    prev = positive / total
  )

#### Step 1: convert UTM to WGS84 ####

# Gambia is in UTM zone 28N (32628)
village_sf <- st_as_sf(village_prev, 
                       coords = c("x", "y"), 
                       crs = 32628)  # UTM Zone 28N

# Projecting coordinates to WGS84
village_latlong <- st_transform(village_sf, crs = 4326)

##### Step 2: Extract coordinates ####
coords_latlong <- st_coordinates(village_latlong)
village_prev$lon <- coords_latlong[, 1]
village_prev$lat <- coords_latlong[, 2]

#### Step 3: Create color palette ####
pal <- colorNumeric(palette = viridis(100), 
                    domain = village_prev$prevalence)

#### Step 4: Create leaflet map ####
#?leaflet
# If you are curious about the different leaflet tiles, check out: https://leaflet-extras.github.io/leaflet-providers/preview/

leaflet(village_prev) %>%
  addTiles() %>%  # this is where you would change your background map
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    radius = 8,
    color = ~pal(prev),
    fillColor = ~pal(prev),
    fillOpacity = 0.8,
    stroke = TRUE,
    weight = 1,
    popup = ~paste0("Prevalence: ", round(prev * 100, 1), "%")
  ) %>%
  addLegend(
    position = "topright",
    pal = pal,
    values = ~prev,
    title = "Malaria<br>Prevalence",
    labFormat = labelFormat(transform = function(x) round(x * 100, 0),
                           suffix = "%")
  ) %>%
  addScaleBar(position = "bottomleft")

# Extract unique village coordinates
coords_utm <- unique(gambia[, c("x", "y")])

# Check the scale of your data to use in creating the mesh
range(coords_utm$x)  
range(coords_utm$y)
# creating the mesh
mesh <- inla.mesh.2d(loc = coords_utm,
                     max.edge = c(10000, 50000),  # 10km inner, 50km outer
                     cutoff = 5000,                # 5km minimum distance
                     offset = c(10000, 30000))     # Extensions

# Visualize the mesh
plot(mesh, main = "Spatial Mesh for Gambia Data (UTM)")
points(coords_utm$x, coords_utm$y, pch = 19, col = "red", cex = 0.5)
##IF THIS DOESN"T PLOT: paste these commands into the consule and run it

# Check mesh size
mesh$n

# What makes a good mesh?

# We can get the elevation raster for Gambia using the geodata package
# for a 1km resolution, we can use elevation_30s

r <- elevation_30s(country = "GMB", path = tempdir())
r <- terra::project(r, "EPSG:32628")
pal <- colorNumeric("viridis", values(r),
                    na.color = "transparent"
)
# map elevation raster
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal, values = values(r),
            title = "Altitude (m)"
  ) %>%
  addScaleBar(position = c("bottomleft"))
#What makes a good mesh? Can you include barriers?
###

# extract elevation for all observations and add "alt" to dataframe
village_prev["alt"] <- 
  terra::extract(r, village_prev[, c("x", "y")], 
                 list=T, ID=F, method="bilinear")
head(village_prev)


# What is this step doing? Why is it important?
##

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
#alpha = 2 is the smoothing parameter (default choice)


# Create index for spatial field...a list object INLA needs to run
spatial.field <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

##Code for adding covariate to a mesh

##slope_on_mesh <- eval_spatial( my_raster , mesh$loc[,1:2])
##spde.slope <- inla.spde2.matern(mesh, alpha = 2,          ##                B.tau = cbind(0, 1, slope_on_mesh, 0, 0), 
##                B.kappa = cbind(0, 0, 0, 1,slope_on_mesh))
#indexs.slope <- inla.spde.make.index("s", spde.slope$n.spde)


A <- inla.spde.make.A(mesh = mesh, 
                      loc = as.matrix(village_prev[, c("x", "y")]))

dim(A)

ra <- terra::aggregate(r, fact = 4, fun = mean) # reduce the number of raster cells, factor 4 combines 4x4 cells of raster into one cell
dp <- terra::as.points(ra) # take aggregated raster and turn into vector of points so it is easier to get coordinates from them

# then use the crds() function to get coordinates and put everything into a matrix
dp <- as.matrix(cbind(crds(dp)[,1], crds(dp)[,2], values(dp)))
colnames(dp) <- c("x", "y", "alt")
head(dp)
# Prepare covariates

# Create the data stack
stk.e <- inla.stack(
  tag = "est", # lets INLA know this is our estimation stack
  data = list(y = village_prev$positive, numtrials = village_prev$total), 
  A = list(1, A), # A, projection matrix
  effects = list(data.frame(b0 = 1, 
                            altitude = village_prev$alt), 
                 s = spatial.field)
)

dim(dp)
coop <- dp[, c("x", "y")] # prediction coordinates from the raster

# make the prediction matrix
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)



#estimation stack
stk.e <- inla.stack(
  tag = "est", # lets INLA know this is our estimation stack
  data = list(y = village_prev$positive, numtrials = village_prev$total), 
  A = list(1, A), # A, projection matrix
  effects = list(data.frame(b0 = 1, 
                            altitude = village_prev$alt), s = spatial.field)
)
# prediction stack
stk.p <- inla.stack(
  tag = "pred", # tag it for prediction, so INLA knows
  data = list(y = NA, numtrials = NA),
  A = list(1, Ap), # Ap, prediction matrix
  effects = list(data.frame(b0 = 1, 
                            altitude = dp[, 3]),
                 s = spatial.field
  )
)

stk.full <- inla.stack(stk.e, stk.p)  # assembles the data for INLA, similar to how we compile a STAN model
formula <- y ~ 0 + b0 + altitude + f(spatial.field, model = spde) 
#s is the name we gave the spatial model; spde is the object name with the model

res <- inla(formula,
            family = "binomial",
            Ntrials = numtrials,
            control.family = list(link = "logit"),
            control.compute=list(return.marginals.predictor=TRUE, waic=TRUE),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, # this computes the posteriors of the predictions
              link = 1,
              A = inla.stack.A(stk.full)
            )
)

summary(res)


index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

# Convert UTM coordinates to lat/long for leaflet
coop_sf <- st_as_sf(data.frame(x = coop[, 1], y = coop[, 2]), 
                    coords = c("x", "y"), 
                    crs = 32628)  # UTM Zone 28N
coop_latlong <- st_transform(coop_sf, crs = 4326)
coop_coords <- st_coordinates(coop_latlong)

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    lng = coop_coords[, 1],  # Use converted longitude
    lat = coop_coords[, 2],  # Use converted latitude
    color = pal(prev_mean),
    fillColor = pal(prev_mean),
    fillOpacity = 0.7
  ) %>%
  addLegend("bottomright",
            pal = pal, values = prev_mean,
            title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))

### Rasterize the Prediction ----
r_prev_mean <- terra::rasterize(
  x = coop, y = ra, values = prev_mean,
  fun = mean
)

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_mean), title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))
# Transform spatial hyperparameters to interpretable scale
# Range: distance at which correlation drops to ~0.1
# Variance: magnitude of spatial variation

spde.result <- inla.spde2.result(res, "spatial.field", spde)

# Posterior mean of range (in coordinate units)
range.mean <- inla.emarginal(function(x) sqrt(8)/exp(x), 
                              spde.result$marginals.log.kappa[[1]])

# Posterior mean of variance
var.mean <- inla.emarginal(function(x) 1/exp(x),
                            spde.result$marginals.log.tau[[1]])

cat("\nSpatial Range (posterior mean):", round(range.mean, 3), "units\n")
cat("Spatial Variance (posterior mean):", round(var.mean, 3), "\n")


#help(package="geodata")

gdata <- worldclim_country(country = "GMB",
                               var="prec", # other options are tmin, tmax, tavg, prec, wind, and bio
                               path = tempdir(),
                               version=2.1,
                               mask = TRUE)
# a quick google tells us that malaria is most present during and immediately after the rainy season, from June-October,
## so lets subset our data!

my_raster <- mean(subset(x=gdata, subset=c(6,7,8,9,10))) #subsetting June-October and taking the mean 
my_raster <- terra::project(my_raster, "EPSG:32628")
# Now extract using UTM coordinates
village_prev["prec"] <- 
  terra::extract(my_raster, village_prev[, c("x", "y")],  # Use x, y instead of lon, lat
                 list = TRUE, ID = FALSE, method = "bilinear")
head(village_prev)
pal <- colorNumeric("viridis", values(my_raster),
                    na.color = "transparent"
)
# map my_raster
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(my_raster, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal, values = values(my_raster),
            title = "Precip (mm)" # change to the units for your variable
  ) %>%
  addScaleBar(position = c("bottomleft"))


# extract measurement for all observations and add this variable to dataframe
village_prev["prec"] <- # change this to reflect the variable you chose
  terra::extract(my_raster, village_prev[, c("lon", "lat")],
                 list=T, ID=F, method="bilinear")
head(village_prev)

# stk.p <- inla.stack(
#   tag = "pred", # tag it for prediction, so INLA knows
#   data = list(y = NA, numtrials = NA),
#   A = list(1, Ap),
#   effects = list(data.frame(b0 = 1, 
#                             altitude = raster.coord[,3]),
#                  s = spatial.field
#   )
# )

# First, extract precipitation for the prediction grid points
prec_raster <- terra::aggregate(my_raster, fact = 4, fun = mean)
dp_prec <- terra::extract(prec_raster, dp[, c("x", "y")], method = "bilinear")
dp <- cbind(dp, prec = dp_prec[, 1])

# Estimation stack (uses village data - this part was correct)
stk.e1 <- inla.stack(
  tag = "est",
  data = list(y = village_prev$positive, numtrials = village_prev$total), 
  A = list(1, A),
  effects = list(data.frame(b0 = 1, 
                            altitude = village_prev$alt,
                            precip = village_prev$prec), 
                 s = spatial.field)
)

# Prediction stack (now uses dp data - this was the problem)
stk.p1 <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, 
                            altitude = dp[, "alt"],
                            precip = dp[, "prec"]),
                 s = spatial.field)
)

# Combine stacks
stk.full1 <- inla.stack(stk.e1, stk.p1)
