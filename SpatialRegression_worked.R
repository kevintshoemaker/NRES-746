
# Spatial regression topic... 


# Install packages --------------------

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)


# Load packages --------------------

library(geoR) # geoR may be obsolete soon
library(dplyr)
library(sf)
library(leaflet)
library(viridis)
library(terra)
library(geodata)
library(INLA)
library(geostats)


# Load data -------------------------

data(gambia)
head(gambia)

dim(gambia)   # 2035 observations at individual level
dim(unique(gambia[, c("x", "y")]))   # 65 unique locations

rm(gambia.borders)

# Process data -------------------

# For the lab-part 1, different variables, consider other models----
# dplyr and summary statistics to create the new columns
df <- group_by(gambia, x, y) %>%
  summarize(
    total = n(),
    positive = sum(pos),
    prev = positive / total,
    netuse = mean(netuse),   # for the lab, you could include these 
    greenness = mean(green)
  )
head(df)


# create geometric points from XY
pts <- 
  st_multipoint(x = as.matrix(df[,1:2]), dim="XY")
# create spatially referenced points: these are in UTM Zone 28 N
sp <- 
  st_sfc(pts, crs = "+proj=utm +zone=28")
# transform to long lat crs: WGS84
sp_tr <- 
  st_transform(sp, crs = "+proj=longlat +datum=WGS84")


# add long, lat columns to data frame
df[, c("long", "lat")] <- st_coordinates(sp_tr)[,1:2]
head(df)


# get raster data -----------------

elev <- elevation_30s(country = "GMB", 
                   path = tempdir(), mask = TRUE)

precip <- worldclim_country(country = "GMB",
                               var="prec", # choose anything!
                               path = tempdir(),
                               version=2.1,
                               mask = TRUE)

precip <- precip$GMB_wc2.1_30s_prec_10
# plot(precip)


# visualize spatial data -------------

pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(df) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~long, lat = ~lat, color = ~ pal(prev)) %>%
  addLegend("bottomright",
            pal = pal, values = ~prev,
            title = "Malaria Prevalence"
  ) %>%
  addScaleBar(position = c("bottomleft"))

pal <- colorNumeric("viridis", values(elev),
                    na.color = "transparent"
)

# map elevation raster
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(elev, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal, values = values(elev),
            title = "Altitude (m)"
  ) %>%
  addScaleBar(position = c("bottomleft"))

pal <- colorNumeric("viridis", values(precip),
                    na.color = "transparent"
)
# map precip raster
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(precip, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal, values = values(precip),
            title = "Precip"
  ) %>%
  addScaleBar(position = c("bottomleft"))

# extract elevation for all observations and add "alt" to dataframe

df["elev"] <- 
  terra::extract(elev, df[, c("long", "lat")], 
                 list=T, ID=F, method="bilinear")
head(df)

# extract precip for all observations and add "alt" to dataframe

df["precip"] <- 
  terra::extract(precip, df[, c("long", "lat")], 
                 list=T, ID=F, method="bilinear")
head(df)


# start the INLA spatial regression process -------------------

## Create the mesh network (nodes) --------

A <- inla.spde.make.A(mesh = mesh, loc = coo) # coo is coordinates of our observations


### Prediction data ----

dp <- terra::as.points(r) # this makes the r raster into a vector of points
dim(dp)

ra <- terra::aggregate(r, fact = 4, fun = mean) # reduce the number of raster cells, factor 4 combines 4x4 cells of raster into one cell
dp <- terra::as.points(ra) # take aggregated raster and turn into vector of points so it is easier to get coordinates from them
# then use the crds() function to get coordinates and put everything into a matrix
dp <- as.matrix(cbind(crds(dp)[,1], crds(dp)[,2], values(dp)))
colnames(dp) <- c("x", "y", "alt")

## For the lab-part 3 ----
# run the next 2 lines to include another environmental covariate raster
#dp <- as.matrix(cbind(dp, terra::extract(my_raster, dp[, c("x", "y")], method="bilinear")))

# for the lab-part 3: # change "wind" to whatever you used
#colnames(dp) <- c("x", "y", "alt", "wind") 

dim(dp)

coop <- dp[, c("x", "y")] # prediction coordinates from the raster

# make the prediction matrix
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
stk.full <- inla.stack(stk.e, stk.p)  # assembles the data for INLA, similar to how we "package" the data for JAGS      
formula <- y ~ 0 + b0 + altitude + f(s, model = spde) # for my example, I would add 'wind'
res <- inla(formula,
            family = "binomial",
            Ntrials = numtrials,
            control.family = list(link = "logit"),
            control.compute=list(return.marginals.predictor=TRUE),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, # this computes the posteriors of the predictions
              link = 1,
              A = inla.stack.A(stk.full)
            )
)
summary(res)
## Mapping Malaria Prevalence ----
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    lng = coop[, 1], lat = coop[, 2],
    color = pal(prev_mean)
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
### Lower limits of prediction
r_prev_ll <- terra::rasterize(
  x = coop, y = ra, values = prev_ll,
  fun = mean
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ll, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_ll), title = "LL"
  ) %>%
  addScaleBar(position = c("bottomleft"))
### Lower limits of prediction
r_prev_ul <- terra::rasterize(
  x = coop, y = ra, values = prev_ul,
  fun = mean
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ul, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_ul), title = "UL"
  ) %>%
  addScaleBar(position = c("bottomleft"))
## For the lab-part 4 ----
### Mapping exceedance probabilities ----
threshold <- 0.2 # exceedance threshold 20%

index <- inla.stack.index(stack = stk.full, tag = "pred")$data
marg <- res$marginals.fitted.values[index][[1]]
1 - inla.pmarginal(q = threshold, marginal = marg) # probability
excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = threshold, marginal = marg)})

head(excprob)
r_excprob <- terra::rasterize(
  x = coop, y = ra, values = excprob,
  fun = mean
)
# map
pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_excprob, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_excprob), title = "P(p>0.2)"
  ) %>%
  addScaleBar(position = c("bottomleft"))
