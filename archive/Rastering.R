
### Rastering an image with the Terra package-----------------

##Load correct packages
library(terra)

##download the photo you would like to raster 
drone<-terra::rast("SWRS_Plot771_1.JPG")

###Plot-------------------------
plot(drone, y=1, main="Elevational Gradient", grid=FALSE, breaks=NULL, legend=TRUE, maxcell=500000)


### Different forms of Rastering -----------------
#Information from: https://rspatial.github.io/terra/reference/plot.html

## SpatRaster
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
plot(r)

#Plotting at intervals 
plot(r, type="interval")

#Plotting with coordinates 
e <- c(6.37, 6.41, 49.9, 50.1)
plot(r, plg=list(ext=e, title="Legend\nTitle", title.cex=0.9), 
    pax=list(side=1:4, retro=TRUE))
north(cbind(5.8, 50.1))  


### Extra cool things that can be done with the new terra package -----------------
#Information from: https://dahtah.github.io/imager/imager.html
library(imager)

file <- system.file('extdata/parrots.png',package='imager')
#system.file gives the full path for a file that ships with a R package
#if you already have the full path to the file you want to load just run:
#im <- load.image("/somedirectory/myfile.png")
im <- load.image(file)

plot(im) #Parrots!

im.blurry <- isoblur(im,10) #Blurry parrots!
plot(im.blurry)

im.xedges <- deriche(im,2,order=2,axis="x") #Edge detector along x-axis
plot(im.xedges)

im.yedges <- deriche(im,2,order=2,axis="y") #Edge detector along y-axis
plot(im.yedges)


file <- system.file('extdata/parrots.png',package='imager')
parrots <- load.image(file)
#The file format is defined by the extension. Here we save as JPEG
imager::save.image(parrots,"parrots.jpeg")
#We call imager::save.image to avoid ambiguity, as base R already has a save.image function 

