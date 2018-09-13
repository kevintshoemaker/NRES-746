##################################################
####                                          ####  
####  R Bootcamp #2, Module X                 ####
####                                          #### 
####   University of Nevada, Reno             ####
####                                          #### 
##################################################

###################################################
####  Advanced data visualization with ggplot  ####
###################################################

#Load up the necessary libraries
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(carData)
library(DAAG)
library(RColorBrewer)
library(leaflet)

#Load the dataset and see what variables it contains
soil <- data.frame(Soils)
head(soil)

#Basic boxplot
ggplot(soil) +
  geom_boxplot(aes(x=Contour, y=pH))

#Basic scatterplot
ggplot(soil) +
  geom_point(aes(x=pH, y=Ca))

#Add color
ggplot(soil) +
  geom_point(aes(x=pH, y=Ca, color=Depth))

#More modifications to points
ggplot(soil) +
  geom_point(aes(x=pH, y=Ca, fill=Depth), shape=21, color="black", size=4, stroke=1.5)

#Plot multiple series
ggplot(soil, aes(x=pH)) +
  geom_point(aes(y=Ca), shape=21, fill="red", color="black", size=4, stroke=1.5) +
  geom_point(aes(y=Mg), shape=21, fill="blue", color="black", size=4, stroke=1.5) +
  geom_point(aes(y=Na), shape=21, fill="gray30", color="black", size=4, stroke=1.5)

#Another way
soil.nut <- gather(soil, nutrient, value, c(10,11,13))
ggplot(soil.nut) +
  geom_point(aes(x=pH, y=value, fill=nutrient), shape=21, color="black", size=4, stroke=1.5)

#Using different set of variables
soil.nut2 <- gather(soil, nutrient, value, c(10,11,12))
ggplot(soil.nut2) +
  geom_point(aes(x=pH, y=value, fill=nutrient), shape=21, color="black", size=4, stroke=1.5)

#Add facets and change theme elements
ggplot(soil.nut2) +
  geom_point(aes(x=pH, y=value, fill=nutrient), 
             shape=21, color="black", size=4, stroke=1.5) +
  facet_wrap(~nutrient, scales="free_y") +
  ylab("mg / 100 g soil") +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face="bold"))

#check out RColorBrewer options
display.brewer.all()

#Plot with RColorBrewer palette
ggplot(soil) +
  geom_point(aes(x=pH, y=Ca, fill=Depth), shape=21, color="black", size=4, stroke=1.5) +
  theme_bw() +
  ylab("Ca (mg/100g soil)") +
  scale_fill_brewer(palette="YlOrBr", name="Depth (cm)")

#Plot with my hand selected palette
ggplot(soil) +
  geom_point(aes(x=pH, y=Ca, fill=Depth), shape=21, color="black", size=4, stroke=1.5) +
  theme_bw() +
  ylab("Ca (mg/100g soil)") +
  scale_fill_manual(values=c("#FFF0BF","#FFC300","#BF9200","#604900"), name="Depth (cm)")

#Plot with trendlines
ggplot(soil.nut2) +
  geom_point(aes(x=pH, y=value, fill=nutrient), 
             shape=21, color="black", size=4, stroke=1.5) +
  geom_smooth(aes(x=pH, y=value), method="lm", color="black") +
  facet_wrap(~nutrient, scales="free_y") +
  ylab("mg / 100 g soil") +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face="bold"))

#Histogram
ggplot(soil.nut) +
  geom_histogram(aes(x=value), color="black", fill="white", bins=15) +
  facet_wrap(~nutrient, scales="free") +
  xlab("mg / 100g soil") +
  theme_dark() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face="bold"))

#Histogram with density curves overlaid
ggplot(soil.nut) +
  geom_histogram(aes(x=value, y=..density..), color="black", fill="white", bins=15) +
  geom_density(aes(x=value,color=nutrient), size=1.5) +
  facet_wrap(~nutrient, scales="free") +
  xlab("mg / 100g soil") +
  theme_dark() +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face="bold"))

#Histogram with normal distribution overlaid
ggplot(soil.nut) +
  geom_histogram(aes(x=value, y=..density..), color="black", fill="white", bins=15) +
  stat_function(fun = dnorm, color = "blue", size = 1.5,
                args=list(mean=mean(soil.nut$value), sd=sd(soil.nut$value))) +
  facet_wrap(~nutrient, scales="free") +
  xlab("mg / 100g soil") +
  theme_dark() +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face="bold"))

#Boxplot with means and error bars in place of capless whiskers
ggplot(soil, aes(x=Contour, y=pH)) +
  stat_boxplot(geom="errorbar", width=0.2) +  
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", size=5, color="black")

#Leaflet showing possum sites
leaflet(possumsites) %>%
  addTiles() %>% #Adds map tiles from OpenStreetMap
  addMarkers(lng=c(possumsites$Longitude), lat=c(possumsites$Latitude), 
             popup=c(as.character(possumsites$altitude))) #Adds markers for the sites