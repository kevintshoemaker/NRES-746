# gam lab- worked

# sage: default- thin plate spline (for GAMs) (s(x, bs="tp"))
   # cubic regression splines- basis functions are cubic polynomial. 
   # different splines fit a little different. 


# Q: what if you want to model a random effect
# p-splines are a little more flexible in some cases.  


library(mgcv) #load the mgcv package


# demo -------------------

mcycle <- MASS::mcycle #pull in some data; we'll be using this dataset for this example but not the lab

# ?gam #take a look at the documentation page for the function to see what you can specify - it's a lot!

#now we'll go ahead and fit our spline using the most basic call of the function possible
gam_mod <- gam(accel ~ s(times), data = mcycle) #here, we predict acceleration based on a smooth, nonlinear function of times

plot(gam_mod, residuals = T, pch = 1) #setting residuals=T will include the CIs on the plot
coef(gam_mod)
gam_mod_s1 <- gam(accel ~ s(times), data = mcycle, sp = 0.1) #smoothing by setting a fixed smoothing parameter
gam_mod_s2 <- gam(accel ~ s(times), data = mcycle, method = "REML")


par(mfrow = c(2, 1)) #we'll plot both models side by side to compare
plot(gam_mod_s1, residuals = TRUE, pch = 1)
plot(gam_mod_s2, residuals = TRUE, pch = 1)

gam_mod_tp <- gam(accel ~ s(times, bs="tp"), data = mcycle)
plot(gam_mod_tp, residuals = T, pch = 1)
gam_mod_cr <- gam(accel ~ s(times, bs="cr"), data = mcycle)
plot(gam_mod_cr, residuals = T, pch = 1)

#you can leave the knots argument blank, as above, or specify it:
gam_mod_cr_3 <- gam(accel ~ s(times, bs="cr",k=3), data = mcycle) #(in this example this is clearly a terrible choice of number of knots)
plot(gam_mod_cr_3, residuals = T, pch=1)  
gam_mod_ps <- gam(accel ~ s(times, bs="ps"), data = mcycle)
plot(gam_mod_ps, residuals = T, pch = 1)
library(sp)
data(meuse, package="sp") #load in and inspect dataset
head(meuse)

mod <- gam(cadmium ~ s(x, y) + s(elev), #first, we fit a model where x and y coordinates interact directly and
           data = meuse, method = "REML") #elevation is separate
plot(mod)
tensor_mod <- gam(cadmium ~ te(x, y, elev),  #this model allows x and y coordinates and elevation to interact 
                  data = meuse, method = "REML")  #despite being on different scales
plot(tensor_mod)
tensor_mod2 <- gam(cadmium ~ s(x, y) + s(elev) + ti(x, y, elev), #notice how two regular smooths are fit in addition to
                   data = meuse, method = "REML")                   #the tensor interaction
mod2da <- gam(cadmium ~ s(dist)+landuse, 
              data = meuse, method = "REML")
mod_sep <- gam(copper ~ s(dist, by = landuse) + landuse, #notice how landuse is included both inside and outside the smooth
               data = meuse, method = "REML") #this gives us different smooths and different intercepts for each value of                                                     the variable
summary(mod_sep)
plot(mod_sep,pages=2)
mod2d <- gam(cadmium ~ s(x,y), data = meuse, method = "ML")
mod2d_elev <- gam(cadmium ~ s(x, y) + s(elev), 
                  data = meuse, method = "ML")

library(lmtest) #we need to load this package to do the likelihood ratio test
lrtest(mod2d_elev,mod2d) 
mod <- gam(cadmium ~ s(x, y) + s(elev), #this model considers elevation separately from x and y
           data = meuse, method = "ML")

tensor_mod <- gam(cadmium ~ te(x, y, elev), #this model considers elevation interacting with x and y
                  data = meuse, method = "ML")

AIC(mod,tensor_mod) # we can calculate both aic values simultaneously to make comparison even easier!
mod2d <- gam(cadmium ~ s(x,y), data = meuse, method = "REML")

vis.gam(mod2d, view = c("x", "y"), #the view argument allows you to define which variables make the axes
        plot.type = "contour", too.far = 0.05)
points(meuse) 

#you may need to specify which parts of the data to pass to points, as below (sometimes it just picks the first two columns which is an issue if that's not the part of your model you're plotting)
points(meuse$x,meuse$y)


# exercise ---------

df <- read.csv("GAMs_SnakeSpeeds.csv")

library(mgcv)

df <- subset(df,Speed>0) # clean data
df <- subset(df,SVL>20)

hist(df$Speed)

names(df)
plot(Speed~Mass,df) # weak pos
plot(Speed~Tail,df) # pos
plot(Speed~Lat,df)  # pos
plot(Speed~Long,df)  # little effect
plot(Speed~SVL,df)    # outliers in the svl? what's going on here?
plot(Speed~MAMU,df)

table(df$Age)  # age is categorical

tapply(df$Speed,df$Age,mean) # not much effect of age


plot(Speed~Age,df)


library(sf)
df2 <- st_as_sf(df,coords=c("Lat","Long"),crs=4326)
plot(df2)

names(df)



# run GAMs ----------------

mod1 <- gam(Speed~te(Lat,Long),data=df)

plot(mod1)
vis.gam(mod1)
vis.gam(mod1,view = c("Lat", "Long"), #the view argument allows you to define which variables make the axes
     plot.type = "contour")

mod2 <- gam(Speed~s(Tail),data=df)
plot(mod2)
summary(mod2)
coef(mod2) # quite a few knots

names(df)
mod3 <- gam(Speed~s(SVL) + s(Mass) + s(Tail) + Lat + MAMU,data=df)

plot(mod3)
hist(residuals.gam(mod3))

summary(mod3)  # tail not significant, lat not significant

mod4 <- gam(Speed~s(SVL) + s(Mass) + s(MAMU),data=df)

plot(mod4)


mod5 <- gam(Speed~te(SVL,Mass) + s(MAMU) ,data=df )
plot(mod5)

summary(mod5)

AIC(mod4,mod5)


mod5 <- gam(Speed~te(SVL,Mass) + s(MAMU) ,data=df )
summary(mod5)



vis.gam(mod5,view = c("SVL", "Mass"), #the view argument allows you to define which variables make the axes
        plot.type = "contour")


