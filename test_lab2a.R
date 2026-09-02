

df = emdbook::ReedfrogSizepred
df

ricker = function(x,a,b) a*x*exp(-b*x)
sse = function(o,e) sum((o-e)^2)


apars = seq(0.1,3,length=100)     # seq(0.1,3,by=0.05) # 
bpars = seq(0.02,0.4,length=100)  # seq(0.02,0.4,by=0.01)  # 

gd = expand.grid(apars,bpars)

names(gd) = c("a","b")

gd$sse  =  numeric(nrow(gd))

i=100
for(i in 1:nrow(gd)){
  gd$sse[i] = sse(df$Kill,ricker(df$TBL,gd$a[i],gd$b[i])) 
}

ndx=which.min(gd$sse)

gd[ndx,]

plot(df$TBL,df$Kill,xlim=c(0,50))
curve(ricker(x,gd$a[ndx],gd$b[ndx]),add=T)

gd$logsse = log(gd$sse)
gd2 <- subset(gd, logsse <= 3.8)

library(ggplot2)
ggplot(gd2, aes(x = a, y = b, z = log(sse))) +
  # Create the heatmap background
  geom_raster(aes(fill = log(sse)), interpolate = TRUE) +
  # Add the contour lines
  geom_contour(color = "white", linewidth = 0.5) +
  # Apply a smooth, modern color palette
  scale_fill_viridis_c(option = "magma", direction = -1) +
  # Clean up the presentation
  theme_minimal() +
  labs(
    title = "Tadpole example",
    x = "a",
    y = "b",
    fill = "SSR"
  )




