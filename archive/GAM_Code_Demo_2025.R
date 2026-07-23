#install.packages("mgcv")

library(ggplot2)
library(mgcv)
isit <- read.csv("ISIT.csv")
#focusing on Season 2 to start
isit2 <- subset(isit, Season == 2)
head(isit2)
linear_model <- gam(Sources ~ SampleDepth, data = isit2)
summary(linear_model)
data_plot <- ggplot(data = isit2, aes(y = Sources, x = SampleDepth)) +
    geom_point() + geom_line(aes(y = fitted(linear_model)), colour = "red",
    linewidth = 1.2) + theme_bw()
data_plot
?smooth.terms
?mgcv::gam
gam_model <- gam(Sources ~ s(SampleDepth), bs = 'tp', method = "REML", data = isit2) #s() is the smooth term

data_plot <- data_plot + geom_line(aes(y = fitted(gam_model)),
    colour = "blue", linewidth = 1.2)
data_plot
plot(gam_model)
linear_model <- gam(Sources ~ SampleDepth, data = isit2)
smooth_model <- gam(Sources ~ s(SampleDepth), data = isit2)
AIC(linear_model, smooth_model)
head(isit)
isit$Season <- as.factor(isit$Season)
levels(isit$Season)
basic_model <- gam(Sources ~ Season + s(SampleDepth), data = isit,
    method = "REML") # notice the s() for SampleDepth, but not for Season
basic_summary <- summary(basic_model)
basic_summary$p.table
basic_summary$s.table
two_term_model <- gam(Sources ~ Season + s(SampleDepth) + RelativeDepth,
    data = isit, method = "REML")
two_term_summary <- summary(two_term_model)
two_term_summary$p.table
two_term_summary$s.table
par(mfrow = c(2, 2))
plot(two_term_model, all.terms = TRUE)
two_smooth_model <- gam(Sources ~ Season + s(SampleDepth) + s(RelativeDepth),
    data = isit, method = "REML")
two_smooth_summary <- summary(two_smooth_model)
two_smooth_summary$p.table
two_smooth_summary$s.table
par(mfrow = c(2, 2))
plot(two_smooth_model, page = 1, all.terms = TRUE)
AIC(basic_model, two_term_model, two_smooth_model)
factor_interact <- gam(Sources ~ Season + s(SampleDepth, by = Season), data = isit, method = "REML")

summary(factor_interact)$s.table
par(mfrow = c(1, 2))
plot(factor_interact)
vis.gam(factor_interact, theta = 120, n.grid = 50, lwd = 0.4)
smooth_interact <- gam(Sources ~ Season + s(SampleDepth, RelativeDepth),
    data = isit, method = "REML")
summary(smooth_interact)$s.table
plot(smooth_interact, page = 1, scheme = 2)
vis.gam(smooth_interact, view = c("SampleDepth", "RelativeDepth"),
    theta = 50, n.grid = 50, lwd = 0.4)
AIC(two_smooth_model, factor_interact, smooth_interact)
smooth_interact <- gam(Sources ~ Season + s(SampleDepth, RelativeDepth),
    data = isit, method = "REML")

summary(smooth_interact)$s.table
par(mfrow=c(1, 3))
plot(gam(Sources ~ s(SampleDepth, k=4), bs = 'tp', method = "REML", data = isit2))
title(main="k=4")
plot(gam(Sources ~ s(SampleDepth, k=6), bs = 'tp', method = "REML", data = isit2))
title(main="k=6")
plot(gam(Sources ~ s(SampleDepth, k=10), bs = 'tp', method = "REML", data = isit2))
title(main="k=10")
k.check(smooth_interact)
smooth_interact_k60 <- gam(Sources ~ Season + s(SampleDepth,
    RelativeDepth, k = 60), data = isit, method = "REML")
k.check(smooth_interact_k60)
par(mfrow = c(2, 2))
gam.check(smooth_interact_k60)
'?'(family.mgcv)
smooth_interact_k100 <- gam(Sources ~ Season + s(SampleDepth, RelativeDepth, k = 60), family=gaussian(), data = isit, method = "REML")
summary(smooth_interact_k100)
