# Random Forest PET demo, worked ------------------------

# load packages ---------------------

library(randomForest)
library(rfUtilities)
library(faux)

# read in data ---------------------

df = read.csv("ET_data.csv")

df = df[order(df$DOY),]   # make sure data is sorted by day of year

# define response and predictors ----------------------------

names(df)

predictors <- c("DOY","wind_speed","Rs","Rs_reflected","Rn","RH","air_temp","random0","random10","random25","random50")
response <- "ET"


# generate random variables ------------------

df$random0 <- rnorm(nrow(df))
df$random10 <- faux::rnorm_pre(df[,response],r=0.1)
df$random25 <- faux::rnorm_pre(df[,response],r=0.25)
df$random50 <- faux::rnorm_pre(df[,response],r=0.5)


# separate test and training data ----------------

testndx <- seq(1,nrow(df),5)   # rows to withhold
trainndx <- setdiff(1:nrow(df),testndx)

testdat <- df[testndx,]
traindat <- df[trainndx,]


# tune random forest analysis ---------------

rf.tune <- tuneRF(y = traindat[,response],x = traindat[,predictors],ntreeTry=500,
       mtryStart = 11,stepFactor = 2,trace = TRUE,plot = TRUE )    # find optimal mtry parameter

# run random forest analysis -----------------

rf <- randomForest(x = traindat[,predictors],y = traindat[,response],
                   ntree = 500,mtry = 4,importance = TRUE,proximity = TRUE)

rf

plot(rf)  # Plotting OOB error rate


# Cross validation for fitted model -------------------------

cv <- rf.crossValidation(x= rf, xdata = traindat[,predictors], ydata = traindat[,response], p = 0.1, n = 50)
summary(cv)

plot(cv, stat = "mse")

cv.v <- rf.crossValidation(x= rf, xdata = testdat[,predictors], ydata = testdat[,response], p = 0.5, n = 25)
plot(cv.v, stat = "mse")

summary(cv.v)   # note higher MSE for validation set...


#Variable Importance Plot -----------------------
varImpPlot(rf)



