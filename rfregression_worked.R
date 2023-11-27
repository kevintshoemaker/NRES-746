

library(randomForest)
library(rfUtilities)
my_data = read.csv("ET_data.csv")
my_data = na.omit(my_data)
# Initialize empty training and test sets
train_set <- data.frame()
test_set <- data.frame()

# Specify the number of rows to include in each set alternately
rows_per_set <- 3

# Create alternating sets
for (i in seq(1, nrow(my_data), by = rows_per_set * 2)) {
  test_indices <- i:(i + rows_per_set - 1)
  train_indices <- (i + rows_per_set):(i + rows_per_set * 2 - 1)

  test_set <- rbind(test_set, my_data[test_indices, , drop = FALSE])
  train_set <- rbind(train_set, my_data[train_indices, , drop = FALSE])
}
train_set <- na.omit(train_set)
set.seed(123)
rf <- randomForest(data = train_set ,x = train_set[,c(1,3:8)],y = train_set$ET,ntree = 600,mtry = 2,importance = TRUE,proximity = TRUE)
print(rf)
plot(rf)
set.seed(123)
tuneRF(y = train_set$ET,x = train_set[,c(1,3:8)],,mtryStart = 2,stepFactor = 3,trace = TRUE,plot = TRUE, ntreeTry = 600 )
cv <- rf.crossValidation(x= rf, xdata = my_data[,c(1,3:8)],ydata = my_data$ET,p = 0.2, n = 99, seed = 123)
mean(cv$fit.var.exp)
mean(cv$fit.mse)
mean(cv$y.rmse)
mean(cv$y.mbe)
mean(cv$y.mae)
mean(cv$D)
mean(cv$p.val)
varImpPlot(rf)

# Start of exercise ------------------------

library(randomForest)
library(rfUtilities)

my_data = read.csv("ET_data.csv")   # replace with your specific filepath if needed
my_data = na.omit(my_data)
# Initialize empty training and test sets
train_set <- data.frame()
test_set <- data.frame()

# Specify the number of rows to include in each set alternately
rows_per_set <- 3

# Create alternating sets
for (i in seq(1, nrow(my_data), by = rows_per_set * 2)) {
  test_indices <- i:(i + rows_per_set - 1)
  train_indices <- (i + rows_per_set):(i + rows_per_set * 2 - 1)

  test_set <- rbind(test_set, my_data[test_indices, , drop = FALSE])
  train_set <- rbind(train_set, my_data[train_indices, , drop = FALSE])
}

train_set <- na.omit(train_set)
set.seed(123)
rf <- randomForest(data = train_set ,x = train_set[,c(1,3:8)],y = train_set$ET,ntree = 600,mtry = 2,importance = TRUE,proximity = TRUE)
print(rf)

# Plotting OOB error rate

plot(rf)

# Tuning RF

set.seed(123)
tuneRF(y = train_set$ET,x = train_set[,c(1,3:8)],,mtryStart = 2,stepFactor = 3,trace = TRUE,plot = TRUE, ntreeTry = 600 )

# Cross validation

#- We trained the model with only train dataset.
#- we are using the complete data set to do cross validation.


cv <- rf.crossValidation(x= rf, xdata = train_set[,c(1,3:8)],ydata = train_set$ET,p = 0.2, n = 99, seed = 123)

#fit.var.exp
mean(cv$fit.var.exp)

#fit.mse
mean(cv$fit.mse)

#y.rmse
mean(cv$y.rmse)

#y.mbe
mean(cv$y.mbe)

#y.mae
mean(cv$y.mae)

# D (Kolmogorov-Smirnov distribution test)
mean(cv$D)

#p.val(p-value for the Kolmogorov-Smirnov distribution test)
mean(cv$p.val)

#Variable Importance Plot
varImpPlot(rf)

