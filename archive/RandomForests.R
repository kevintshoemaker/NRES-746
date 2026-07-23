############################################
#### Random Forests                     ####
#### William Richardson, Wade Lieurance ####
#### 2019-11-27                         ####
############################################

### Thinking about Partitioning ###
## Simple Example: Quantiles ##
library(dplyr)
library(tibble)

# lets look at example data mtcars
data(mtcars)
mtcars <- as_tibble(rownames_to_column(mtcars, var = "model"))

# a function to split the data into quatiles and determine the means and SSE
quantile.range <- function(data, return.var, split.var, split.n = 10){
  q <- quantile(data[[split.var]], probs = seq(0, 1, 1/split.n))
  y.store <- vector(mode="numeric", length=0)
  min.store <- vector(mode="numeric", length=0)
  max.store <- vector(mode="numeric", length=0)
  sse.store <- vector(mode="numeric", length=0)
  for (i in 1:length(q - 1)){
    min <- q[i]
    max <- q[i+1]
    restrict <- filter(data, .data[[split.var]] >= min & .data[[split.var]] < max)
    restrict.col <- restrict[[return.var]]
    y.pred <- mean(restrict.col)
    y.store <-  c(y.store, y.pred)
    min.store <- c(min.store, min)
    max.store <- c(max.store, max)
    sse <- 0
    for (j in  restrict.col){
      res.err <- (j - y.pred)^2 
      sse <-  sse + res.err
    }
    sse.store <- c(sse.store, sse)
  }
  df <- as.data.frame(cbind(y.store, min.store, max.store, sse.store))
  names(df) <- c(return.var, "min", "max", "sse")
  # df <- filter(df, !is.na(.data[[return.var]]))
  return(df)
}

# let look at the type of SSE we get with a simple linear model
plot(mpg ~ disp, data = mtcars, main = "Linear regression")
lm1 <- lm(mpg ~ disp, data = mtcars)
abline(lm1, col = "yellow", lwd = 1.5)
SSE.lm <- with(summary(lm1), df[2] * sigma^2)

# what if we divide displacement in quartiles?
means.quartile <- quantile.range(mtcars, "mpg", "disp", split.n = 4)
SSE.quartile <- sum(means.quartile$sse) 
plot(mpg ~ disp, data = mtcars, main = "Quartile means")
segments(means.quartile$min, means.quartile$mpg, means.quartile$max, 
         means.quartile$mpg, col = "blue", lwd = 2)
abline(v = min(means.quartile$min), col = "gray")
abline(v = means.quartile$max, col = "gray")

# what if we divide displacement in deciles?
means.decile <- quantile.range(mtcars, "mpg", "disp", split.n = 10)
SSE.decile <- sum(means.decile$sse) 
plot(mpg ~ disp, data = mtcars, main = "Decile means")
segments(means.decile$min, means.decile$mpg, means.decile$max, means.decile$mpg, 
         col = "red", lwd = 2)
abline(v = min(means.decile$min), col = "gray")
abline(v = means.decile$max, col = "gray")

# what if we go overboard an divide displacement into n = unique(disp)
means.n <- quantile.range(mtcars, "mpg", "disp", 
                          split.n = length(unique(mtcars$disp)))
SSE.n <- sum(means.n$sse) 
plot(mpg ~ disp, data = mtcars, main = "n means")
segments(means.n$min, means.n$mpg, means.n$max, means.n$mpg, 
         col = "violet", lwd = 2)
abline(v = min(means.n$min), col = "gray")
abline(v = means.n$max, col = "gray")


### Decision Trees - A partial solution... ###
## Binary Splitting with Continuous Response (Regression Trees) ##

# binary splitting to reduce SSE example.
find_split <- function(data, y.var, x.var){
  require(foreach)
  require(dplyr)
  data <- arrange(data, .data[[x.var]])  # puts our input in order of x.var 
  se <- function(y, y.bar) (y-y.bar)^2
  sse.breaks <- foreach (i = 1:(nrow(data)-1), .combine = bind_rows) %do% {
    x.mean <- mean(c(data[[x.var]][i], data[[x.var]][i+1]))
    left.data <- head(data[[y.var]], i)
    right.data <- tail(data[[y.var]], -i)
    left.mean <- mean(left.data)
    right.mean <- mean(right.data)
    se.left <- sapply(left.data, se, y.bar = mean(left.data))
    se.right <- sapply(right.data, se, y.bar = mean(right.data))
    sse = sum(c(se.left, se.right))
    out <- tibble(!!x.var := x.mean, sse = sse)
    out
  }
  return(sse.breaks)
}

# plot the first partitioning into regions R1 and R2 to minimize SSE for each 
# region
SSEs <- find_split(mtcars, "mpg", "disp")
plot(sse ~ disp, data = SSEs, ylim = c(min(SSEs$sse)-200, max(SSEs$sse)+200),
     main = "SSE by disp break value")
abline(h = min(SSEs$sse), col = "red", lty = 3)
x.break <- filter(SSEs, sse == min(SSEs$sse))
abline(v = x.break$disp, col = "red")

#lets see where this breaks off
plot(mpg ~ disp, data = mtcars, main = "Binary breaks")
abline(v = x.break$disp, col = "red")
axis(1, at=x.break$disp, labels=round(x.break$disp, 0))


# lets develop a recursive function to get n splits.
recursive.split <- function(data, y.var, x.var, n.splits, r = 1){
  SSEs <- find_split(data, y.var, x.var)
  x.break <- filter(SSEs, sse == min(SSEs$sse)) %>% 
    mutate(n = r)
  left.df <- filter(data, .data[[x.var]] <= x.break[[x.var]])
  right.df <- filter(data, .data[[x.var]] > x.break[[x.var]])
  if ((n.splits - 1) > 0){
    bl <- recursive.split(left.df, y.var, x.var, n.splits -1, r +1)
    br <- recursive.split(right.df, y.var, x.var, n.splits -1, r +1)
    x.break <- rbind(x.break, bl, br)
  }
  return(arrange(x.break, n, .data[[x.var]]))
}
R <- recursive.split(data = mtcars, y.var = "mpg", x.var = "disp", n.splits = 3)
plot(mpg ~ disp, data = mtcars, main = "Binary breaks", xaxt="n" )
pal <- colorRampPalette(c("red", "blue"))(max(R$n))
R <- R %>% mutate(col = pal[n])
abline(v = R$disp, lwd = max(R$n)/R$n, col = R$col)
axis(1, R$disp, labels=round(R$disp, 0), las = 2)


# lets produce a simple tree using the rpart library and the above examples
library(rpart)
library(rpart.plot)
tree <- rpart(formula = mpg~disp, data = mtcars, 
              control = rpart.control(maxdepth = 3, minbucket = 1, cp = 0))
prp(tree, faclen = 0, cex = 0.8, extra = 1)


## Binary Splitting with Multiple Variables ##
# what about a multivariate example?
tree <- rpart(formula = mpg~disp+wt, data = mtcars, 
              control = rpart.control(maxdepth = 3, minbucket = 1, cp = 0))

# lets use another package rattle to get a fancier looking output
library(rattle)
library(RColorBrewer)
fancyRpartPlot(tree, caption = NULL)
# splits <- data.frame(var = row.names(tree$splits), tree$splits)
# rownames(splits) <- NULL
wts <- c(1.7, 2.3, 3.3)
disps <- c(79, 267, 450)

pal <- colorRampPalette(c("yellow", "blue"))(length(unique(mtcars$mpg)))
colors <- tibble(color = pal, mpg = sort(unique(mtcars$mpg)))
d <- mtcars %>% left_join(colors, by = c("mpg" = "mpg"))
plot(wt~disp, data = d, pch = 19, col = color, 
     main = "Mpg by Displacement and Weight")
abline(h = 2.3)
segments(79, 0, 79, 2.3)
abline(h = 1.7)
segments(267, 2.3, 267, 10)
segments(450, 2.3, 450, 10)
segments(0, 3.3, 267, 3.3)

tree.m <- rpart(formula = mpg~disp+wt+cyl+drat+vs+am+gear+carb, data = mtcars, 
              control = rpart.control(maxdepth = 3, minbucket = 1, cp = 0))
fancyRpartPlot(tree.m, caption = NULL)

## Classification Trees ##
mtcars.class <- mutate(mtcars,
                mpg = case_when(mpg < 16.7 ~ "low",
                                mpg >= 16.7 & mpg < 21.4 ~ "medium",
                                mpg >= 21.4 ~ "high"))

# gini impurity
set <- c("horse", "horse", "cart", "cart", "cart")
set
p.horse <- length(set[set == "horse"])/length(set) # 2/5
p.cart <-  length(set[set == "cart"])/length(set)  # 3/5
I <- (p.horse * (1-p.horse))+ (p.cart * (1-p.cart))
I  # 0.48

# a bad split
set.1 <- c("horse", "cart")
set.2 <- c("horse", "cart", "cart")

I.1 <- (0.5*(1-0.5)) + (0.5*(1-0.5))  # 0.5
I.2 <- (1/3 * (1 - 1/3)) + (2/3 * (1 - 2/3))  # 0.44
I.new <- (I.1 * length(set.1)/length(set)) + (I.2 * length(set.2)/length(set))
I.new
I.gain <- I - I.new
I.gain

# classification tree pruning
tree.class <- rpart(formula = mpg~disp+wt+hp+vs+gear, data = mtcars.class, 
              control = rpart.control(maxdepth = 3, minbucket = 1, cp = 0))
fancyRpartPlot(tree.class, caption = NULL)

printcp(tree)
tree.p <- prune(tree, cp = 0.013)
fancyRpartPlot(tree.p, caption = NULL)


### Regression Tree Example ###
library(rsample) 
library(dplyr)
library(rpart)
library(rpart.plot)
library(ModelMetrics)
library(AmesHousing)

## Split Data ##

set.seed(123)
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

## Basic regression tree model ##

m1 <- rpart(
  formula = Sale_Price ~ .,
  data    = ames_train,
  method  = "anova"
)

m1
rpart.plot(m1)
plotcp(m1)

# Model when cp=0
m2 <- rpart(
  formula = Sale_Price ~ .,
  data    = ames_train,
  method  = "anova", 
  control = list(cp = 0, xval = 10)
)

plotcp(m2)
abline(v = 12, lty = "dashed")

## Tuning ##

m3 <- rpart(
  formula = Sale_Price ~ .,
  data    = ames_train,
  method  = "anova", 
  control = list(minsplit = 10, maxdepth = 12, xval = 10)
)

m3$cptable

# Create grid and do for loop
hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(8, 15, 1)
)

head(hyper_grid)
nrow(hyper_grid)
models <- list()

for (i in 1:nrow(hyper_grid)) {
  
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  
  models[[i]] <- rpart(
    formula = Sale_Price ~ .,
    data    = ames_train,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}

#get top models
get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}

get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}

hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models, get_cp),
    error = purrr::map_dbl(models, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)

## Predict using optimal model

optimal_tree <- rpart(
  formula = Sale_Price ~ .,
  data    = ames_train,
  method  = "anova",
  control = list(minsplit = 11, maxdepth = 8, cp = 0.01)
)

pred <- predict(optimal_tree, newdata = ames_test)
rmse(pred = pred, actual = ames_test$Sale_Price)


### Bagging Example ###
library(rsample) 
library(dplyr)
library(ipred)       
library(caret)
library(ModelMetrics)
library(AmesHousing)

## Split data
set.seed(123)

ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

## Make bagged regression model
set.seed(123)
bagged_m1 <- bagging(
  formula = Sale_Price ~ .,
  data    = ames_train,
  coob    = TRUE
)
bagged_m1

## Mess with number of trees
ntree <- 10:70
rmse <- vector(mode = "numeric", length = length(ntree))

for (i in seq_along(ntree)) {
  set.seed(123)
  model <- bagging(
    formula = Sale_Price ~ .,
    data    = ames_train,
    coob    = TRUE,
    nbagg   = ntree[i]
  )
  rmse[i] <- model$err
}

plot(ntree, rmse, type = 'l', lwd = 2)

## Use caret
ctrl <- trainControl(method = "cv",  number = 10) 

bagged_cv <- train(
  Sale_Price ~ .,
  data = ames_train,
  method = "treebag",
  trControl = ctrl,
  importance = TRUE
)

bagged_cv
plot(varImp(bagged_cv), 20)

## Predict using training data
pred <- predict(bagged_cv, ames_test)
rmse(pred, ames_test$Sale_Price)



### Random Forest Example ###
library(rsample)      
library(randomForest)
library(ranger)      
library(caret)        
library(h2o)
library(AmesHousing)

## Split Data
set.seed(123)
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

## Basic Model
set.seed(123)
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)
m1

## Plot it out
plot(m1)

## Find min mse
which.min(m1$mse)
sqrt(m1$mse[which.min(m1$mse)])

## Tune your model
features <- setdiff(names(ames_train), "Sale_Price")
set.seed(123)
m2 <- tuneRF(
  x          = ames_train[features],
  y          = ames_train$Sale_Price,
  ntreeTry   = 500,
  mtryStart  = 5,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time progress 
)

## Check system times
system.time(
  ames_randomForest <- randomForest(
    formula = Sale_Price ~ ., 
    data    = ames_train, 
    ntree   = 500,
    mtry    = floor(length(features) / 3)
  )
)

system.time(
  ames_ranger <- ranger(
    formula   = Sale_Price ~ ., 
    data      = ames_train, 
    num.trees = 500,
    mtry      = floor(length(features) / 3)
  )
)

## Create grid of parameters

hyper_grid <- expand.grid(
  mtry       = seq(20, 30, by = 2),
  node_size  = seq(3, 9, by = 2),
  sampe_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

nrow(hyper_grid)

## Fill it in
for(i in 1:nrow(hyper_grid)) {
  
  model <- ranger(
    formula         = Sale_Price ~ ., 
    data            = ames_train, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sampe_size[i],
    seed            = 123
  )
  
  hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

## Let's Use the best Model
OOB_RMSE <- vector(mode = "numeric", length = 100)

for(i in seq_along(OOB_RMSE)) {
  optimal_ranger <- ranger(
    formula         = Sale_Price ~ ., 
    data            = ames_train, 
    num.trees       = 500,
    mtry            = 28,
    min.node.size   = 3,
    sample.fraction = .8,
    importance      = 'impurity'
  )
  OOB_RMSE[i] <- sqrt(optimal_ranger$prediction.error)
}

hist(OOB_RMSE, breaks = 20)

## Measure impurity
optimal_ranger <- ranger(
  formula         = Sale_Price ~ ., 
  data            = ames_train, 
  num.trees       = 500,
  mtry            = 28,
  min.node.size   = 3,
  sample.fraction = .8,
  importance      = 'impurity'
)

plot(optimal_ranger$variable.importance)
which.max(optimal_ranger$variable.importance)
which.min(optimal_ranger$variable.importance)

## Run h2o package
h2o.no_progress()
h2o.init(max_mem_size = "5g")

y <- "Sale_Price"
x <- setdiff(names(ames_train), y)
train.h2o <- as.h2o(ames_train)

hyper_grid.h2o <- list(
  ntrees      = seq(200, 500, by = 100),
  mtries      = seq(20, 30, by = 2),
  sample_rate = c(.55, .632, .70, .80)
)

### WARNING: CPU INTENSIVE ###
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = list(strategy = "Cartesian")
)

grid_perf <- h2o.getGrid(
  grid_id = "rf_grid", 
  sort_by = "mse", 
  decreasing = FALSE
)
print(grid_perf)

## Faster H2o method
hyper_grid.h2o <- list(
  ntrees      = seq(200, 500, by = 150),
  mtries      = seq(15, 35, by = 10),
  max_depth   = seq(20, 40, by = 5),
  min_rows    = seq(1, 5, by = 2),
  nbins       = seq(10, 30, by = 5),
  sample_rate = c(.55, .632, .75)
)

search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid2",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

grid_perf2 <- h2o.getGrid(
  grid_id = "rf_grid2", 
  sort_by = "mse", 
  decreasing = FALSE
)
print(grid_perf2)

## Calculate best models's RMSE
best_model_id <- grid_perf2@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

ames_test.h2o <- as.h2o(ames_test)
best_model_perf <- h2o.performance(model = best_model, newdata = ames_test.h2o)

h2o.mse(best_model_perf) %>% sqrt()

## Check all models
pred_randomForest <- predict(ames_randomForest, ames_test)
head(pred_randomForest)

pred_ranger <- predict(ames_ranger, ames_test)
head(pred_ranger$predictions)

pred_h2o <- predict(best_model, ames_test.h2o)
head(pred_h2o)
