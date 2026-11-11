# NRES 746, Bonus: Machine learning with random forests (code from the web page)
# Packages: ranger, party

# STATISTICS AND MACHINE LEARNING ==============================================
# (Concepts and a widget only; no code.)
#  [PAUSE] (widget) Add trees: why does the error on new data level off instead of reaching 0?

# EXAMPLE: WHO SURVIVED THE TITANIC? ==========================================

# 1. Data and a single tree -----------------------------------------------------

titanic <- read.csv(
  "https://kevintshoemaker.github.io/NRES-746/data/titanic.csv"
)
features <- c("Pclass", "Sex", "Age", "SibSp", "Parch", "Fare")
titanic <- na.omit(titanic[, c("Survived", features)])   # 714 passengers with known age
titanic$Survived <- factor(titanic$Survived, levels = 0:1, labels = c("died", "survived"))
titanic$Sex <- factor(titanic$Sex)
titanic$Pclass <- factor(titanic$Pclass)       # ticket class 1, 2, 3
str(titanic)

library(party)
tree <- ctree(Survived ~ ., data = titanic,
              controls = ctree_control(mincriterion = 0.85, maxdepth = 3))
plot(tree)

# 2. Fitting a forest -----------------------------------------------------------

library(ranger)
rf <- ranger(Survived ~ ., data = titanic, num.trees = 500, probability = TRUE,
             importance = "permutation", seed = 746)
rf

# 3. Variable importance ---------------------------------------------------------

imp <- sort(importance(rf))
par(mar = c(4, 7, 0.5, 0.5))
barplot(imp, horiz = TRUE, las = 1, col = "steelblue",
        xlab = "Permutation importance (increase in prediction error)")

# 4. Partial dependence ----------------------------------------------------------

# mean predicted survival when 'feature' is set to each value in 'values'
partial_dependence <- function(feature, values) {
  sapply(seq_along(values), function(i) {
    nd <- titanic
    nd[[feature]] <- values[i]
    mean(predict(rf, nd)$predictions[, "survived"])
  })
}

par(mfrow = c(2, 3), mar = c(4, 4, 1, 0.5))
for (f in names(sort(importance(rf), decreasing = TRUE))) {
  x <- titanic[[f]]
  values <- if (is.factor(x)) factor(levels(x), levels = levels(x)) else
    seq(min(x), quantile(x, 0.98), length.out = 30)    # skip the extreme fares
  pd <- partial_dependence(f, values)
  if (is.factor(x)) {
    barplot(pd, names.arg = levels(x), ylim = c(0, 1), col = "steelblue",
            xlab = f, ylab = "Predicted survival")
  } else {
    plot(values, pd, type = "l", lwd = 2, ylim = c(0, 1), xlab = f,
         ylab = "Predicted survival")
    rug(jitter(x[x <= max(values)]))
  }
}
par(mfrow = c(1, 1))

# 5. Interactions ----------------------------------------------------------------

# two-way partial dependence on a grid, and its departure from an additive fit
grid_values <- function(x) {
  if (is.factor(x)) factor(levels(x), levels = levels(x)) else
    seq(min(x), quantile(x, 0.98), length.out = 8)
}
pd_2way <- function(f1, f2) {
  grid <- expand.grid(v1 = grid_values(titanic[[f1]]), v2 = grid_values(titanic[[f2]]))
  nd <- titanic[rep(seq_len(nrow(titanic)), times = nrow(grid)), ]   # every passenger at every grid point
  nd[[f1]] <- rep(grid$v1, each = nrow(titanic))
  nd[[f2]] <- rep(grid$v2, each = nrow(titanic))
  grid$pd <- tapply(predict(rf, nd)$predictions[, "survived"],
                    rep(seq_len(nrow(grid)), each = nrow(titanic)), mean)
  grid
}
pairs <- as.data.frame(t(combn(features, 2)))
names(pairs) <- c("f1", "f2")
pairs$strength <- sapply(seq_len(nrow(pairs)), function(i) {
  g <- pd_2way(pairs$f1[i], pairs$f2[i])
  sqrt(mean(resid(lm(pd ~ factor(v1) + factor(v2), data = g))^2))   # RMS departure from additive
})
pairs <- pairs[order(pairs$strength, decreasing = TRUE), ]
head(pairs, 5)

# two-way partial dependence: x axis = the feature with more grid values, one line per value of the other
plot_2way <- function(f1, f2) {
  g <- pd_2way(f1, f2)
  if (nlevels(factor(g$v1)) < nlevels(factor(g$v2))) {
    g <- g[, c("v2", "v1", "pd")]
    names(g) <- c("v1", "v2", "pd")
    tmp <- f1; f1 <- f2; f2 <- tmp
  }
  x <- as.numeric(g$v1)
  groups <- split(seq_len(nrow(g)), factor(g$v2))
  plot(NA, xlim = range(x), ylim = c(0, 1), xlab = f1, ylab = "Predicted survival",
       xaxt = if (is.factor(g$v1)) "n" else "s")
  if (is.factor(g$v1)) axis(1, at = seq_along(levels(g$v1)), labels = levels(g$v1))
  for (k in seq_along(groups)) lines(x[groups[[k]]], g$pd[groups[[k]]], lwd = 2, lty = k)
  legend("topright", paste(f2, "=", names(groups)), lwd = 2, lty = seq_along(groups),
         bty = "n", cex = 0.8)
}

par(mfrow = c(1, 2), mar = c(4, 4, 1, 0.5))
plot_2way(pairs$f1[1], pairs$f2[1])
plot_2way(pairs$f1[2], pairs$f2[2])
par(mfrow = c(1, 1))

# 6. Performance -----------------------------------------------------------------

y <- as.numeric(titanic$Survived == "survived")
auc <- function(p, y) {
  r <- rank(p)
  n1 <- sum(y)
  n0 <- sum(1 - y)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
dev <- function(p) {
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)       # avoid log(0) for predictions of exactly 0 or 1
  -2 * sum(dbinom(y, 1, p, log = TRUE))
}
mcfadden <- function(p) 1 - dev(p) / dev(rep(mean(y), length(y)))

set.seed(746)
folds <- sample(rep(1:10, length.out = nrow(titanic)))
p_cv_rf <- p_cv_glm <- numeric(nrow(titanic))
for (k in 1:10) {
  train <- titanic[folds != k, ]
  test <- titanic[folds == k, ]
  m_rf <- ranger(Survived ~ ., data = train, num.trees = 500, probability = TRUE, seed = k)
  p_cv_rf[folds == k] <- predict(m_rf, test)$predictions[, "survived"]
  m_glm <- glm(Survived ~ ., data = train, family = binomial)
  p_cv_glm[folds == k] <- predict(m_glm, test, type = "response")
}
m_glm <- glm(Survived ~ ., data = titanic, family = binomial)

preds <- list(forest_training_data = predict(rf, titanic)$predictions[, "survived"],
              forest_out_of_bag = rf$predictions[, "survived"],
              forest_10fold_cv = p_cv_rf,
              logistic_training_data = fitted(m_glm),
              logistic_10fold_cv = p_cv_glm)
perf <- t(sapply(preds, function(p) c(AUC = auc(p, y), McFadden = mcfadden(p))))
round(perf, 3)

#  [PAUSE] Which model would you use to predict survival for a new passenger list?
#          Which would you use to explain why some passengers survived?
