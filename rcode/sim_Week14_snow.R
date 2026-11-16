#  Simulates the snow course data for the Week 14 take-home WOD (Neighbors).
#  Writes data/SnowCourses.csv: 60 stations in a 60 x 40 km mountain range.
#  Truth: SWE = 40 + 45 (elev_km - 2.4) + field + noise; field SD 8 cm, exponential, ell = 10 km; noise SD 3 cm.

if (!exists("seed")) seed <- 65
set.seed(seed)
snow <- data.frame(station = sprintf("S%02d", 1:60), x_km = round(runif(60, 0, 60), 1), y_km = round(runif(60, 0, 40), 1))
snow$elev_m <- round(2100 + 8 * snow$x_km + 150 * sin(snow$y_km / 9) + runif(60, -250, 250), -1)
D <- as.matrix(dist(snow[, c("x_km", "y_km")]))
field <- as.vector(t(chol(8^2 * exp(-D / 10))) %*% rnorm(60))
snow$swe_cm <- round(40 + 45 * (snow$elev_m - 2400) / 1000 + field + rnorm(60, 0, 3), 1)
if (!exists("no_write")) write.csv(snow, "data/SnowCourses.csv", row.names = FALSE)
