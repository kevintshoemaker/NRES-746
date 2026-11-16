#  Simulates the cheatgrass survey used on the Week 14 page (spatial regression).
#  Writes data/CheatgrassSurvey.csv (200 plots) and data/CheatgrassGrid.csv (0.5-km prediction grid).

# Landscape: 40 x 30 km, coordinates in km -------------------------------------------------

elev_fun <- function(x, y) 1500 + 18 * x + 120 * sin(y / 6) + 80 * cos(x / 7 + y / 9)    # m

# three fire perimeters (ellipses): center, radii, rotation
fires <- data.frame(cx = c(10, 27, 33), cy = c(20, 9, 24), rx = c(9, 10, 6), ry = c(6, 5.5, 5),
                    rot = c(0.3, -0.4, 0))
burned_fun <- function(x, y) {
  inside <- rep(FALSE, length(x))
  for (i in seq_len(nrow(fires))) {
    dx <- x - fires$cx[i]
    dy <- y - fires$cy[i]
    u <- dx * cos(fires$rot[i]) + dy * sin(fires$rot[i])
    v <- -dx * sin(fires$rot[i]) + dy * cos(fires$rot[i])
    inside <- inside | (u / fires$rx[i])^2 + (v / fires$ry[i])^2 <= 1
  }
  as.integer(inside)
}

# Truth ----------------------------------------------------------------------------------

truth <- c(b0 = 3.0, b_burn = 0.5, b_elev = -2.0, sigma_s = 0.7, ell = 6, sigma_e = 0.3)
if (!exists("seed")) seed <- 18

# Plots: 200 random locations, none in the roadless block (x 13-22 km, y 1-10 km) ----------

set.seed(seed)
xy <- matrix(numeric(0), 0, 2)
while (nrow(xy) < 200) {
  p <- c(runif(1, 0, 40), runif(1, 0, 30))
  if (!(p[1] > 13 & p[1] < 22 & p[2] > 1 & p[2] < 10)) xy <- rbind(xy, p)
}
plots <- data.frame(plot = sprintf("P%03d", 1:200), x_km = round(xy[, 1], 2), y_km = round(xy[, 2], 2))
plots$elev_m <- round(elev_fun(plots$x_km, plots$y_km))
plots$burned <- burned_fun(plots$x_km, plots$y_km)

# Grid: 0.5-km cells -------------------------------------------------------------------------

grid <- expand.grid(x_km = seq(0.25, 39.75, by = 0.5), y_km = seq(0.25, 29.75, by = 0.5))
grid$elev_m <- round(elev_fun(grid$x_km, grid$y_km))
grid$burned <- burned_fun(grid$x_km, grid$y_km)

# Spatial field at plots and grid cells together: exponential covariance -----------------

eps <- rnorm(200, 0, truth["sigma_e"])                    # plot-level noise (nugget)
all_xy <- rbind(as.matrix(plots[, c("x_km", "y_km")]), as.matrix(grid[, c("x_km", "y_km")]))
D <- as.matrix(dist(all_xy))
if (!exists("plots_only")) plots_only <- FALSE          # TRUE skips the grid (for quick checks)
keep <- if (plots_only) 1:200 else seq_len(nrow(D))
L <- chol(truth["sigma_s"]^2 * exp(-D[keep, keep] / truth["ell"]) + diag(1e-8, length(keep)))
s_all <- as.vector(t(L) %*% rnorm(length(keep)))
s_plot <- s_all[1:200]
if (!plots_only) grid$field_true <- round(s_all[-(1:200)], 3)

# Response: log cheatgrass biomass (g per square meter) ------------------------------------

elev_km <- (plots$elev_m - 1800) / 1000
plots$log_biomass <- round(truth["b0"] + truth["b_burn"] * plots$burned + truth["b_elev"] * elev_km +
                             s_plot + eps, 2)

# A count response on the same field (bonus section): seedlings in a 0.1 square-meter frame
plots$seedlings <- rpois(200, exp(0.8 + truth["b_burn"] * plots$burned + truth["b_elev"] * elev_km + s_plot))

if (!plots_only) {
  grid$log_biomass_true <- round(truth["b0"] + truth["b_burn"] * grid$burned +
                                   truth["b_elev"] * (grid$elev_m - 1800) / 1000 + grid$field_true, 3)
  write.csv(plots, "data/CheatgrassSurvey.csv", row.names = FALSE)
  write.csv(grid, "data/CheatgrassGrid.csv", row.names = FALSE)
}
