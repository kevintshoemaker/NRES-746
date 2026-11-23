#  Simulates the pika patch surveys used on the Week 15 page (dynamic occupancy).
#  Writes data/PikaPatches.csv: 100 talus patches, 8 summers, 3 visits per summer.

# Truth -------------------------------------------------------------------------------------
#   occupied in year 1 with probability psi (rises with elevation)
#   occupied -> empty with probability eps[t] (rises with that summer's temperature)
#   empty -> occupied with probability gam;  detected on a visit to an occupied patch with probability p

truth <- c(psi0 = 0.5, psi_elev = 1.0, eps0 = -1.6, eps_temp = 0.8, gam = 0.12, p = 0.6)
years <- 2018:2025
n_patch <- 100
n_visit <- 3
if (!exists("seed")) seed <- 3
if (!exists("write_csv")) write_csv <- TRUE

set.seed(seed)
elev_m <- round(runif(n_patch, 2300, 3300))
temp_z <- round(as.numeric(scale(rnorm(length(years)))), 2)     # summer temperature, SD units

z <- matrix(0L, n_patch, length(years))
z[, 1] <- rbinom(n_patch, 1, plogis(truth["psi0"] + truth["psi_elev"] * (elev_m - 2800) / 500))
for (t in 2:length(years)) {
  eps_t <- plogis(truth["eps0"] + truth["eps_temp"] * temp_z[t - 1])
  z[, t] <- rbinom(n_patch, 1, ifelse(z[, t - 1] == 1, 1 - eps_t, truth["gam"]))
}

pikas <- do.call(rbind, lapply(seq_along(years), function(t) {
  det <- matrix(rbinom(n_patch * n_visit, 1, truth["p"] * z[, t]), n_patch)
  data.frame(patch = sprintf("T%03d", seq_len(n_patch)), elev_m = elev_m, year = years[t], temp_z = temp_z[t],
             visit1 = det[, 1], visit2 = det[, 2], visit3 = det[, 3], occupied_true = z[, t])
}))
if (write_csv) write.csv(pikas, "data/PikaPatches.csv", row.names = FALSE)
