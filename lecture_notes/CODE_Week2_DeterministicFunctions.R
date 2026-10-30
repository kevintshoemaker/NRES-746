#  NRES 746, Week 2: Deterministic functions and the bestiary --------------------
#  Code from Week2_DeterministicFunctions.html. Base R only.

# READING A DETERMINISTIC FUNCTION ==============================================

# 1. Big picture ----------------------------------------------------------------

frogs <- read.csv(
  "https://kevintshoemaker.github.io/NRES-746/data/ReedfrogSizepred.csv"
)
set.seed(1)
plot(jitter(frogs$TBL, amount = 0.5), frogs$Kill, pch = 21, bg = "grey70",
     cex = 1.3, xlab = "Tadpole size, TBL (mm)",
     ylab = "Number killed (of 10)")

# 2. WOD: R syntax and a seasonal function --------------------------------------

#  bug 1: unnamed arguments match by position
ricker_abx <- function(a, b, x) {
  a * x * exp(-b * x)
}
ricker_abx(10, 2, 0.5)                 # R reads a = 10, b = 2, x = 0.5
ricker_abx(a = 2, b = 0.5, x = 10)     # what was intended

#  bug 2: R has no constant e; use exp()
a <- 2
b <- 0.3
try(curve(a * x * e^(-b * x), from = 0, to = 15))
curve(a * x * exp(-b * x), from = 0, to = 15, lwd = 2,
      xlab = "x", ylab = "f(x)")

#  bug 3: pass the function, don't call it
ricker <- function(x, a, b = 0.2) a * x * exp(-b * x)
avals <- c(0.5, 1, 1.5)
try(sapply(avals, ricker(), x = 5))
sapply(avals, ricker, x = 5)           # each aval fills a, the first unnamed argument

#  y(t) = 5 + 5 sin(2 pi (t - 91) / 365): period 365, range 0-10, steepest at day 91
seasonal <- function(t) 5 + 5 * sin(2 * pi * (t - 91) / 365)
curve(seasonal(x), from = 1, to = 365, lwd = 2,
      xlab = "Day of year", ylab = "y(t)")
abline(v = 91, lty = 2, col = "grey50")
abline(h = 5, lty = 3, col = "grey70")

# 3. Writing and plotting functions in R ----------------------------------------
#  x first, then parameters with defaults (Bolker 3.3.1)

micmen <- function(x, a = 2, b = 1) a * x / (b + x)
micmen(3)                    # uses the defaults a = 2, b = 1
micmen(3, b = 5)             # override one parameter by name

curve(micmen(x), from = 0, to = 8, ylim = c(0, 9), lwd = 2,
      xlab = "x", ylab = "f(x)")
curve(micmen(x, b = 3), add = TRUE, lwd = 2, col = "steelblue")
curve(micmen(x, a = 8), add = TRUE, lwd = 2, col = "darkorange")
abline(h = c(2, 8), lty = 3, col = "grey60")
legend("topleft", c("a = 2, b = 1", "a = 2, b = 3", "a = 8, b = 1"),
       col = c("black", "steelblue", "darkorange"), lwd = 2, bty = "n")

#  'if' is not vectorized
hockey_if <- function(x, a = 1, s = 5) {
  if (x < s) a * x else a * s
}
hockey_if(3)                 # one value: fine
try(hockey_if(c(3, 7)))      # a vector: error

sapply(c(3, 7), hockey_if)
hockey <- function(x, a = 1, s = 5) ifelse(x < s, a * x, a * s)
hockey(c(3, 7))
#  [PAUSE] why does curve(hockey_if(x), from = 0, to = 10) fail?

# 4. Limits: throw out the lower-order terms ------------------------------------
#  x -> infinity: a x / (b + x) ~ a x / x = a;   x -> 0: ~ (a / b) x

a_mm <- 10
b_mm <- 5
curve(micmen(x, a = a_mm, b = b_mm), from = 0, to = 60, lwd = 2.5,
      ylim = c(0, 12), xlab = "x", ylab = "f(x)")
curve(a_mm / b_mm * x, add = TRUE, lty = 2, lwd = 2, col = "steelblue")
abline(h = a_mm, lty = 3, lwd = 2, col = "darkorange")
legend("bottomright", c("f(x)", "small-x line, ax/b", "asymptote, a"),
       col = c("black", "steelblue", "darkorange"), lty = c(1, 2, 3),
       lwd = 2, bty = "n")

micmen(c(b_mm, 9 * b_mm), a = a_mm, b = b_mm) / a_mm   # half at x = b, 90% at 9b
#  [PAUSE] widget: Holling III on log-log axes. why slope 2 at small x?

# 5. Derivatives and critical points --------------------------------------------
#  f'(x) = a (1 - b x) e^(-b x) = 0  ->  x* = 1/b, height a / (b e)

D(quote(a * x * exp(-b * x)), "x")

a_rk <- 0.9
b_rk <- 1 / 12
ricker_ab <- function(x, a, b) a * x * exp(-b * x)
optimize(ricker_ab, c(0, 40), a = a_rk, b = b_rk, maximum = TRUE)
c(peak = 1 / b_rk, height = a_rk / (b_rk * exp(1)))
#  [PAUSE] what does doubling b do to the peak's location and height?

#  discussion 2: quadratic seedling survival at 40% and 60% soil moisture
0.2 + 0.05 * c(40, 60) - 0.001 * c(40, 60)^2

# REPEATING A CALCULATION IN R ==================================================

# 2. Loops, vectorization, and sapply -------------------------------------------

xvec <- seq(0, 10, by = 0.5)

y_loop <- numeric(length(xvec))          # set aside storage first
for (i in seq_along(xvec)) {
  y_loop[i] <- micmen(xvec[i])
}

y_vec <- micmen(xvec)                    # vectorized
y_sap <- sapply(xvec, micmen)            # apply to each element

all.equal(y_loop, y_vec)
all.equal(y_sap, y_vec)

#  sweep over a parameter: one column per b value
bvals <- c(0.5, 1, 2, 4, 8)
curves <- sapply(bvals, function(b) micmen(xvec, a = 2, b = b))
dim(curves)
matplot(xvec, curves, type = "l", lty = 1, lwd = 2,
        col = hcl.colors(5, "Viridis"), xlab = "x", ylab = "f(x)")
legend("bottomright", paste("b =", bvals), col = hcl.colors(5, "Viridis"),
       lwd = 2, bty = "n")
#  [PAUSE] which dimension of 'curves' goes with x, and which with b?

# 3. Functions of two inputs ----------------------------------------------------
#  per-capita risk with a Holling type II response: alpha / (1 + alpha h N)

risk <- function(N, alpha, h = 0.84) alpha / (1 + alpha * h * N)
N_vals <- seq(1, 40, length.out = 40)
alpha_vals <- seq(0.01, 0.5, length.out = 40)
risk_grid <- outer(N_vals, alpha_vals, risk)
dim(risk_grid)

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
persp(N_vals, alpha_vals, risk_grid, theta = 40, phi = 25,
      xlab = "Density N", ylab = "Attack rate", zlab = "Risk",
      col = "lightskyblue1", border = "grey40")
image(N_vals, alpha_vals, risk_grid, col = hcl.colors(30, "YlOrRd", rev = TRUE),
      xlab = "Density N", ylab = "Attack rate")
contour(N_vals, alpha_vals, risk_grid, add = TRUE)
par(mfrow = c(1, 1))
#  [PAUSE] why does risk fall with density? where does density matter least?

# 4. Repeating a random calculation ---------------------------------------------

counts <- c(0, 2, 1, 4, 3, 0, 7, 2, 5, 1, 3, 0, 2, 6, 1, 4, 2, 3, 0, 9)
mean(counts)

set.seed(2)
sample(counts, 5)                          # one random set of 5 quadrats
means5 <- replicate(2000, mean(sample(counts, 5)))
quantile(means5, c(0.025, 0.975))
hist(means5, breaks = 30, col = "grey80", border = "white",
     xlab = "Mean count in 5 random quadrats", main = "")
abline(v = mean(counts), col = "darkorange", lwd = 2)
#  [PAUSE] what happens with all 20 quadrats drawn without replacement?

# THE BESTIARY OF FUNCTIONS =====================================================

# 2. WOD: limits and a derivative -----------------------------------------------
#  logistic f'(x) = b e^(-bx) / (1 + e^(-bx))^2, so f'(0) = b/4

logis_b <- function(x, b = 1) 1 / (1 + exp(-b * x))
D(quote(1 / (1 + exp(-b * x))), "x")
b <- 2
(logis_b(1e-6, b) - logis_b(-1e-6, b)) / 2e-6     # numerical slope at 0
b / 4

# 3. Critical points and inflection points --------------------------------------
#  Ricker inflection at 2/b: second derivative is zero there

d2 <- D(D(quote(a * x * exp(-b * x)), "x"), "x")
a <- 2
b <- 0.5
eval(d2, list(x = 2 / b, a = a, b = b))       # second derivative at 2/b
#  [PAUSE] widget: why does Michaelis-Menten have no critical or inflection point?

# 5. Where shapes come from -----------------------------------------------------
#  dn/dt = r n (1 - n/K), stepped forward with Euler's method

r <- 0.5
K <- 100
n0 <- 2
dt <- 0.1
times <- seq(0, 25, by = dt)
n <- numeric(length(times))
n[1] <- n0
for (i in seq_along(times)[-1]) {
  n[i] <- n[i - 1] + r * n[i - 1] * (1 - n[i - 1] / K) * dt
}
plot(times, n, type = "l", lwd = 4, col = "grey70",
     xlab = "Time", ylab = "Population size")
curve(K / (1 + (K / n0 - 1) * exp(-r * x)), add = TRUE, lwd = 2, lty = 2)
abline(h = K / 2, lty = 3, col = "grey50")
legend("bottomright", c("Euler loop", "logistic formula"),
       col = c("grey70", "black"), lwd = c(4, 2), lty = c(1, 2), bty = "n")

#  closest Michaelis-Menten to a monomolecular curve (least squares via optim)
mono_curve <- function(x) 5 * (1 - exp(-0.4 * x))
x_grid <- seq(0, 15, by = 0.01)
gap_ss <- function(p) sum((micmen(x_grid, a = p[1], b = p[2]) - mono_curve(x_grid))^2)
mm_fit <- optim(c(6, 3), gap_ss)$par        # closest Michaelis-Menten: a, b
mm_fit
max_gap <- max(abs(micmen(x_grid, mm_fit[1], mm_fit[2]) - mono_curve(x_grid)))
max_gap

curve(mono_curve(x), from = 0, to = 15, lwd = 2, ylim = c(0, 6.2),
      xlab = "x", ylab = "f(x)")
curve(micmen(x, mm_fit[1], mm_fit[2]), add = TRUE, lwd = 2, lty = 2,
      col = "steelblue")
abline(h = c(5, mm_fit[1]), lty = 3, col = c("black", "steelblue"))
legend("bottomright", c("monomolecular", "closest Michaelis-Menten"),
       col = c("black", "steelblue"), lty = c(1, 2), lwd = 2, bty = "n")
#  [PAUSE] both fit equally well up to x = 15. what is the maximum feeding rate?

# Bonus: Taylor approximations --------------------------------------------------
#  Ricker ~ a x - a b x^2 near 0; logistic ~ 1/2 + (b/4) x near its midpoint

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
curve(ricker_ab(x, a = 2, b = 0.5), from = 0, to = 4, lwd = 2.5,
      ylim = c(0, 2), xlab = "x", ylab = "f(x)", main = "Ricker near 0")
curve(2 * x - 2 * 0.5 * x^2, add = TRUE, lty = 2, lwd = 2, col = "steelblue")
curve(logis_b(x, b = 1), from = -5, to = 5, lwd = 2.5, ylim = c(0, 1),
      xlab = "x", ylab = "f(x)", main = "Logistic near its midpoint")
curve(0.5 + x / 4, add = TRUE, lty = 2, lwd = 2, col = "steelblue")
par(mfrow = c(1, 1))

# Bonus: changing parameters and transforming axes ------------------------------
#  shift right by c, up by d (Bolker 3.4.6)

mm_shift <- function(x, a = 4, b = 2, c = 3, d = 1) {
  ifelse(x < c, d, a * (x - c) / (b + (x - c)) + d)
}
curve(mm_shift(x), from = 0, to = 25, lwd = 2, ylim = c(0, 5.5),
      xlab = "x", ylab = "f(x)")
abline(h = c(1, 3, 5), v = c(3, 5), lty = 3, col = "grey60")
