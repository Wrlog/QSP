# Check the base-R Dormand-Prince engine against mrgsolve running
# models/tmdd.cpp. If they agree, the browser build computes what the
# mrgsolve model computes.
#
# Run from the repository root: Rscript tests/test_engine_vs_mrgsolve.R

source(file.path("R", "tmdd_engine.R"))
source(file.path("R", "population.R"))
source(file.path("reference", "mrgsolve_engine.R"))
if (!mrgsolve_ready()) stop("mrgsolve is not installed")

# Off-grid output times, so no comparison lands exactly on a dose time.
times <- seq(0, 112, by = 0.5) + 0.0137
rel <- function(x, y) max(abs(x - y)) / max(abs(y))

# The default target, one that accumulates when bound (KINT < KDEG) and a
# high-affinity, high-capacity one, so every branch of the QSS solution
# and both directions of target change are exercised.
targets <- list(
  default      = TYPICAL,
  accumulating = modifyList(TYPICAL, list(KINT = 0.05, KDEG = 0.5, R0 = 0.5)),
  high_capacity = modifyList(TYPICAL, list(KSS = 0.05, R0 = 20, KINT = 3, GAMMA = 2))
)

worst <- 0
cases <- 0
for (nm in names(targets)) {
  P <- tmdd_individual(5, c(45, 110), targets[[nm]],
                       cv = list(CL = 30, V1 = 20, R0 = 30, KA = 20), seed = 11)
  for (route in c("iv", "sc")) {
    for (dose in c(0.03, 0.3, 3)) {
      reg <- list(route = route, amt = dose * P$WT, load = 2 * dose * P$WT,
                  interval = 14, n_doses = 6)
      a <- tmdd_simulate(P, reg, times)
      b <- tmdd_simulate_mrgsolve(P, reg, times)
      d <- c(conc = rel(a$conc, b$conc), ro = max(abs(a$ro - b$ro)),
             rtot = rel(a$rtot, b$rtot), bio = rel(a$bio_rel, b$bio_rel))
      cat(sprintf("  %-13s %-2s %5.2f mg/kg  conc %.1e  RO %.1e  Rtot %.1e  biomarker %.1e\n",
                  nm, route, dose, d[1], d[2], d[3], d[4]))
      worst <- max(worst, d)
      cases <- cases + 1
    }
  }
}

cat(sprintf("cases compared: %d\n", cases))
cat(sprintf("worst difference vs mrgsolve: %.3e\n", worst))
# The app's solver runs at rtol 1e-6; this is well inside what a plot or a
# summary statistic can show.
if (!is.finite(worst) || worst > 1e-4) {
  stop(sprintf("Dormand-Prince engine disagrees with mrgsolve (%.3e)", worst))
}
cat("PASS\n")
