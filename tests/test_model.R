# Properties the model must have whatever the parameters. Base R only.
#
# Run from the repository root: Rscript tests/test_model.R

source(file.path("R", "tmdd_engine.R"))
source(file.path("R", "population.R"))

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}

P <- tmdd_individual(1)
tt <- seq(0, 56, by = 0.5)

# --- Quasi-steady-state solution ----------------------------------------------

set.seed(3)
ctot <- 10^runif(200, -4, 3); rtot <- 10^runif(200, -2, 2); kss <- 10^runif(200, -2, 1)
cf <- tmdd_free(ctot, rtot, kss)
cplx <- ctot - cf
check(all(cf >= 0 & cf <= ctot), "free drug is between 0 and total drug")
check(max(abs(cplx - (rtot - cplx) * cf / kss) / pmax(cplx, 1e-12)) < 1e-8,
      "free drug satisfies the binding equilibrium Cplx = Rfree * C / Kss")

# --- Baseline -------------------------------------------------------------------

zero <- tmdd_simulate(P, list(route = "iv", amt = 0, interval = 7, n_doses = 1), tt)
check(max(abs(zero$rtot - P$R0)) < 1e-10 && max(abs(zero$bio_rel - 1)) < 1e-10,
      "with no drug, target and biomarker stay at baseline")

# --- Linear limit -----------------------------------------------------------------

# Without target the model is a linear two-compartment model: compare
# against its closed form after an IV bolus.
lin <- P; lin$R0 <- 1e-9
dose_nmol <- 70 * 1e6 / lin$MW
r <- tmdd_simulate(lin, list(route = "iv", amt = 70, interval = 1e3, n_doses = 1), tt, rtol = 1e-10, atol = 1e-12)
k10 <- lin$CL / lin$V1; k12 <- lin$Q / lin$V1; k21 <- lin$Q / lin$V2
s <- k10 + k12 + k21
alpha <- (s + sqrt(s^2 - 4 * k10 * k21)) / 2; beta <- (s - sqrt(s^2 - 4 * k10 * k21)) / 2
exact <- dose_nmol / lin$V1 * ((alpha - k21) / (alpha - beta) * exp(-alpha * tt) +
                                 (k21 - beta) / (alpha - beta) * exp(-beta * tt))
check(max(abs(r$free_nM[, 1] - exact)) / max(exact) < 1e-7,
      "with no target the solver reproduces the two-compartment closed form")

# --- TMDD ---------------------------------------------------------------------------

auc <- function(dose) {
  r <- tmdd_simulate(P, list(route = "iv", amt = dose, interval = 1e3, n_doses = 1), seq(0, 400, by = 0.5))
  sum(diff(r$time) * (head(r$free_nM[, 1], -1) + tail(r$free_nM[, 1], -1)) / 2) / dose
}
check(auc(1) < 0.7 * auc(100), "dose-normalised exposure rises with dose (target-mediated clearance saturates)")

sc <- tmdd_simulate(P, list(route = "sc", amt = 700, interval = 1e3, n_doses = 1), seq(0, 1500, by = 1))
iv <- tmdd_simulate(P, list(route = "iv", amt = 700, interval = 1e3, n_doses = 1), seq(0, 1500, by = 1))
ro_sc <- max(sc$ro); ro_iv <- max(iv$ro)
check(ro_sc > 0.99 && ro_iv > 0.99, "a 10 mg/kg dose saturates the target by either route")

# --- Regimen handling ------------------------------------------------------------------

load <- tmdd_simulate(P, list(route = "sc", amt = 21, load = 70, interval = 14, n_doses = 4), tt)
noload <- tmdd_simulate(P, list(route = "sc", amt = 21, interval = 14, n_doses = 4), tt)
check(load$ro[tt == 13, 1] > noload$ro[tt == 13, 1], "a loading dose raises occupancy in the first interval")

# --- Population --------------------------------------------------------------------------

Pp <- tmdd_individual(500, c(40, 120), cv = list(CL = 30, V1 = 0, R0 = 0, KA = 0), seed = 2)
fit <- stats::lm(log(Pp$CL) ~ log(Pp$WT / 70))
check(abs(coef(fit)[2] - TYPICAL$WT_EXP_CL) < 0.1, "clearance scales with weight at the set exponent")

dr <- dose_ranging(c(0.03, 0.3, 3), Pp[1:50, ], list(route = "sc", amt = 1, interval = 14, n_doses = 4),
                   seq(0, 56, by = 1), 0.9)
check(all(diff(dr$ro_med) > 0) && all(diff(dr$pta) >= 0), "occupancy and target attainment rise with dose")

cat("PASS\n")
