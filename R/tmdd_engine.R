# ============================================================================
# mAb target-mediated drug disposition: quasi-steady-state (QSS) model.
#
# Two-compartment monoclonal antibody PK with IV or subcutaneous dosing,
# binding to a membrane target that is synthesised, degraded and, when
# bound, internalised with the drug. Binding is fast relative to everything
# else, so free drug and free target are in quasi-equilibrium with the
# complex (Gibiansky et al., J Pharmacokinet Pharmacodyn 2008;35:573):
#
#   C = 0.5 * ((Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4 Kss Ctot))
#
# A downstream biomarker is an indirect response to free target: its
# production rate scales as 1 - IMAX + IMAX * (Rfree / R0)^GAMMA, so full
# target suppression lowers it by at most IMAX.
#
# States (per subject):
#   DEPOT  subcutaneous depot (nmol)
#   CENT   total drug in the central compartment, free + bound (nmol)
#   PERI   free drug in the peripheral compartment (nmol)
#   RTOT   total target, free + bound (nM)
#   BIO    biomarker (units of its baseline)
#
# The same ODEs are in models/tmdd.cpp for mrgsolve. mrgsolve can't run in a
# browser, so the app integrates them here with an adaptive Dormand-Prince
# 5(4) scheme written in base R and vectorised across subjects.
# tests/test_engine_vs_mrgsolve.R checks the two against each other.
#
# Units: time in days, drug in nmol and nM, target in nM.
# ============================================================================

TMDD_STATES <- c("DEPOT", "CENT", "PERI", "RTOT", "BIO")

#' Free drug concentration from total drug and total target (QSS)
#'
#' Written in the form that avoids cancellation when free drug is a tiny
#' fraction of the total (low doses, excess target).
tmdd_free <- function(ctot, rtot, kss) {
  b <- ctot - rtot - kss
  disc <- sqrt(b * b + 4 * kss * ctot)
  ifelse(b >= 0, 0.5 * (b + disc), 2 * kss * ctot / (disc - b))
}

#' Right-hand side for every subject at once
#'
#' @param Y matrix, one row per subject, columns TMDD_STATES
#' @param P data frame of individual parameters, one row per subject
tmdd_rhs <- function(Y, P) {
  depot <- Y[, 1]; cent <- Y[, 2]; peri <- Y[, 3]; rtot <- Y[, 4]; bio <- Y[, 5]

  ctot <- cent / P$V1
  cf <- tmdd_free(ctot, rtot, P$KSS)
  cplx <- ctot - cf
  rfree <- rtot - cplx

  ksyn <- P$KDEG * P$R0

  cbind(
    -P$KA * depot,
    P$F * P$KA * depot - P$CL * cf - P$Q * cf + P$Q * peri / P$V2 - P$KINT * cplx * P$V1,
    P$Q * cf - P$Q * peri / P$V2,
    ksyn - P$KDEG * rtot - (P$KINT - P$KDEG) * cplx,
    P$KOUT * P$B0 * (1 - P$IMAX + P$IMAX * (pmax(rfree, 0) / P$R0)^P$GAMMA) - P$KOUT * bio
  )
}

tmdd_initial <- function(P) {
  cbind(0, 0, 0, P$R0, P$B0)
}

# Dormand-Prince 5(4) coefficients.
DP <- list(
  c = c(0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1),
  a = list(
    c(1 / 5),
    c(3 / 40, 9 / 40),
    c(44 / 45, -56 / 15, 32 / 9),
    c(19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
    c(9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
    c(35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
  ),
  e = c(71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)
)

#' Integrate from t0 to t1 with one shared adaptive step for all subjects
#'
#' The step is accepted only when every subject's scaled error is below 1,
#' so each subject is solved at least as accurately as on its own.
dopri5_segment <- function(Y, t0, t1, P, h, rtol, atol, max_steps = 1e5) {
  t <- t0
  k1 <- tmdd_rhs(Y, P)
  steps <- 0
  while (t < t1 - 1e-12) {
    h <- min(h, t1 - t)
    k <- vector("list", 7)
    k[[1]] <- k1
    for (i in 2:7) {
      a <- DP$a[[i - 1]]
      Yi <- Y
      for (j in seq_along(a)) if (a[j] != 0) Yi <- Yi + h * a[j] * k[[j]]
      if (i == 7) Ynew <- Yi
      k[[i]] <- tmdd_rhs(Yi, P)
    }
    err <- DP$e[1] * k[[1]]
    for (j in 3:7) err <- err + DP$e[j] * k[[j]]
    err <- h * err
    scale <- atol + rtol * pmax(abs(Y), abs(Ynew))
    enorm <- max(sqrt(rowMeans((err / scale)^2)))

    if (!is.finite(enorm)) {
      h <- h / 10
    } else if (enorm <= 1) {
      t <- t + h
      Y <- Ynew
      k1 <- k[[7]]                     # first-same-as-last
      h <- h * min(5, 0.9 * max(enorm, 1e-10)^(-0.2))
    } else {
      h <- h * max(0.2, 0.9 * enorm^(-0.2))
    }
    steps <- steps + 1
    if (steps > max_steps) stop("ODE solver exceeded the maximum number of steps")
  }
  list(Y = Y, h = h)
}

#' Dose times and per-subject amounts (nmol)
#'
#' @param regimen list(route = "iv" | "sc", amt (mg per subject, maintenance),
#'   load (mg per subject for the first dose, NA/0 = same as maintenance),
#'   interval (days), n_doses)
expand_tmdd_regimen <- function(regimen, P, t_end) {
  n <- nrow(P)
  times <- (seq_len(regimen$n_doses) - 1) * regimen$interval
  times <- times[times <= t_end + 1e-9]
  maint <- rep_len(regimen$amt, n)
  load <- if (is.null(regimen$load)) maint else rep_len(regimen$load, n)
  load <- ifelse(is.na(load) | load <= 0, maint, load)
  to_nmol <- 1e6 / P$MW
  amounts <- lapply(seq_along(times), function(i) (if (i == 1) load else maint) * to_nmol)
  list(times = times, amounts = amounts,
       cmt = if (identical(regimen$route, "sc")) 1 else 2)
}

#' Simulate a population
#'
#' @param P data frame of individual parameters (see tmdd_individual())
#' @param times output times (days)
#' @return list(time, and matrices [time x subject] of free drug (nM and
#'   ug/mL), receptor occupancy, total and free target, biomarker)
tmdd_simulate <- function(P, regimen, times, rtol = 1e-6, atol = 1e-9) {
  n <- nrow(P)
  times <- sort(unique(times))
  t_end <- max(times)
  reg <- expand_tmdd_regimen(regimen, P, t_end)
  grid <- sort(unique(round(c(0, times, reg$times), 10)))
  grid <- grid[grid <= t_end + 1e-9]
  out_idx <- match(round(times, 10), grid)

  store <- array(NA_real_, c(length(times), 5, n))
  Y <- tmdd_initial(P)
  h <- 0.01
  for (j in seq_along(grid)) {
    d <- which(abs(reg$times - grid[j]) < 1e-9)
    for (i in d) Y[, reg$cmt] <- Y[, reg$cmt] + reg$amounts[[i]]
    o <- which(out_idx == j)
    if (length(o)) for (oo in o) store[oo, , ] <- t(Y)
    if (j < length(grid)) {
      # Restart with a small step after a dose: the solution has a kink.
      if (length(d)) h <- min(h, 0.01)
      seg <- dopri5_segment(Y, grid[j], grid[j + 1], P, h, rtol, atol)
      Y <- seg$Y
      h <- seg$h
    }
  }
  tmdd_outputs(times, store, P)
}

#' Derived quantities from stored states [time x state x subject]
tmdd_outputs <- function(times, store, P) {
  nt <- length(times)
  rep_p <- function(x) matrix(rep(x, each = nt), nrow = nt)
  cent <- store[, 2, , drop = TRUE]; rtot <- store[, 4, , drop = TRUE]; bio <- store[, 5, , drop = TRUE]
  if (is.null(dim(cent))) {
    cent <- matrix(cent, nt); rtot <- matrix(rtot, nt); bio <- matrix(bio, nt)
  }
  ctot <- cent / rep_p(P$V1)
  cf <- tmdd_free(ctot, rtot, rep_p(P$KSS))
  cplx <- ctot - cf
  list(
    time = times,
    free_nM = cf,
    conc = cf * rep_p(P$MW) / 1e6,                  # ug/mL
    total_nM = ctot,
    ro = ifelse(rtot > 0, pmin(pmax(cplx / rtot, 0), 1), 0),
    rtot = rtot,
    rfree = rtot - cplx,
    rfree_rel = (rtot - cplx) / rep_p(P$R0),
    bio_rel = bio / rep_p(P$B0)
  )
}
