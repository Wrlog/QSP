# ============================================================================
# Vaccine settings for the mRNA vaccine app: fitted parameter sets and
# standard schedules. Derivative of the COSBI "Multiscale QSP model for mRNA
# vaccines" (COSBI-SSLA, non-commercial; see ../LICENSE). Modified 2026-09-28.
# ============================================================================

AB_NAMES <- c("CC_N", "SigmaMB", "RhoAB_N", "RhoAB_M", "LambdaSP", "BetaA", "expAg", "SigmaNT", "SigmaMT",
              "RhoAT", "BetaLP", "BetaSP", "mRNA_max", "ksat", "delayB")

# opt_param vectors fitted to antibody data (BNT162b2_fit.mat,
# BNT162b2_fit_over60.mat, mRNA-1273_fit.mat)
VACCINES <- list(
  bnt = list(label = "BNT162b2 (Pfizer-BioNTech)", dose = 30, interval = 21,
             fit = c(6.6178599775e+01, 8.6211997917e+00, 6.0612000850e+00, 6.3442991786e+00, 4.0482099838e+01,
                     9.1799901351e-02, 3.3713258652e+05, 2.9465880005e+02, 9.3584379878e+02, 4.9808998282e+00,
                     1.4300002628e-02, 1.3830001542e-01, 3.1416002938e+00, 5.5832999644e+01, 7.9505908077e-07)),
  bnt60 = list(label = "BNT162b2, adults over 60", dose = 30, interval = 21,
               fit = c(6.6996197660e+01, 8.6646607499e+00, 6.0705847520e+00, 6.2592698731e+00, 4.0230694929e+01,
                       9.3121383247e-02, 3.3713201787e+05, 2.9420361824e+02, 9.3609829468e+02, 4.9806218220e+00,
                       1.4296138729e-02, 1.4480144748e-01, 3.1416000000e+00, 5.5833000000e+01, 3.7519530260e-08)),
  m1273 = list(label = "mRNA-1273 (Moderna)", dose = 100, interval = 28,
               fit = c(9.6660835849e+01, 1.4252521283e+01, 9.7872570305e+00, 2.9111966684e+00, 2.3888070078e+01,
                       1.4244131319e-01, 3.3713258652e+05, 3.0559394301e+02, 6.3759266444e+02, 6.2399530409e+00,
                       6.9960468437e-03, 6.5068664864e-02, 3.1416002938e+00, 5.5832999644e+01, 8.1199799530e-04))
)

MRNA_MW <- 1377479.8   # g/mol, BNT162b2 mRNA sequence (used for both vaccines in the source)
ug_to_pmol <- function(ug) 1e12 * ug * 1e-6 / MRNA_MW
OUT_TIMES <- c(seq(0.25, 2, by = 0.25), 3:365)

#' Parameters and doses for a schedule; booster = NA for none
vaccine_setup <- function(key, dose, interval, booster = NA) {
  v <- VACCINES[[key]]
  P <- as.data.frame(as.list(stats::setNames(v$fit, AB_NAMES)))
  P$TD2 <- interval
  P$TD3 <- if (is.na(booster)) 1e9 else booster
  days <- c(0, interval, if (!is.na(booster)) booster)
  list(P = P, ev = data.frame(time = days, cmt = "mRNA", amt = ug_to_pmol(dose)))
}

#' The outputs kept for plotting
vaccine_outputs <- function(r) {
  ab <- r$states[, paste0("Ab_BL", 1:17), 1]
  list(time = r$time, IGG = r$IGG[, 1], AFFINITY = r$AFFINITY[, 1], GCBTOT = r$GCBTOT[, 1], MBTOT = r$MBTOT[, 1],
       SPTOT = r$SPTOT[, 1], LPTOT = r$LPTOT[, 1], APC_IS = r$APC_IS[, 1], APC_LN = r$APC_LN[, 1],
       THELP = r$THELP[, 1], TACT = r$TACT[, 1], TMEM = r$TMEM[, 1], AB = ab)
}

#' Standard schedules shipped precomputed (tools/precompute_mrna_presets.R)
PRESET_RUNS <- list(
  list(key = "bnt", dose = 30, interval = 21), list(key = "bnt", dose = 20, interval = 21),
  list(key = "bnt", dose = 10, interval = 21), list(key = "bnt", dose = 1, interval = 21),
  list(key = "bnt60", dose = 30, interval = 21), list(key = "m1273", dose = 100, interval = 28)
)
preset_id <- function(key, dose, interval, booster = NA) {
  sprintf("%s_%g_%g_%s", key, dose, interval, if (is.na(booster)) "none" else format(booster))
}
