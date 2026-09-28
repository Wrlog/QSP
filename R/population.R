# ============================================================================
# Individual parameters, virtual populations and response summaries.
#
# The defaults describe a generic IgG1 against a membrane target. They are
# round numbers in the range typical of published mAb models (linear
# clearance about 0.2 L/day, central volume about 3 L, subcutaneous
# bioavailability about 70%), not the estimates for any real antibody.
# ============================================================================

TYPICAL <- list(
  # Antibody PK, at the reference weight
  CL = 0.20,        # linear clearance (L/day)
  V1 = 3.0,         # central volume (L)
  Q = 0.50,         # intercompartmental clearance (L/day)
  V2 = 2.5,         # peripheral volume (L)
  KA = 0.25,        # SC absorption rate (1/day)
  F = 0.70,         # SC bioavailability
  MW = 150000,      # molecular weight (g/mol)
  # Target
  KSS = 1.0,        # quasi-steady-state constant (nM)
  KINT = 1.0,       # complex internalisation rate (1/day)
  KDEG = 0.2,       # free target degradation rate (1/day)
  R0 = 2.0,         # baseline target concentration (nM)
  # Biomarker driven by free target
  B0 = 100,         # baseline biomarker (arbitrary units)
  KOUT = 0.3,       # biomarker turnover rate (1/day)
  GAMMA = 1.0,      # sensitivity of biomarker production to free target
  IMAX = 0.8,       # maximum fractional suppression of biomarker production
  # Covariates
  WT_REF = 70,
  WT_EXP_CL = 0.8,  # on CL and Q
  WT_EXP_V = 0.6    # on V1 and V2
)

cv_to_sd <- function(cv_percent) sqrt(log(1 + (cv_percent / 100)^2))

#' Individual parameter table
#'
#' Weight is drawn uniformly over wt_range. Between-subject variability is
#' log-normal on linear clearance, central volume, baseline target and the
#' absorption rate.
#'
#' @param cv named list of CVs (%) for CL, V1, R0, KA
tmdd_individual <- function(n, wt_range = c(70, 70), typ = TYPICAL,
                            cv = list(CL = 0, V1 = 0, R0 = 0, KA = 0), seed = 1) {
  set.seed(seed)
  wt <- if (diff(range(wt_range)) > 0) stats::runif(n, min(wt_range), max(wt_range)) else rep(wt_range[1], n)
  eta <- function(k) exp(stats::rnorm(n, 0, cv_to_sd(if (is.null(cv[[k]])) 0 else cv[[k]])))
  e_cl <- eta("CL"); e_v1 <- eta("V1"); e_r0 <- eta("R0"); e_ka <- eta("KA")

  fw_cl <- (wt / typ$WT_REF)^typ$WT_EXP_CL
  fw_v <- (wt / typ$WT_REF)^typ$WT_EXP_V
  data.frame(
    ID = seq_len(n), WT = wt,
    CL = typ$CL * fw_cl * e_cl,
    V1 = typ$V1 * fw_v * e_v1,
    Q = typ$Q * fw_cl,
    V2 = typ$V2 * fw_v,
    KA = typ$KA * e_ka,
    F = typ$F,
    MW = typ$MW,
    KSS = typ$KSS,
    KINT = typ$KINT,
    KDEG = typ$KDEG,
    R0 = typ$R0 * e_r0,
    B0 = typ$B0,
    KOUT = typ$KOUT,
    GAMMA = typ$GAMMA,
    IMAX = typ$IMAX
  )
}

#' Median and prediction-interval bands across subjects at each time
summarise_profiles <- function(time, x) {
  q <- apply(x, 1, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
             na.rm = TRUE, names = FALSE)
  if (is.null(dim(q))) q <- matrix(q, nrow = 5)
  data.frame(time = time, q05 = q[1, ], q25 = q[2, ], med = q[3, ],
             q75 = q[4, ], q95 = q[5, ])
}

#' Response over the final dosing interval, per subject
#'
#' Trough values are taken just before the next dose would be due (or at
#' the end of the simulation for the last interval).
final_interval <- function(res, regimen) {
  t_end <- max(res$time)
  doses <- (seq_len(regimen$n_doses) - 1) * regimen$interval
  doses <- doses[doses <= t_end + 1e-9]
  start <- max(doses)
  end <- if (regimen$n_doses == 1) t_end else min(start + regimen$interval, t_end)
  w <- which(res$time >= start - 1e-9 & res$time <= end + 1e-9)
  last <- w[length(w)]
  list(
    window = c(start, end),
    trough_ro = res$ro[last, ],
    trough_conc = res$conc[last, ],
    min_ro = apply(res$ro[w, , drop = FALSE], 2, min),
    bio_end = res$bio_rel[nrow(res$bio_rel), ],
    rfree_end = res$rfree_rel[nrow(res$rfree_rel), ]
  )
}

#' Dose-ranging: trough occupancy across a grid of doses
#'
#' All doses and subjects go through the solver as one population, so a
#' dose sweep costs one simulation, not one per dose.
dose_ranging <- function(doses_mgkg, P, regimen, times, target_ro,
                         simulate = tmdd_simulate) {
  nd <- length(doses_mgkg)
  n <- nrow(P)
  Pall <- P[rep(seq_len(n), nd), ]
  Pall$ID <- seq_len(nrow(Pall))
  reg <- regimen
  reg$amt <- rep(doses_mgkg, each = n) * Pall$WT
  reg$load <- NA
  res <- simulate(Pall, reg, times)
  fi <- final_interval(res, reg)
  dose_id <- rep(seq_len(nd), each = n)
  do.call(rbind, lapply(seq_len(nd), function(k) {
    ro <- fi$min_ro[dose_id == k]
    bio <- fi$bio_end[dose_id == k]
    data.frame(dose = doses_mgkg[k],
               ro_med = stats::median(ro), ro_lo = stats::quantile(ro, 0.05),
               ro_hi = stats::quantile(ro, 0.95),
               pta = mean(ro >= target_ro),
               bio_med = stats::median(bio),
               ctrough = stats::median(fi$trough_conc[dose_id == k]))
  }))
}
