# ============================================================================
# Wajima 2009 coagulation network: dosing, and the prothrombin time test.
#
# The in-vivo model (models/coagulation_wajima2009.cpp) predicts clotting
# factor concentrations over days. The clinical read-out, PT/INR, is an
# in-vitro test on a plasma sample, and Wajima et al. simulated it the same
# way: take the plasma composition at the sampling time, dilute it 3-fold,
# add tissue factor (300 nM before dilution) and find the time at which the
# integral of fibrin reaches 1500 nM*s. That is what coag_pt() does.
# ============================================================================

COAG_MW_VITK <- 450.7        # phytonadione (g/mol)
COAG_VC_VITK <- 24           # vitamin K central volume in the model (L)

#' Dose records: daily oral warfarin and optional IV vitamin K
#'
#' @param warfarin_mg daily dose (mg); warfarin goes into the gut amount A_warf
#' @param days number of daily doses
#' @param vitk_mg,vitk_day an optional single IV vitamin K dose and its day
coag_events <- function(warfarin_mg, days, vitk_mg = 0, vitk_day = NA) {
  ev <- data.frame(time = numeric(0), cmt = character(0), amt = numeric(0),
                   ii = numeric(0), addl = numeric(0))
  if (warfarin_mg > 0 && days > 0) {
    ev <- rbind(ev, data.frame(time = 0, cmt = "A_warf", amt = warfarin_mg,
                               ii = 24, addl = days - 1))
  }
  if (vitk_mg > 0 && is.finite(vitk_day)) {
    # The model's vitamin K state is a concentration (nM) in a 24 L volume.
    ev <- rbind(ev, data.frame(time = vitk_day * 24, cmt = "VK",
                               amt = vitk_mg / COAG_MW_VITK * 1e6 / COAG_VC_VITK,
                               ii = 0, addl = 0))
  }
  ev
}

#' Starting state for the in-vitro test from a plasma state
coag_pt_start <- function(state) {
  st <- pmax(state, 0) / 3             # 1:3 dilution of the sample
  st["TF"] <- 300 / 3                  # thromboplastin (tissue factor)
  st["AUC_Fibrin"] <- 0
  st
}

#' Prothrombin time (s) with the base-R engine
#'
#' Steps one second at a time and stops as soon as the fibrin integral
#' crosses 1500 nM*s, so a normal sample costs ~11 s of simulated time
#' rather than a fixed window.
#'
#' @param prep from mrg_prepare(model, data.frame(INVITRO = 1), nonneg = TRUE)
coag_pt <- function(prep, state, max_s = 180, rtol = 1e-5, atol = 1e-8) {
  Y <- prep$init
  st <- coag_pt_start(state)
  Y[names(st), 1] <- st
  auc_row <- which(rownames(Y) == "AUC_Fibrin")
  h <- 1e-7
  prev <- 0
  for (s in seq_len(max_s)) {
    seg <- mrg_advance(prep, Y, (s - 1) / 3600, s / 3600, h, rtol, atol)
    Y <- seg$Y
    h <- seg$h
    auc <- Y[auc_row, 1] * 3600
    if (auc >= 1500) return(s - 1 + (1500 - prev) / (auc - prev))
    prev <- auc
  }
  NA_real_
}

#' Prothrombin time (s) with mrgsolve, for the tests
coag_pt_mrgsolve <- function(model, state, max_s = 180, step_s = 0.25) {
  tt <- seq(0, max_s, by = step_s) / 3600
  r <- mrg_solve_mrgsolve(model, data.frame(INVITRO = 1), NULL, tt,
                          init = coag_pt_start(state), rtol_ref = 1e-8, atol_ref = 1e-10)
  auc <- r$states[, "AUC_Fibrin", 1] * 3600
  k <- which(auc >= 1500)[1]
  if (is.na(k)) return(NA_real_)
  (tt[k - 1] + (1500 - auc[k - 1]) * (tt[k] - tt[k - 1]) / (auc[k] - auc[k - 1])) * 3600
}
