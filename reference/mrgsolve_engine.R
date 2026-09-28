# ============================================================================
# mrgsolve engine: runs models/tmdd.cpp through the same interface as
# tmdd_simulate() in R/tmdd_engine.R.
#
# The app sources this file when it runs locally. If mrgsolve is installed
# and the model compiles, simulations go through mrgsolve; otherwise the
# app keeps the base-R solver. The browser build never includes this
# directory, since a browser can't compile C++.
#
# tests/test_engine_vs_mrgsolve.R uses it as the reference.
# ============================================================================

mrgsolve_ready <- function() {
  requireNamespace("mrgsolve", quietly = TRUE)
}

.tmdd_mod <- NULL

tmdd_mrgsolve_model <- function(model_dir = "models") {
  if (is.null(.tmdd_mod)) {
    .tmdd_mod <<- mrgsolve::mread("tmdd", project = model_dir, quiet = TRUE)
  }
  .tmdd_mod
}

#' Same contract as tmdd_simulate()
tmdd_simulate_mrgsolve <- function(P, regimen, times, model_dir = "models",
                                   rtol = 1e-10, atol = 1e-12) {
  mod <- tmdd_mrgsolve_model(model_dir)
  n <- nrow(P)
  times <- sort(unique(times))
  reg <- expand_tmdd_regimen(regimen, P, max(times))
  cmt <- TMDD_STATES[reg$cmt]

  ev_df <- do.call(rbind, lapply(seq_along(reg$times), function(i) {
    data.frame(ID = P$ID, time = reg$times[i], amt = reg$amounts[[i]], cmt = cmt, evid = 1)
  }))
  ev_df <- ev_df[order(ev_df$ID, ev_df$time), ]
  pars <- c("CL", "V1", "Q", "V2", "KA", "F", "MW", "KSS", "KINT", "KDEG", "R0", "B0", "KOUT", "GAMMA", "IMAX")

  out <- mrgsolve::mrgsim_df(
    mrgsolve::update(mod, rtol = rtol, atol = atol, maxsteps = 1e6),
    data = ev_df, idata = P[, c("ID", pars)], add = times, end = -1,
    obsonly = TRUE, carry_out = character(0), recover = character(0)
  )

  nt <- length(times)
  store <- array(NA_real_, c(nt, 5, n))
  for (k in seq_along(TMDD_STATES)) store[, k, ] <- matrix(out[[TMDD_STATES[k]]], nrow = nt)
  tmdd_outputs(times, store, P)
}
