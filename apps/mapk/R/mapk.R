# ============================================================================
# Kirouac 2017 MAPK model: virtual patients, dosing and RECIST responses.
#
# Follows the authors' VPop_Simulator script (Supplementary Material):
#   - virtual patients are drawn from the 1000 parameter sets of Table S10,
#     weighted by their prevalence weights (PW)
#   - drug PK parameters are then re-drawn from the population PK models
#     of Table S9 (log-normal, full covariance), overriding Table S10's PK
#   - responses are read from tumour size (CELLS, relative to baseline) at
#     day 56: CR < 0.1, PR < 0.7, SD <= 1.2, PD > 1.2
# ============================================================================

MAPK_DRUGS <- list(
  CETUX = list(label = "Cetuximab (EGFR)", cmt = "RTK1i_blood", dose = 450, schedule = "weekly"),
  VEMU  = list(label = "Vemurafenib (BRAF)", cmt = "RAFi_gut", dose = 960, schedule = "bid"),
  COBI  = list(label = "Cobimetinib (MEK)", cmt = "MEKi_gut", dose = 60, schedule = "qd21of28"),
  ERKI  = list(label = "GDC-0994 (ERK)", cmt = "ERKi_gut", dose = 400, schedule = "qd21of28")
)

mapk_data <- function(file) utils::read.csv(repo_file_app("data", file), check.names = FALSE)

# data/ sits next to app.R in both layouts
repo_file_app <- function(...) file.path(...)

#' Draw virtual patients
#'
#' @return data frame of model parameters, one row per patient
mapk_patients <- function(n, model, seed = 1) {
  vp <- mapk_data("vpop_kirouac2017.csv")
  set.seed(seed)
  idx <- sample(nrow(vp), n, replace = TRUE, prob = abs(vp$PW))
  P <- vp[idx, intersect(names(vp), names(model$param)), drop = FALSE]
  rownames(P) <- NULL

  # Population PK (Table S9): P = THETA * exp(ETA), ETA ~ MVN(0, OMEGA)
  mvn <- function(n, omega) matrix(stats::rnorm(n * nrow(omega)), n) %*% chol(omega)
  vemu <- exp(mvn(n, matrix(c(0.109, 0, 0.106, 0, 1.02, 0, 0.106, 0, 0.457), 3)))
  P$ke2 <- 32.2 * vemu[, 1] / (114 * vemu[, 3])
  P$ka2 <- 4.48 * vemu[, 2]
  P$V2 <- 114 * vemu[, 3]
  cobi_om <- matrix(c(0.35, 0, 0, 0.264, 0.18,
                      0, 2.81, 0, 0, 0,
                      0, 0, 0.517, 0, 0,
                      0.264, 0, 0, 0.26, 0.0625,
                      0.18, 0, 0, 0.0625, 0.557), 5, byrow = TRUE)
  cobi <- exp(mvn(n, cobi_om))
  P$ke3 <- 327 * cobi[, 1] / (487 * cobi[, 4])
  P$ka3 <- 33.7 * cobi[, 2]
  P$q2 <- 252 * cobi[, 3]
  P$V3 <- 487 * cobi[, 4]
  P$V3b <- 335 * cobi[, 5]
  erki <- exp(mvn(n, diag(c(0.182, 3, 0.149))))
  P$ke4 <- 161 * erki[, 1] / (171 * erki[, 3])
  P$ka4 <- 35 * erki[, 2]
  P$V4 <- 171 * erki[, 3]
  P
}

#' Dose records for a set of drugs over `days`
mapk_events <- function(drugs, doses, days = 56) {
  cyc <- unlist(lapply(seq(0, days - 1, by = 28), function(s) s + 0:20))
  cyc <- cyc[cyc < days]
  ev <- lapply(drugs, function(d) {
    x <- MAPK_DRUGS[[d]]
    t <- switch(x$schedule,
                weekly = seq(0, days - 1, by = 7),
                bid = seq(0, days - 0.5, by = 0.5),
                qd21of28 = cyc)
    data.frame(time = t + 0.001, cmt = x$cmt, amt = doses[[d]])
  })
  if (!length(ev)) return(NULL)
  do.call(rbind, ev)
}

recist <- function(cells) {
  cut(cells, c(-Inf, 0.1, 0.7, 1.2, Inf), right = FALSE,
      labels = c("Complete response", "Partial response", "Stable disease", "Progressive disease"))
}

RECIST_COLS <- c("Complete response" = "#0d7a54", "Partial response" = "#1baf7a",
                 "Stable disease" = "#eda100", "Progressive disease" = "#d1453b")
