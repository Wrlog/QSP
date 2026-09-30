# ============================================================================
# T-DM1 (Singh & Shah, AAPS J 2017): parameter sets and dosing helpers
# ============================================================================

# Table I: tumour growth (exponential doubling time) and killing constant
# (mean and inter-individual variability, CV) estimated for each mouse model;
# HER2 status from Table II, converted to antigen per the paper
# (3+ ~1 million, 2+ ~0.5 million, 1+ ~0.1 million receptors per cell).
MOUSE_MODELS <- data.frame(
  model = c("BT474EEI", "Fo5", "Calu-3", "KPL-4", "N87", "BT474", "MAXF 1162",
            "HBCx-34 HP", "ST313 HP", "HBCx-10", "MAXF 449"),
  type = c("Xenograft", "Orthotopic (BOM)", "Xenograft", "Xenograft", "Xenograft", "Xenograft", "PDX",
           "PDX", "PDX", "PDX", "PDX (TNBC)"),
  her2 = c("3+", "3+", "3+", "3+", "3+", "3+", "3+", "2+", "2+", "1+", "1+"),
  dt = c(21.0, 4.5, 8.99, 8.21, 14.3, 11.2, 10.2, 11.4, 10.5, 6.84, 13.8),
  kkill = c(6.83e-5, 2.38e-4, 4.86e-4, 1.96e-3, 6.38e-4, 9e-4, 5.58e-4, 1.8e-4, 1.93e-4, 1.43e-4, 1.57e-4),
  iiv = c(0.28, 0.34, 0.41, 0.21, 0.32, 0.31, 0.23, 0.18, 0.11, 0.14, 0),
  studied = c("0.3-15 mg/kg Q3W x3; 3.3-18 mg/kg Q1W x9; 6-18 mg/kg Q2W x5", "1-30 mg/kg Q3W x3; 3.3-10 mg/kg Q1W x9",
              "1-7 mg/kg single; 15 mg/kg Q6D x3", "0.3-3 and 15 mg/kg single", "1-10 and 5 mg/kg single",
              "0.2-5 mg/kg single", "1-10 mg/kg single", "3-30 mg/kg single", "3-30 mg/kg single",
              "3-30 mg/kg single", "3-30 mg/kg single"),
  stringsAsFactors = FALSE
)

HER2_AG <- c("3+" = 1660, "2+" = 830, "1+" = 166)

# Table I systemic PK: mouse (Singh et al. 2016) and human (allometrically
# scaled from the monkey fit: exponent 1 for the antibody, 0.75 for DM1
# clearances). Kdec is the same in all species.
PK_SETS <- list(
  mouse = list(BW = 0.025, CLADC = 0.0934, CLDADC = 0.118, V1ADC = 0.043, V2ADC = 0.0948,
               CLDRUG = 11.29, CLDDRUG = 155, V1DRUG = 3.30, V2DRUG = 2.01, KDEC = 0.241),
  human = list(BW = 70, CLADC = 0.0043, CLDADC = 0.014, V1ADC = 0.034, V2ADC = 0.04,
               CLDRUG = 2.23, CLDDRUG = 1.0, V1DRUG = 0.034, V2DRUG = 5.0, KDEC = 0.241)
)

# Clinical killing constants used for translation (Results): HER2 1+ and 3+,
# 60% inter-individual variability
KKILL_CLIN <- c("3+" = 1.8e-4, "2+" = NA, "1+" = 3.8e-5)

#' Doses: T-DM1 into total and conjugated antibody, and the average DAR set
#' back to its initial value at every dose. The paper gives DAR only an
#' initial condition; for repeated doses the circulating ADC is almost all
#' from the latest dose (T-DM1 falls ~50-fold between doses), so its DAR is
#' that of fresh drug.
tdm1_doses <- function(times, mgkg, bw, mw = 148500, dar0 = 3.5, ID = NULL) {
  nmol <- mgkg * bw * 1e6 / mw
  ev <- rbind(data.frame(time = times, cmt = "X1TT", amt = nmol, evid = 1),
              data.frame(time = times, cmt = "X1ADC", amt = nmol, evid = 1),
              data.frame(time = times, cmt = "DAR", amt = dar0, evid = 8))
  ev <- ev[ev$amt > 0 | ev$cmt == "DAR", ]
  ev <- ev[order(ev$time), ]
  if (!is.null(ID)) ev$ID <- ID
  ev
}

#' The regimens compared in the paper (Fig. 6), per 21-day cycle unless noted
REGIMENS <- list(
  "3.6 mg/kg every 3 weeks (approved)" = function(cycles) data.frame(time = 21 * (seq_len(cycles) - 1), mgkg = 3.6),
  "3.6 mg/kg every 4 weeks" = function(cycles) data.frame(time = 28 * (seq_len(cycles) - 1), mgkg = 3.6),
  "1.2 mg/kg weekly" = function(cycles) data.frame(time = 7 * (seq_len(3 * cycles) - 1), mgkg = 1.2),
  "3 + 0.3 + 0.3 mg/kg on days 0, 7, 14" = function(cycles) {
    s <- 21 * (seq_len(cycles) - 1)
    data.frame(time = c(rbind(s, s + 7, s + 14)), mgkg = rep(c(3, 0.3, 0.3), cycles))
  }
)

regimen_events <- function(reg, bw, ID = NULL) {
  do.call(rbind, lapply(seq_len(nrow(reg)), function(i) tdm1_doses(reg$time[i], reg$mgkg[i], bw, ID = ID)))
}

#' Log-normal draws with a given median and CV
draw_ln <- function(n, median, cv) if (cv > 0) median * exp(stats::rnorm(n, 0, sqrt(log(1 + cv^2)))) else rep(median, n)
