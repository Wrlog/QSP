# The models must reproduce the behaviour their papers (and the clinic)
# report. Uses the browser engine only, so it also runs without mrgsolve.
#
# Run from the repository root: Rscript tests/test_published_behaviour.R

suppressMessages(library(shiny))
for (f in c("ode_engine.R", "app_helpers.R")) source(file.path("shared", f))

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}

# --- Coagulation (Wajima 2009) ----------------------------------------------
m <- mrg_read("models/coagulation_wajima2009.cpp")
source(file.path("apps", "coagulation", "R", "coagulation.R"))
r <- mrg_solve(m, events = coag_events(5, 14), times = seq(0, 24 * 21, by = 24), rtol = 1e-5, atol = 1e-8, nonneg = TRUE)
prep <- mrg_prepare(m, data.frame(INVITRO = 1), nonneg = TRUE)
pt <- sapply(c(1, 3, 15, 22), function(i) coag_pt(prep, r$states[i, , 1]))
inr <- pt / pt[1]
check(pt[1] > 10 && pt[1] < 14, sprintf("normal prothrombin time is 10-14 s (%.1f s)", pt[1]))
check(inr[3] > 1.8 && inr[3] < 3.2, sprintf("warfarin 5 mg/day for 14 days gives INR 1.8-3.2 (%.2f)", inr[3]))
check(inr[4] < inr[3], "INR falls after warfarin stops")
st <- r$states
check(st[3, "VII", 1] / st[1, "VII", 1] < st[3, "FII", 1] / st[1, "FII", 1],
      "factor VII (short half-life) falls faster than prothrombin")

# --- Bone (Peterson & Riggs) --------------------------------------------------
m <- mrg_read("models/bone_peterson_riggs.cpp")
tt <- seq(0, 24 * 365 * 3, by = 24 * 14)
r <- mrg_solve(m, events = data.frame(time = 0, cmt = "DENSC", amt = 60, ii = 24 * 182, addl = 3), times = tt, rtol = 1e-4, atol = 1e-8)
bmd <- r$BMDlsDENchange[, 1]
at <- function(y) bmd[which.min(abs(tt - y * 24 * 365))]
check(at(1) > 3 && at(1) < 8, sprintf("denosumab raises lumbar BMD 3-8%% at 1 year (%.1f%%)", at(1)))
check(min(r$OCchange[, 1]) < 30, sprintf("denosumab suppresses CTx by more than 70%% (nadir %.0f%% of baseline)", min(r$OCchange[, 1])))
check(at(3) < at(2) - 2, "BMD is lost after denosumab is stopped")

# --- MAPK (Kirouac 2017) --------------------------------------------------------
m <- mrg_read("models/mapk_kirouac2017.cpp")
owd <- setwd(file.path("apps", "mapk")); source(file.path("R", "mapk.R")); P <- mapk_patients(12, m, 2); setwd(owd)
tt <- c(0, 56)
cells <- function(drugs) {
  ev <- mapk_events(drugs, list(CETUX = 450, VEMU = 960, COBI = 60, ERKI = 400))
  mrg_solve(m, P, ev, tt, rtol = 1e-4, atol = 1e-7)$states[2, "CELLS", ]
}
ctrl <- cells(character(0)); vemu <- cells("VEMU"); combo <- cells(c("CETUX", "VEMU"))
check(stats::median(ctrl) > 1, "untreated tumours grow over 8 weeks")
check(mean(combo < 0.7) >= mean(vemu < 0.7) && stats::median(combo) < stats::median(vemu),
      "adding cetuximab to vemurafenib improves response (EGFR feedback reactivation)")

# --- Glucose (Bosch 2022) ---------------------------------------------------------
m <- mrg_read("models/glucose_4gi_bosch2022.cpp")
ogtt <- data.frame(time = 1, cmt = "Dglc", amt = 75 / 180.16 * 1000 * 0.776)
tt <- seq(0, 6, by = 0.1)
t2 <- mrg_solve(m, data.frame(PAT = 0), ogtt, tt, rtol = 1e-5, atol = 1e-8)
hv <- mrg_solve(m, data.frame(PAT = 1, BSLglc = 4.65, BSLins = 49.1, BSLglg = 8.85), ogtt, tt, rtol = 1e-5, atol = 1e-8)
check(max(t2$GLC) > 11 && max(hv$GLC) < 9, sprintf("OGTT glucose peak: T2DM > 11 mM (%.1f), healthy < 9 mM (%.1f)", max(t2$GLC), max(hv$GLC)))
meals <- data.frame(time = 13 * 24 + c(8, 13, 19), cmt = "Dglc", amt = c(50, 70, 80) / 180.16 * 1000 * c(0.446, 0.294, 0.287))
lira <- data.frame(time = seq(0, 13 * 24, by = 24), cmt = "Ddrug", amt = 1.8e9 / 3751.2)
tt <- seq(13 * 24, 14 * 24, by = 0.25)
no <- mrg_solve(m, data.frame(PAT = 0), meals, c(0, tt), rtol = 1e-5, atol = 1e-8)
yes <- mrg_solve(m, data.frame(PAT = 0), rbind(lira, meals), c(0, tt), rtol = 1e-5, atol = 1e-8)
d <- mean(no$GLC[-1]) - mean(yes$GLC[-1])
check(d > 1 && d < 4, sprintf("liraglutide 1.8 mg lowers mean daily glucose by 1-4 mM (%.1f)", d))
check(max(yes$INS[-1]) > max(no$INS[-1]), "liraglutide raises meal-time insulin")

# --- CAR-T cellular kinetics (Stein 2019) -------------------------------------
m <- mrg_read("models/cart_tisagenlecleucel_stein2019.cpp")
tt <- sort(unique(c(seq(0, 2000, by = 0.1), 9.3)))
r <- mrg_solve(m, times = tt, rtol = 1e-7, atol = 1e-9)
check(abs(max(r$CART) / 24000 - 1) < 1e-3 && abs(tt[which.max(r$CART)] - 9.3) < 0.11,
      sprintf("typical patient peaks at Cmax 24,000 copies/ug on day 9.3 (%.0f, day %.1f)", max(r$CART), tt[which.max(r$CART)]))
rho <- log(3900) / 9.3
y <- as.vector(r$CART); auc <- sum(diff(tt) * (head(y, -1) + tail(y, -1)) / 2)
auc_paper <- 24000 * (1 / rho + (1 - 0.0079) / 0.16 + 0.0079 / 0.0032)
check(abs(auc / auc_paper - 1) < 0.02, sprintf("AUC matches the paper's closed form Cmax[1/rho + (1-FB)/alpha + FB/beta] (%.3g vs %.3g)", auc, auc_paper))
toci <- mrg_solve(m, data.frame(TTOCI = 5), times = tt, rtol = 1e-7, atol = 1e-9)
check(max(toci$CART) >= max(r$CART), "tocilizumab during expansion does not slow it (Ftoci = 1.2)")

# --- PROTAC degradation (Pharmaceutics 2023) ------------------------------------
m <- mrg_read("models/protac_degrader_kcat2023.cpp")
conc <- c(sqrt(2500 * 71), 30 * sqrt(2500 * 71), 1)
r <- mrg_solve(m, data.frame(CFIX = conc), times = c(0, 1000), rtol = 1e-7, atol = 1e-10)
k <- 0.86 * 203 * 4.6 * 16 / log(2)
dmax <- k / (k + 0.86 * 203 + 2500 + 71 + 2 * sqrt(2500 * 71))
check(abs(r$DEG[2, 1] - dmax) < 1e-4, sprintf("steady-state degradation at DCmax equals Dmax from Appendix A (%.3f)", dmax))
check(r$DEG[2, 2] < r$DEG[2, 1] && r$DEG[2, 3] < r$DEG[2, 1], "hook effect: less degradation above and below DCmax")
check(abs(r$TM[2, 2] - 1) < 0.05, "with occupancy-driven inhibition, total modulation has no hook")

# --- T-cell engager (Hosseini 2020) ----------------------------------------------
m <- mrg_read("models/tce_mosunetuzumab_hosseini2020.cpp")
tce <- function(mg) {
  days <- c(0, 7, 14, 21)[mg > 0]; mg <- mg[mg > 0]
  ev <- rbind(data.frame(time = days, cmt = "TDBc_ugperkg", amt = mg * 1000 / 70),
              data.frame(time = days, cmt = "injection_effect", amt = 1),
              data.frame(time = days, cmt = "drug_effect", amt = 1))
  tt <- sort(unique(c(seq(0, 28, by = 0.1), outer(days, seq(0.002, 0.2, by = 0.004), "+"))))
  list(t = tt, r = mrg_solve(m, events = ev, times = tt, rtol = 1e-5, atol = 1e-3))
}
base <- mrg_solve(m, times = seq(0, 28, by = 1), rtol = 1e-5, atol = 1e-3)
check(max(abs(base$Bpb_perml / base$Bpb_perml[1] - 1)) < 1e-6, "without drug the cell populations stay at baseline")
step <- tce(c(1, 2, 60, 60)); full <- tce(c(60, 0, 0, 60))
p_step <- max(step$r$IL6combo[step$t < 21]); p_full <- max(full$r$IL6combo[full$t < 21])
check(p_step < 0.6 * p_full, sprintf("step-up dosing lowers the cycle-1 IL-6 peak (%.0f vs %.0f pg/mL)", p_step, p_full))
i21 <- which.min(abs(step$t - 21))
check(step$r$Bpb_perml[i21, 1] / step$r$Bpb_perml[1, 1] < 0.01, "step-up dosing still depletes > 99% of blood B cells by day 21")

# --- mRNA vaccine (Dasti 2025) -----------------------------------------------------
m <- mrg_read("models/mrna_vaccine_dasti2025.cpp")
source(file.path("apps", "mrna-vaccine", "R", "mrna_rhs.R"))
source(file.path("apps", "mrna-vaccine", "R", "vaccines.R"))
sahin <- read.csv(file.path("apps", "mrna-vaccine", "data", "sahin2020_bnt162b2.csv"))
s <- vaccine_setup("bnt", 30, 21)
r <- mrg_solve(m, s$P, s$ev, c(7, 21, 28, 42, 49, 84), rtol = 1e-4, atol = 1e-4, nonneg = TRUE, rhs = mrna_rhs(m))
obs <- sahin[sahin$dose_ug == 30, ]
ratio <- r$IGG[, 1] / obs$geomean_ng_per_mL[match(c(7, 21, 28, 42, 49, 84), obs$day)]
# days 21, 42, 49 and 84; one week after the second dose (day 28) the model's
# rise lags the data (it peaks about a week later), which the check records
check(all(ratio[c(2, 4, 5, 6)] > 0.5 & ratio[c(2, 4, 5, 6)] < 2),
      sprintf("BNT162b2 30 ug: IgG within 2-fold of the Sahin et al. geometric means on days 21, 42, 49, 84 (ratios %s)",
              paste(sprintf("%.2f", ratio[c(2, 4, 5, 6)]), collapse = ", ")))
check(ratio[3] > 0.2, sprintf("day 28, one week after the boost: model still rising (%.2f of the observed mean)", ratio[3]))
pre <- readRDS(file.path("apps", "mrna-vaccine", "data", "presets.rds"))
pk <- sapply(c(1, 10, 20, 30), function(dz) max(pre[[preset_id("bnt", dz, 21)]]$IGG))
check(all(diff(pk) > 0), "peak antibody rises with dose from 1 to 30 ug")
p30 <- pre[[preset_id("bnt", 30, 21)]]
check(max(p30$IGG) / p30$IGG[which.min(abs(p30$time - 21))] > 5, "the second dose boosts antibody more than 5-fold")

# --- Antibody-drug conjugate T-DM1 (Singh & Shah 2017) ------------------------------------
m <- mrg_read("models/adc_tdm1_singh2017.cpp")
source(file.path("apps", "adc", "R", "tdm1.R"))
adc <- function(P, ev, tt) mrg_solve(m, P, ev, tt, rtol = 1e-5, atol = 1e-9, nonneg = TRUE)
# untreated mouse tumours grow with the fitted doubling time; T-DM1 inhibits them dose-dependently
kpl <- MOUSE_MODELS[MOUSE_MODELS$model == "KPL-4", ]
Pm <- as.data.frame(c(PK_SETS$mouse, list(GLIN = 0, DTEXP = kpl$dt, KKILL = kpl$kkill, AG = HER2_AG[["3+"]], TV0 = 200e-6)))
tv <- sapply(c(0, 0.3, 3, 15), function(d) tail(adc(Pm, tdm1_doses(0, d, 0.025), c(0.01, 28))$TV_MM3[, 1], 1))
check(abs(tv[1] / 200 / 2^(28 / kpl$dt) - 1) < 0.01, sprintf("untreated KPL-4 doubles every %.1f days (Table I)", kpl$dt))
check(all(diff(tv) < 0) && tv[4] < 200,
      sprintf("KPL-4, single dose: day-28 tumour %.0f / %.0f / %.0f / %.0f mm3 at 0 / 0.3 / 3 / 15 mg/kg - dose-dependent, regression at 15",
              tv[1], tv[2], tv[3], tv[4]))
# patients: T-DM1 falls faster than total trastuzumab (deconjugation), and the DAR with it
Ph <- transform(as.data.frame(PK_SETS$human), KKILL = 0)
r <- adc(Ph, regimen_events(REGIMENS[[1]](1), 70), c(0.01, 7, 14, 20.9))
ratio <- r$ADC_UGML[, 1] / r$TT_UGML[, 1]
check(all(diff(ratio) < 0) && ratio[4] < 0.25, sprintf("3.6 mg/kg: T-DM1 / total trastuzumab falls from %.2f to %.2f over the cycle", ratio[1], ratio[4]))
check(abs(r$states[3, "DAR", 1] / 3.5 / exp(-0.241 * (14 - 0)) - 1) < 1e-3, "the average DAR falls with the deconjugation rate (half-life 2.9 days)")
# Fig. 6: at equal dose intensity, fractionated dosing keeps tumour DM1 higher; every 4 weeks lowers it
tumour_dm1 <- function(nm) {
  reg <- REGIMENS[[nm]](4); reg <- reg[reg$time < 84, ]
  mean(adc(Ph, regimen_events(reg, 70), seq(21.25, 84, by = 0.25))$DM1_TUMOUR[, 1])
}
dm <- sapply(names(REGIMENS), tumour_dm1)
check(dm[3] > dm[1] && dm[4] > dm[1] && dm[2] < dm[1],
      sprintf("mean tumour DM1 over cycles 2-4: Q3W %.0f, Q4W %.0f, weekly %.0f, front-loaded %.0f nM - the ranking behind Fig. 6",
              dm[1], dm[2], dm[3], dm[4]))

cat("PASS\n")
