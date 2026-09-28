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

cat("PASS\n")
