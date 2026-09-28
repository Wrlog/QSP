# Every model, solved by the browser engine (shared/ode_engine.R, at the
# tolerances the apps use) and by mrgsolve compiling the same model file,
# must agree. If they do, the browser apps compute what mrgsolve computes.
#
# Run from the repository root: Rscript tests/test_engine_vs_mrgsolve.R

suppressMessages(library(shiny))
for (f in c("ode_engine.R", "app_helpers.R")) source(file.path("shared", f))
source(file.path("reference", "mrgsolve_engine.R"))
if (!mrgsolve_ready()) stop("mrgsolve is not installed")

rel <- function(a, b) max(abs(a - b)) / max(abs(b))
worst <- c()
report <- function(name, d) {
  cat(sprintf("  %-40s %.2e\n", name, d))
  worst[name] <<- d
}

# Coagulation: 14 days of warfarin 5 mg, and the in-vitro PT test
m <- mrg_read("models/coagulation_wajima2009.cpp")
source(file.path("apps", "coagulation", "R", "coagulation.R"))
tt <- seq(0, 24 * 21, by = 24)
ev <- coag_events(5, 14)
a <- mrg_solve(m, events = ev, times = tt, rtol = 1e-5, atol = 1e-8, nonneg = TRUE)
b <- mrg_solve_mrgsolve(m, events = ev, times = tt)
for (v in c("FII", "VII", "IX", "X", "PC", "C_warf")) report(paste("coagulation", v), rel(a$states[, v, 1], b$states[, v, 1]))
prep <- mrg_prepare(m, data.frame(INVITRO = 1), nonneg = TRUE)
pt_r <- sapply(c(1, 15), function(i) coag_pt(prep, a$states[i, , 1]))
pt_m <- sapply(c(1, 15), function(i) coag_pt_mrgsolve(m, b$states[i, , 1]))
report("coagulation prothrombin time", rel(pt_r, pt_m))

# Bone: denosumab 60 mg every 6 months, 2 years
m <- mrg_read("models/bone_peterson_riggs.cpp")
tt <- seq(0, 24 * 730, by = 24 * 14)
ev <- data.frame(time = 0, cmt = "DENSC", amt = 60, ii = 24 * 182, addl = 3)
a <- mrg_solve(m, events = ev, times = tt, rtol = 1e-4, atol = 1e-8)
b <- mrg_solve_mrgsolve(m, events = ev, times = tt)
report("bone BMD (percentage points / 100)", max(abs(a$BMDlsDENchange - b$BMDlsDENchange)) / 100)
report("bone CTx", rel(a$OCchange, b$OCchange))
report("bone serum calcium", rel(a$CaC, b$CaC))

# MAPK: 8 virtual patients on cetuximab + vemurafenib + cobimetinib
m <- mrg_read("models/mapk_kirouac2017.cpp")
owd <- setwd(file.path("apps", "mapk")); source(file.path("R", "mapk.R"))
P <- mapk_patients(8, m, 3)
setwd(owd)
ev <- mapk_events(c("CETUX", "VEMU", "COBI"), list(CETUX = 450, VEMU = 960, COBI = 60))
tt <- seq(0, 56, by = 2)
a <- mrg_solve(m, P, ev, tt, rtol = 1e-4, atol = 1e-7)
b <- mrg_solve_mrgsolve(m, P, ev, tt)
report("mapk tumour size", rel(a$states[, "CELLS", ], b$states[, "CELLS", ]))
report("mapk ERK", rel(a$ERK, b$ERK))

# Glucose: type 2 diabetes, liraglutide 1.8 mg daily, three meals on day 14
m <- mrg_read("models/glucose_4gi_bosch2022.cpp")
ev <- rbind(data.frame(time = 0, cmt = "Ddrug", amt = 1.8e9 / 3751.2, ii = 24, addl = 13),
            data.frame(time = 13 * 24 + c(8, 13, 19), cmt = "Dglc", amt = c(50, 70, 80) / 180.16 * 1000 * c(0.446, 0.294, 0.287), ii = 0, addl = 0))
tt <- seq(0, 14 * 24, by = 0.5)
for (pat in 0:1) {
  a <- mrg_solve(m, data.frame(PAT = pat), ev, tt, rtol = 1e-5, atol = 1e-8)
  b <- mrg_solve_mrgsolve(m, data.frame(PAT = pat), ev, tt)
  for (v in c("GLC", "INS", "GLP1", "GLG", "GIP")) report(sprintf("glucose %s (PAT=%d)", v, pat), rel(a[[v]], b[[v]]))
}

cat(sprintf("comparisons: %d, worst: %.2e\n", length(worst), max(worst)))
if (!all(is.finite(worst)) || max(worst) > 2e-3) stop("browser engine disagrees with mrgsolve")
cat("PASS\n")
