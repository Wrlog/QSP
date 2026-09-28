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

# CAR-T: 40 virtual patients with between-patient variability, one year
m <- mrg_read("models/cart_tisagenlecleucel_stein2019.cpp")
set.seed(7)
om <- c(FOLDX = 2.4, TMAX = 0.38, CMAX = 0.65, ALPHA = 0.91, FB = 0.8, BETA = 0.86)
P <- as.data.frame(lapply(names(om), function(k) m$param[[k]] * exp(rnorm(40, 0, om[[k]]))))
names(P) <- names(om)
P$FB <- pmin(P$FB, 0.95); P$TTOCI <- 5
tt <- c(seq(0.1, 28, by = 0.3), seq(30, 365, by = 5))
a <- mrg_solve(m, P, times = tt, rtol = 1e-6, atol = 1e-8)
b <- mrg_solve_mrgsolve(m, P, times = tt)
report("cart transgene (log scale)", max(abs(log(a$CART) - log(b$CART))))

# PROTAC: in-vitro concentration sweep, and daily oral dosing for 2 weeks
m <- mrg_read("models/protac_degrader_kcat2023.cpp")
tt <- c(0, 6, 24, 72)
P <- data.frame(CFIX = 10^seq(-1, 5, by = 0.5))
a <- mrg_solve(m, P, times = tt, rtol = 1e-6, atol = 1e-9)
b <- mrg_solve_mrgsolve(m, P, times = tt)
report("protac degradation, in vitro", max(abs(a$DEG - b$DEG)))
ev <- data.frame(time = 0, cmt = "GUT", amt = 200, ii = 24, addl = 13)
tt <- seq(0.5, 400, by = 1)
a <- mrg_solve(m, events = ev, times = tt, rtol = 1e-6, atol = 1e-9)
b <- mrg_solve_mrgsolve(m, events = ev, times = tt)
report("protac degradation, oral dosing", max(abs(a$DEG - b$DEG)))
report("protac downstream response", max(abs(a$PDR - b$PDR)))

# T-cell engager: step-up dosing 1 / 2 / 60 / 60 mg (70 kg), outputs off the dose times
m <- mrg_read("models/tce_mosunetuzumab_hosseini2020.cpp")
days <- c(0, 7, 14, 21); mg <- c(1, 2, 60, 60)
ev <- rbind(data.frame(time = days, cmt = "TDBc_ugperkg", amt = mg * 1000 / 70),
            data.frame(time = days, cmt = "injection_effect", amt = 1),
            data.frame(time = days, cmt = "drug_effect", amt = 1))
tt <- seq(0.05, 35, by = 0.25)
a <- mrg_solve(m, events = ev, times = tt, rtol = 1e-5, atol = 1e-3)
b <- mrg_solve_mrgsolve(m, events = ev, times = tt)
for (v in c("TDBc_ugperml", "Bpb_perml", "totTpb_perml", "IL6combo", "Btiss_perml")) report(paste("tce", v), rel(a[[v]], b[[v]]))

# mRNA vaccine: the hand-vectorised right-hand side the browser uses must equal
# the translated model file, and the solution must match mrgsolve
m <- mrg_read("models/mrna_vaccine_dasti2025.cpp")
source(file.path("apps", "mrna-vaccine", "R", "mrna_rhs.R"))
source(file.path("apps", "mrna-vaccine", "R", "vaccines.R"))
P <- data.frame(TD2 = c(21, 28), TD3 = c(1e9, 150), expAg = c(3.37e5, 2e5))
mn <- mrg_main(m, P)
V <- c(lapply(stats::setNames(names(m$param), names(m$param)), function(nm) if (nm %in% names(P)) P[[nm]] else rep(m$param[[nm]], 2)), mn$vars)
V <- V[!duplicated(names(V), fromLast = TRUE)]
f_file <- mrg_rhs(m, names(V)); f_hand <- mrna_rhs(m)
set.seed(3)
d_rhs <- 0
for (k in 1:10) {
  Y <- mn$init * exp(rnorm(length(mn$init), 0, 2)) + matrix(rexp(length(mn$init), 1 / 50), nrow(mn$init))
  tk <- runif(1, 0, 60)
  a <- f_file(tk, Y, V); b <- f_hand(tk, Y, V)
  d_rhs <- max(d_rhs, max(abs(a - b) / (abs(a) + 1e-9 * max(abs(a)))))
}
report("mrna hand-vectorised vs model-file RHS", d_rhs)
s <- vaccine_setup("bnt", 30, 21)
tt <- c(seq(0.5, 20.5, by = 1), seq(21.5, 287, by = 3))
a <- mrg_solve(m, s$P, s$ev, tt, rtol = 1e-5, atol = 1e-6, nonneg = TRUE, rhs = mrna_rhs(m))
b <- mrg_solve_mrgsolve(m, s$P, s$ev, tt)
for (v in c("IGG", "LPTOT", "THELP")) report(paste("mrna", v), rel(a[[v]], b[[v]]))
presets <- readRDS(file.path("apps", "mrna-vaccine", "data", "presets.rds"))
pr <- presets[[preset_id("bnt", 30, 21)]]
b2 <- mrg_solve_mrgsolve(m, s$P, s$ev, OUT_TIMES)
report("mrna shipped preset vs mrgsolve (IgG)", rel(pr$IGG, b2$IGG[, 1]))

cat(sprintf("comparisons: %d, worst: %.2e\n", length(worst), max(worst)))
if (!all(is.finite(worst)) || max(worst) > 2e-3) stop("browser engine disagrees with mrgsolve")
cat("PASS\n")
