# QSP mAb TMDD Simulator

A quantitative systems pharmacology (QSP) model of a monoclonal antibody that
binds a membrane target, with an interactive Shiny dashboard. It follows the
chain behind most antibody dose selection, from dose to exposure, exposure to
target engagement, and target engagement to a downstream biomarker. The main
question it answers is which dose and interval keep receptor occupancy above
target at trough across a population.

The model is written for [mrgsolve](https://mrgsolve.org) (`models/tmdd.cpp`).
The browser version integrates the same equations with an adaptive
Dormand-Prince solver in base R, and the test suite checks the two against each
other before every deploy.

This is for research and teaching only. The defaults are round numbers typical
of an IgG1, not estimates for any real antibody or target, and nothing here is
validated for clinical use.

**Browser version:** <https://wrlog.github.io/QSP/> (runs entirely in the
browser through WebAssembly, so there's nothing to install).

## Model

- **Antibody PK.** Two compartments with linear clearance. Doses are IV
  bolus, or subcutaneous with first-order absorption and bioavailability.
- **Target.** A membrane target synthesised at kdeg × R0, degraded at kdeg
  when free, and internalised with the drug at kint when bound. If kint > kdeg,
  total target falls under treatment; if kint < kdeg, it accumulates.
- **Binding.** Quasi-steady-state approximation (Gibiansky et al., J
  Pharmacokinet Pharmacodyn 2008;35:573). Free drug comes from total drug and
  total target through the QSS quadratic, written in a form that stays
  accurate when almost all the drug is bound.
- **Occupancy.** Complex / total target. At quasi-steady state this equals
  C / (Kss + C), so 90% occupancy needs free drug at 9 × Kss.
- **Biomarker.** An indirect response to free target. Production scales as
  1 − Imax + Imax × (Rfree/R0)^γ, so full target suppression lowers the
  biomarker by at most Imax.
- **Population.** Weight effects on CL, Q (exponent 0.8) and V1, V2 (0.6), and
  log-normal variability on CL, V1, R0 and ka.

| Parameter | Default | |
|---|---|---|
| CL, V1, Q, V2 | 0.2 L/day, 3 L, 0.5 L/day, 2.5 L | at 70 kg |
| ka, F (SC) | 0.25 /day, 0.7 | |
| Kss | 1 nM | |
| R0, kdeg, kint | 2 nM, 0.2 /day, 1 /day | |
| kout, Imax, γ | 0.3 /day, 0.8, 1 | biomarker |

## Dashboard

| Tab | What's on it |
|---|---|
| Simulation | Free drug, receptor occupancy and biomarker over time (median with 50%/90% prediction intervals); trough occupancy, share of subjects at target, trough concentration and biomarker change |
| Dose selection | Minimum occupancy over the final interval, and the share of subjects at target, across a 0.01–10 mg/kg dose grid for the current route and schedule |
| TMDD mechanism | Dose-normalised profiles after a single dose, which show the TMDD bend, and free and total target under the current regimen |
| Model setup | Every PK, target, binding, biomarker and variability parameter |

The sidebar holds the population weight range, route, dose (mg/kg or flat), an
optional loading dose, interval, number of doses and the occupancy target.

## Two engines

mrgsolve compiles C++, which a browser can't do. The model is nonlinear, so
there's no closed-form shortcut. Instead, `R/tmdd_engine.R` implements an
adaptive Dormand-Prince 5(4) integrator in base R. It's vectorised across
subjects, so a population, or a whole dose sweep, is integrated as one system
with a step size that satisfies every subject's error tolerance. A 7-dose × 200
subject sweep takes about a second in desktop R.

Run locally with mrgsolve installed, the app simulates through mrgsolve. The
About tab shows which engine is in use. Individual parameters are computed in R
(`R/population.R`) and passed to both engines, so they get identical inputs.

`tests/test_engine_vs_mrgsolve.R` compares the engines across three target
scenarios (the defaults, an accumulating target, and a high-affinity,
high-capacity target), both routes, three dose levels and a loading dose. It
fails if they differ by more than 1e-4. The current worst case is about 3e-7.

`tests/test_model.R` checks the QSS algebra, baseline stability, the
two-compartment closed form when there's no target, saturation of
target-mediated clearance, loading doses, the weight exponents and the
dose-response ordering.

## Running it locally

```r
install.packages(c("shiny", "shinydashboard", "DT", "ggplot2", "scales"))
install.packages("mrgsolve")   # optional; needs a C++ toolchain (Rtools on Windows)
shiny::runApp()
```

Tests, from the repository root:

```sh
Rscript tests/test_model.R                # base R only
Rscript tests/test_engine_vs_mrgsolve.R   # needs mrgsolve
```

## Using it from the console

```r
source("R/tmdd_engine.R"); source("R/population.R")

P <- tmdd_individual(200, wt_range = c(50, 100),
                     cv = list(CL = 30, V1 = 20, R0 = 30, KA = 20), seed = 1)

# 0.3 mg/kg SC every 2 weeks, 6 doses
reg <- list(route = "sc", amt = 0.3 * P$WT, interval = 14, n_doses = 6)
res <- tmdd_simulate(P, reg, times = seq(0, 84, by = 0.5))
fi  <- final_interval(res, reg)
mean(fi$min_ro >= 0.9)            # share of subjects at >= 90% occupancy

# Dose sweep
dose_ranging(c(0.03, 0.1, 0.3, 1), P, reg, seq(0, 84, by = 0.5), target_ro = 0.9)

# The same simulation through mrgsolve
source("reference/mrgsolve_engine.R")
res_mrg <- tmdd_simulate_mrgsolve(P, reg, seq(0, 84, by = 0.5))
```

## Files

- `app.R`: the Shiny dashboard
- `models/tmdd.cpp`: the mrgsolve model
- `R/tmdd_engine.R`: QSS equations and the vectorised Dormand-Prince solver
- `R/population.R`: default parameters, virtual populations, response summaries and dose ranging
- `R/theme.R`: styling
- `reference/mrgsolve_engine.R`: runs `models/tmdd.cpp` through the same interface
- `tests/`: model properties and the engine-vs-mrgsolve check
- `.github/workflows/shinylive.yml`: tests, WebAssembly export and GitHub Pages deploy

## License

MIT
