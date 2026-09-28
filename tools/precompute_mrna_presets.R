# Precompute the mRNA vaccine app's standard schedules with mrgsolve, so the
# app shows them at once instead of waiting for the browser engine. Other
# settings are simulated live. tests/test_engine_vs_mrgsolve.R checks the
# stored results against a fresh browser-engine solve.
#
# Run from the repository root: Rscript tools/precompute_mrna_presets.R
#
# Part of a derivative of the COSBI "Multiscale QSP model for mRNA vaccines";
# COSBI-SSLA licence (non-commercial), see apps/mrna-vaccine/LICENSE.

source("shared/ode_engine.R")
source("reference/mrgsolve_engine.R")
source("apps/mrna-vaccine/R/vaccines.R")
m <- mrg_read("models/mrna_vaccine_dasti2025.cpp")
out <- list()
for (run in PRESET_RUNS) {
  s <- vaccine_setup(run$key, run$dose, run$interval)
  r <- mrg_solve_mrgsolve(m, s$P, s$ev, OUT_TIMES)
  out[[preset_id(run$key, run$dose, run$interval)]] <- vaccine_outputs(r)
  cat(preset_id(run$key, run$dose, run$interval), " peak IgG", signif(max(r$IGG), 4), "ng/mL\n")
}
saveRDS(out, "apps/mrna-vaccine/data/presets.rds", version = 2)
