# Setup -------------------------------------------------------------------

set.seed(543)

rm(list = ls())
gc()

source("R/helper_fxns.R")

load_latest("Simulation/Data", "^simulation_setup_.*\\.Rdata$")

# Generate test inputs --------------------------------------------------

n_inputs <- 100

# Test inputs at least min_dist from every design run and from each other
test_inputs <- generate_test_inputs(n_inputs, design, lower = 0.05, upper = 0.95, min_dist = 0.05, weights = c(1, 1, 0.5))

# Ground truth ------------------------------------------------------------

# Simulated once and shared by every method
truth <- simulate_for_design(
  test_inputs, coord_lf, coord_hf,
  year_ids, month_ids,
  start_lf = 1, end_lf = 5, k_interp = 4
)

# Save results ------------------------------------------------------------

dir.create("Simulation/Data/PostPreds", recursive = TRUE, showWarnings = FALSE)

save(n_inputs, test_inputs, truth, file = "Simulation/Data/PostPreds/test_inputs.Rdata")

# Clean up ------------------------------------------------------------------

rm(list = ls())
gc()
