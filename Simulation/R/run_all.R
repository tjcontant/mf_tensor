# Runs the full simulation study end to end, from the repository root.

# Generate data -----------------------------------------------------------

source("Simulation/R/generate_data.R")
source("Simulation/R/generate_test_inputs.R")

# Emulators ---------------------------------------------------------------

source("Simulation/R/hf_tensor.R")
source("Simulation/R/lf_tensor.R")
source("Simulation/R/mf_tensor.R")

# Compare methods -------------------------------------------------------

source("Simulation/R/compare_marginal_metrics.R")
source("Simulation/R/compare_variogram_score.R")
