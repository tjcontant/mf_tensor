# Compares marginal predictive performance (absolute error, coverage, CRPS)
# across all emulators

# Setup -------------------------------------------------------------------

source("R/helper_fxns.R")
source("R/scoring_fxns.R")

load_latest("Simulation/Data", "^simulation_setup_.*\\.Rdata$")

# Produced by Simulation/R/generate_test_inputs.R -- run that first.
load("Simulation/Data/PostPreds/test_inputs.Rdata")

alpha_level <- 0.05
z <- stats::qnorm(1 - alpha_level / 2)

# in_shrinkwrap() is defined in R/helper_fxns.R (sourced above).

# Which test inputs to include below
idx_test_keep <- which(in_shrinkwrap(test_inputs, design[runs_lf, ]))

# Load each method's latest saved predictions --------------------------------

load_latest("Simulation/Data/Tensor", "^hf_tensor_.*\\.Rdata$")
load_latest("Simulation/Data/Tensor", "^lf_tensor_.*\\.Rdata$")
load_latest("Simulation/Data/Tensor", "^mf_tensor_.*\\.Rdata$")

methods <- list(
  `LF-TEN` = list(mean = pred_mean_lf, sd = pred_sd_lf),
  `HF-TEN` = list(mean = pred_mean_hf, sd = pred_sd_hf),
  `MF-TEN` = list(mean = pred_mean_mf, sd = pred_sd_mf)
)

# Compute per-(method, test input) averages -----------------------------------

results <- do.call(rbind, lapply(names(methods), function(method_name) {
  pred <- methods[[method_name]]
  mean_k <- pred$mean[, , , idx_test_keep, drop = FALSE]
  sd_k <- pred$sd[, , , idx_test_keep, drop = FALSE]
  truth_k <- truth$high_fi[, , , idx_test_keep, drop = FALSE]

  abs_err <- abs(mean_k - truth_k)
  covered <- abs_err <= z * sd_k
  crps_arr <- crps_gaussian(truth_k, mean_k, sd_k)

  data.frame(
    method = method_name,
    input = idx_test_keep,
    abs_error = apply(abs_err, 4, mean),
    coverage = apply(covered, 4, mean),
    crps = apply(crps_arr, 4, mean)
  )
}))

results$method <- factor(results$method, levels = names(methods))

# Save results ------------------------------------------------------------

dir.create("Simulation/Data/Compare", recursive = TRUE, showWarnings = FALSE)

save(
  results, alpha_level,
  file = paste0("Simulation/Data/Compare/marginal_metrics_", get_date(), ".Rdata")
)

# Clean up ------------------------------------------------------------------

rm(list = ls())
gc()
