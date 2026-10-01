# Compares SPATIOTEMPORAL correlation-awareness across all emulators
# via the Variogram Score (Scheuerer & Hamill 2015)

# Setup -------------------------------------------------------------------

library(pbmcapply)

source("R/helper_fxns.R")
source("R/scoring_fxns.R")

n_workers <- 4
n_pairs <- 10000

# Short-lag pair-sampling scale, shared between sample_pairs_restricted() and
# spatiotemporal_pair_weights() so both treat "nearby" the same way
max_spatial_dist <- 1.1
max_temporal_lag <- 1

load_latest("Simulation/Data", "^simulation_setup_.*\\.Rdata$")
load("Simulation/Data/PostPreds/test_inputs.Rdata")

load_latest("Simulation/Data/Tensor", "^lf_tensor_.*\\.Rdata$")
load_latest("Simulation/Data/Tensor", "^hf_tensor_.*\\.Rdata$")
load_latest("Simulation/Data/Tensor", "^mf_tensor_.*\\.Rdata$")

# Every field is (space, month, year) -- s_hf * n_months * n_years locations.
n_loc <- prod(dim(truth$high_fi)[1:3])

n_draws <- dim(draws_mf)[4]

# Draw-generating function per method --------------------------------------

# Each function returns the FULL [n_draws x n_loc] spatiotemporal field for one
# test input `j`, flattened in the same (space fastest, then month, then year)
# order as as.vector(truth$high_fi[, , , j])
methods <- list(
  `LF-TEN` = function(j) t(matrix(draws_lf[, , , , j], nrow = n_loc, ncol = n_draws)),
  `HF-TEN` = function(j) t(matrix(draws_hf[, , , , j], nrow = n_loc, ncol = n_draws)),
  `MF-TEN` = function(j) t(matrix(draws_mf[, , , , j], nrow = n_loc, ncol = n_draws))
)

# Compute Variogram Score per (method, test input) ---------------------------

# One test input is a single forked unit of work
# Seeded per test input so pair sampling is reproducible across forked workers
vs_per_input <- pbmcapply::pbmclapply(seq_len(n_inputs), function(j) {
  set.seed(2718 + j)
  truth_field <- as.vector(truth$high_fi[, , , j])
  pairs <- sample_pairs_restricted(
    coord_hf,
    n_time = n_months * n_years, n_pairs,
    max_spatial_dist = max_spatial_dist, max_temporal_lag = max_temporal_lag
  )
  weights <- spatiotemporal_pair_weights(
    coord_hf,
    n_time = n_months * n_years, pairs$i, pairs$j,
    spatial_scale = max_spatial_dist, temporal_scale = max_temporal_lag
  )

  sapply(names(methods), function(method_name) {
    draws_field <- methods[[method_name]](j)
    variogram_score_subsampled(draws_field, truth_field, pairs$i, pairs$j, weights = weights)
  })
}, mc.cores = n_workers)

# Reassemble into one [n_inputs x n_methods] matrix.
vs_mat <- do.call(rbind, vs_per_input)

# Collect per-(method, test input) scores -----------------------------------

results <- do.call(rbind, lapply(names(methods), function(method_name) {
  data.frame(
    method = method_name,
    input  = seq_len(n_inputs),
    vs     = vs_mat[, method_name]
  )
}))
results$method <- factor(results$method, levels = names(methods))

# Save results ------------------------------------------------------------

dir.create("Simulation/Data/Compare", recursive = TRUE, showWarnings = FALSE)

save(
  vs_mat, results, n_pairs, n_loc, max_spatial_dist, max_temporal_lag,
  file = paste0("Simulation/Data/Compare/variogram_score_", get_date(), ".Rdata")
)

# Clean up ------------------------------------------------------------------

rm(list = ls())
gc()
