# Setup -------------------------------------------------------------------

source("R/helper_fxns.R")
source("R/interpolation_fxns.R")
source("R/tensor_rank_fxns.R")
source("R/tensor_decomp_fxns.R")
source("R/tensor_gp_fxns.R")

load_latest("Simulation/Data", "^simulation_setup_.*\\.Rdata$")

# Produced by Simulation/R/generate_test_inputs.R -- run that first.
load("Simulation/Data/PostPreds/test_inputs.Rdata")

# Priors: half-normal(sigma_lambda_w) on lambda_w, log-normal(a_rho_w, b_rho_w)
# on rho, gamma(a_eta, b_eta) on lambda_eta
prior <- list(sigma_lambda_w = 1, a_rho_w = 0, b_rho_w = 1)
a_eta <- 1
b_eta <- 0.5

set.seed(543)

design_train_hf <- design[runs_hf, ]

# Set Tucker ranks ----------------------------------------------------------

ranks_hf <- select_rank_variance(
  tnsr = as.tensor(output_hf),
  target_var = 0.99
)

r_d_hf <- ranks_hf[length(ranks_hf)]

# Tucker decomposition on training outputs -----------------------------------

tuck_hf <- apply_tucker_decomp(
  output = output_hf,
  ranks = ranks_hf,
  method = "hooi",
  calculate_bases = FALSE
)

# Prepare structures for MCMC -------------------------------------------

mcmc_prep_hf <- prepare_sf_mcmc(
  output = output_hf,
  tucker_outputs = tuck_hf,
  ranks = ranks_hf,
  a_eta = a_eta,
  b_eta = b_eta
)

# Run MCMC sampling -----------------------------------------------------

n_chains <- 4 # many chains helps diagnose convergence of multivariate chains (e.g. rho_mat)
n_iter_hf <- 10000
burn_in_hf <- 5000

run_one_chain_hf <- function(i, progress_file) {
  set.seed(1000 + i) # different seed for each chain

  init_lambda_eta_i <- mcmc_prep_hf$a_eta_prime / mcmc_prep_hf$b_eta_prime
  init_lambda_w_vec_i <- runif(r_d_hf, 0.5, 5)
  init_rho_mat_i <- matrix(runif(r_d_hf * ncol(design_train_hf), 0.5, 3),
    nrow = r_d_hf, ncol = ncol(design_train_hf)
  )

  perform_sf_mcmc(
    output = output_hf,
    design = design_train_hf,
    mcmc_setup_outputs = mcmc_prep_hf,
    n_iter = n_iter_hf,
    burn_in = burn_in_hf,
    r3 = r_d_hf,
    init_lambda_eta = init_lambda_eta_i,
    init_lambda_w_vec = init_lambda_w_vec_i,
    init_rho_mat = init_rho_mat_i,
    learning_rate = 50, # updating jump sizes during burn_in
    prior = prior,
    progress_file = progress_file
  )
}

# Chains run in parallel (forked processes)
chains_hf <- run_chains_parallel(run_one_chain_hf, n_chains = n_chains, n_iter_total = n_iter_hf)

# Posterior predictions over held-out test inputs ----------------------------

n_post_draws <- 50

# Predict at all test inputs
pred <- evaluate_sf_mcmc(
  untested_input = test_inputs,
  mcmc_setup_outputs = mcmc_prep_hf,
  mcmc_outputs = chains_hf,
  tucker_outputs = tuck_hf,
  n_post_draws = n_post_draws,
  design_train = design_train_hf,
  output_dim = dim(output_hf)[1:3],
  ranks = ranks_hf,
  hf = TRUE,
  include_noise = TRUE,
  verbose = TRUE
)

pred_mean_hf <- pred$mean
pred_sd_hf <- pred$sd
draws_hf <- pred$draws

# Save outputs to file --------------------------------------------------

dir.create("Simulation/Data/Tensor", recursive = TRUE, showWarnings = FALSE)

save(
  ranks_hf, tuck_hf, mcmc_prep_hf, chains_hf,
  pred_mean_hf, pred_sd_hf, draws_hf,
  file = paste0("Simulation/Data/Tensor/hf_tensor_", get_date(), ".Rdata")
)

# Clean up ----------------------------------------------------------------

rm(list = ls())
gc()
