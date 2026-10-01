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

design_train_lf <- design[runs_lf, ]

# Set Tucker ranks ----------------------------------------------------------

ranks_lf <- select_rank_variance(
  tnsr = as.tensor(output_lf),
  target_var = 0.99
)

r_d_lf <- ranks_lf[length(ranks_lf)]

# Tucker decomposition on training outputs -----------------------------------

tuck_lf <- apply_tucker_decomp(
  output = output_lf,
  ranks = ranks_lf,
  method = "hooi",
  calculate_bases = FALSE
)

# Prepare structures for MCMC -------------------------------------------

mcmc_prep_lf <- prepare_sf_mcmc(
  output = output_lf,
  tucker_outputs = tuck_lf,
  ranks = ranks_lf,
  a_eta = a_eta,
  b_eta = b_eta
)

# Run MCMC sampling -----------------------------------------------------

n_chains <- 4 # many chains helps diagnose convergence of multivariate chains (e.g. rho_mat)
n_iter_lf <- 10000
burn_in_lf <- 5000

run_one_chain_lf <- function(i, progress_file) {
  set.seed(1000 + i) # different seed for each chain

  init_lambda_eta_i <- mcmc_prep_lf$a_eta_prime / mcmc_prep_lf$b_eta_prime
  init_lambda_w_vec_i <- runif(r_d_lf, 0.5, 5)
  init_rho_mat_i <- matrix(runif(r_d_lf * ncol(design_train_lf), 0.5, 3),
    nrow = r_d_lf, ncol = ncol(design_train_lf)
  )

  perform_sf_mcmc(
    output = output_lf,
    design = design_train_lf,
    mcmc_setup_outputs = mcmc_prep_lf,
    n_iter = n_iter_lf,
    burn_in = burn_in_lf,
    r3 = r_d_lf,
    init_lambda_eta = init_lambda_eta_i,
    init_lambda_w_vec = init_lambda_w_vec_i,
    init_rho_mat = init_rho_mat_i,
    learning_rate = 50, # updating jump sizes during burn_in
    prior = prior,
    progress_file = progress_file
  )
}

# Chains run in parallel (forked processes)
chains_lf <- run_chains_parallel(run_one_chain_lf, n_chains = n_chains, n_iter_total = n_iter_lf)

# Posterior predictions over held-out test inputs ----------------------------

n_post_draws <- 50

# Predict at all test inputs
pred <- evaluate_sf_mcmc(
  untested_input = test_inputs,
  mcmc_setup_outputs = mcmc_prep_lf,
  mcmc_outputs = chains_lf,
  tucker_outputs = tuck_lf,
  n_post_draws = n_post_draws,
  design_train = design_train_lf,
  output_dim = dim(output_lf)[1:3],
  ranks = ranks_lf,
  hf = FALSE,
  coord_lf = coord_lf,
  coord_hf = coord_hf,
  include_noise = TRUE,
  verbose = TRUE
)

pred_mean_lf <- pred$mean
pred_sd_lf <- pred$sd
draws_lf <- pred$draws

# Save outputs to file --------------------------------------------------

dir.create("Simulation/Data/Tensor", recursive = TRUE, showWarnings = FALSE)

save(
  ranks_lf, tuck_lf, mcmc_prep_lf, chains_lf,
  pred_mean_lf, pred_sd_lf, draws_lf,
  file = paste0("Simulation/Data/Tensor/lf_tensor_", get_date(), ".Rdata")
)

# Clean up ----------------------------------------------------------------

rm(list = ls())
gc()
