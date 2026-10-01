# Setup -------------------------------------------------------------------

source("R/helper_fxns.R")
source("R/interpolation_fxns.R")
source("R/tensor_rank_fxns.R")
source("R/tensor_decomp_fxns.R")
source("R/tensor_gp_fxns.R")

load_latest("Simulation/Data", "^simulation_setup_.*\\.Rdata$")

# Produced by Simulation/R/generate_test_inputs.R -- run that first.
load("Simulation/Data/PostPreds/test_inputs.Rdata")

# Priors: half-normal on lambda_w/lambda_v, log-normal on rho_w/rho_v, gamma on
# lambda_eta/lambda_delta. The rho_v prior favors larger lengthscales given the
# sparse HF design.
prior <- list(
  sigma_lambda_w = 1, a_rho_w = 0, b_rho_w = 1,
  sigma_lambda_v = 1, a_rho_v = 0.75, b_rho_v = 0.3
)
a_eta <- 1
b_eta <- 0.5
a_delta <- 1
b_delta <- 0.5

set.seed(543)

design_train_lf <- design[runs_lf, ]
design_train_hf <- design[runs_hf, ]

# Set Tucker ranks ----------------------------------------------------------

ranks_lf <- select_rank_variance(
  tnsr = as.tensor(output_lf),
  target_var = 0.99
)

r_d_lf <- ranks_lf[length(ranks_lf)]

# LF Tucker decomposition -----------------------------------------------

# calculate_bases = TRUE gives prepare_sf_mcmc() its faster path
tuck_lf_mf <- apply_tucker_decomp(
  output = output_lf,
  ranks = ranks_lf,
  method = "hooi",
  calculate_bases = TRUE
)

# Prepare structures for LF MCMC -----------------------------------------

mcmc_prep_lf_mf <- prepare_sf_mcmc(
  output = output_lf,
  tucker_outputs = tuck_lf_mf,
  ranks = ranks_lf,
  a_eta = a_eta,
  b_eta = b_eta
)

# Interpolate LF bases onto the HF grid ----------------------------------

# Use same k as in generate_data.R
bases_interp <- interpolate_bases(tuck_lf_mf, coord_lf, coord_hf, k = 4)

# Discrepancy Tucker decomposition ---------------------------------------

# Discrepancy relative to the Tucker reconstruction of LF, interpolated to the
# HF grid
low_interp <- bases_interp$low_interp
discrep <- output_hf - low_interp[, , , runs_hf]

ranks_discrep <- select_rank_variance(as.tensor(discrep), target_var = 0.99)

r_d_discrep <- ranks_discrep[length(ranks_discrep)]

# calculate_bases = FALSE: the full basis matrices use a lot of memory
tuck_discrep <- apply_tucker_decomp(
  output = discrep,
  ranks = ranks_discrep,
  method = "hooi",
  calculate_bases = FALSE
)

# Prepare MF MCMC ---------------------------------------------------------

mcmc_prep_mf <- prepare_mf_mcmc(
  output_hf = output_hf,
  tucker_outputs_lf = tuck_lf_mf,
  tucker_outputs_discrep = tuck_discrep,
  interp_bases_outputs = bases_interp,
  ranks_lf = ranks_lf,
  ranks_discrep = ranks_discrep,
  a_delta = a_delta,
  b_delta = b_delta
)

# Run MCMC sampling -----------------------------------------------------

n_chains <- 4
n_iter_mf <- 10000
burn_in_mf <- 5000

run_one_chain_mf <- function(i, progress_file) {
  set.seed(1000 + i) # different seed for each chain

  init_lambda_eta_i <- mcmc_prep_lf_mf$a_eta_prime / mcmc_prep_lf_mf$b_eta_prime
  init_lambda_w_vec_i <- runif(r_d_lf, 0.5, 5)
  init_rho_w_mat_i <- matrix(runif(r_d_lf * ncol(design_train_lf), 0.5, 3),
    nrow = r_d_lf, ncol = ncol(design_train_lf)
  )

  init_lambda_delta_i <- mcmc_prep_mf$a_delta_prime / mcmc_prep_mf$b_delta_prime
  init_lambda_v_vec_i <- runif(r_d_discrep, 0.5, 5)
  init_rho_v_mat_i <- matrix(runif(r_d_discrep * ncol(design_train_hf), 0.5, 3),
    nrow = r_d_discrep, ncol = ncol(design_train_hf)
  )

  perform_mf_mcmc(
    mcmc_lf_setup_outputs = mcmc_prep_lf_mf,
    mcmc_mf_setup_outputs = mcmc_prep_mf,
    n_iter = n_iter_mf,
    burn_in = burn_in_mf,
    design_train_lf = design_train_lf,
    design_train_hf = design_train_hf,
    init_lambda_eta = init_lambda_eta_i,
    init_lambda_w_vec = init_lambda_w_vec_i,
    init_rho_w_mat = init_rho_w_mat_i,
    init_lambda_delta = init_lambda_delta_i,
    init_lambda_v_vec = init_lambda_v_vec_i,
    init_rho_v_mat = init_rho_v_mat_i,
    learning_rate = 125,
    prior = prior,
    progress_file = progress_file
  )
}

# Chains run in parallel (forked processes)
chains_mf <- run_chains_parallel(run_one_chain_mf, n_chains = n_chains, n_iter_total = n_iter_mf)

# Posterior predictions over held-out test inputs ----------------------------

n_post_draws <- 50

pred <- evaluate_mf_mcmc(
  test_inputs,
  tuck_lf_mf, tuck_discrep,
  bases_interp,
  mcmc_prep_lf_mf, mcmc_prep_mf,
  chains_mf,
  n_post_draws = n_post_draws,
  output_dim = dim(output_hf)[1:3],
  design_train_lf = design_train_lf,
  design_train_hf = design_train_hf,
  include_noise = TRUE,
  verbose = TRUE,
  r_d_lf = r_d_lf,
  r_d_discrep = r_d_discrep,
  D_G = NULL
)

D_G <- pred$D_G

pred_mean_mf <- pred$mean
pred_sd_mf <- pred$sd
draws_mf <- pred$draws

pred_mean_lf_mf <- pred$mean_lf
pred_sd_lf_mf <- pred$sd_lf
pred_mean_discrepancy_mf <- pred$mean_discrepancy
pred_sd_discrepancy_mf <- pred$sd_discrepancy

# Save outputs to file ----------------------------------------------------

dir.create("Simulation/Data/Tensor", recursive = TRUE, showWarnings = FALSE)

save(
  ranks_lf, ranks_discrep, tuck_lf_mf, mcmc_prep_lf_mf, tuck_discrep,
  bases_interp,
  mcmc_prep_mf, chains_mf,
  pred_mean_mf, pred_sd_mf, draws_mf,
  pred_mean_lf_mf, pred_sd_lf_mf, pred_mean_discrepancy_mf, pred_sd_discrepancy_mf,
  file = paste0("Simulation/Data/Tensor/mf_tensor_", get_date(), ".Rdata")
)

# Clean up ------------------------------------------------------------------

rm(list = ls())
gc()
