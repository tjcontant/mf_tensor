setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved sea ice data
load("save/mpas-data.Rdata")

# Define prior hyperparameters for high-fidelity emulator
sigma_lambda_w <- 1
a_rho_w <- 0
b_rho_w <- 1
a_eta <- 1
b_eta <- 0.5

# Remove no-ice locations
no_ice_idx_lf <- apply(output_lf_4d_unscaled, 1, function(x) all(x == 0))
output_lf_4d_unscaled <- output_lf_4d_unscaled[!no_ice_idx_lf, , , ]
output_lf_4d_trunc <- output_lf_4d_trunc[!no_ice_idx_lf, , , ]
output_lf_4d_trans <- output_lf_4d_trans[!no_ice_idx_lf, , , ]
coord_lf <- coord_lf[!no_ice_idx_lf, ]
s_lf <- nrow(coord_lf)
n_lf <- s_lf * t_lf

no_ice_idx_hf <- apply(output_hf_4d_unscaled, 1, function(x) all(x == 0))
output_hf_4d_unscaled <- output_hf_4d_unscaled[!no_ice_idx_hf, , , ]
output_hf_4d_trunc <- output_hf_4d_trunc[!no_ice_idx_hf, , , ]
output_hf_4d_trans <- output_hf_4d_trans[!no_ice_idx_hf, , , ]
coord_hf <- coord_hf[!no_ice_idx_hf, ]
s_hf <- nrow(coord_hf)
n_hf <- s_hf * t_hf


output_lf_4d_trans_scaled <- scale_transformed(output_lf_4d_trans, max(output_hf_4d_trans))
output_hf_4d_trans_scaled <- scale_transformed(output_hf_4d_trans, max(output_hf_4d_trans))

output_interp_4d_trans_scaled <- scale_transformed(output_interp_4d_trans[!no_ice_idx_hf, , , ], max(output_hf_4d_trans))

discrep_4d_trans_scaled <- output_hf_4d_trans_scaled - output_interp_4d_trans_scaled


# Clear console
cat("\014")

set.seed(543)


# ------------------------------------------------------------------------------
# Begin Leave-One-Out Cross-Validation loop for high-fidelity runs
# ------------------------------------------------------------------------------

for (loo in runs_hf) {
  message <- paste0("LF: ", which(loo == runs_hf), " / ", length(runs_hf))
  pushover_quiet(message)
  
  loo_hf_idx <- which(loo == runs_hf)
  
  # Prepare training data by removing current LOO index
  d_train_lf <- d_lf - 1
  runs_train_lf <- runs_lf[-loo]
  design_train_lf <- design_unif[runs_train_lf, ]
  
  timings <- list()
  
  cat(paste0("\nRun ", loo, "\n"))
  
  
  # --------------------------------------
  # Set Tucker ranks
  # --------------------------------------  
  cat(paste0("\tRank Select.\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_rank_select <- Sys.time()
  
  ranks_lf <- select_rank_variance(
    tnsr = as.tensor(output_lf_4d_trans_scaled[, , , -loo]),
    target_var = 0.99
  )
  
  ranks_hf <- select_rank_variance(
    tnsr = as.tensor(output_hf_4d_trans_scaled[, , , -loo_hf_idx]),
    target_var = 0.99
  )
  
  ranks_lf <- pmax(ranks_lf, ranks_hf)
  
  r_s_lf <- ranks_lf[1]
  r_m_lf <- ranks_lf[2]
  r_y_lf <- ranks_lf[3]
  r_d_lf <- ranks_lf[4]
  
  
  cat("Ranks:\t")
  cat(ranks_lf)
  cat("\t(Target 99% explained var)\n\n")
  
  end_rank_select <- Sys.time()
  
  timings[["rank_select"]] <- as.numeric(difftime(start_rank_select, end_rank_select, units = "secs"))
  
  
  # --------------------------------------
  # Tucker decomposition on training outputs
  # --------------------------------------
  cat(paste0("\tTucker\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_tucker <- Sys.time()
  
  tuck_lf <- apply_tucker_decomp(output = output_lf_4d_trans_scaled[ , , , -loo], 
                                 ranks = ranks_lf, 
                                 threshold = 95,
                                 method = "hooi",
                                 calculate_bases = F)
  
  end_tucker <- Sys.time()
  
  timings[["tucker"]] <- as.numeric(difftime(start_tucker, end_tucker, units = "secs"))
  
  
  # --------------------------------------
  # Prepare structures for MCMC
  # --------------------------------------
  cat(paste0("\tMCMC Prep\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc_prep <- Sys.time()
  
  mcmc_prep_lf <- prepare_sf_mcmc(output = output_lf_4d_trans_scaled[ , , , -loo],
                                  tucker_outputs = tuck_lf,
                                  ranks = ranks_lf)
  
  print(mcmc_prep_lf$a_eta_prime / mcmc_prep_lf$b_eta_prime)
  
  end_mcmc_prep <- Sys.time()
  
  timings[["mcmc_prep"]] <- as.numeric(difftime(start_mcmc_prep, end_mcmc_prep, units = "secs"))
  
  
  # --------------------------------------
  # Run MCMC sampling
  # --------------------------------------
  cat(paste0("\tMCMC\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc <- Sys.time()
  
  n_chains <- 3 # as there are many chains, good for diagnosing convergence of multivariate chains (e.g. rho_mat)
  chains_lf <- vector("list", n_chains)
  for (i in 1:n_chains) {
    set.seed(1000 + i)  # different seed for each chain
    
    # Generate different initial values for each chain
    init_lambda_eta_i <- mcmc_prep_lf$a_eta_prime / mcmc_prep_lf$b_eta_prime
    init_lambda_w_vec_i <- runif(r_d_lf, 0.5, 5)
    init_rho_mat_i <- matrix(runif(r_d_lf * ncol(design_train_lf), 0.5, 5), nrow = r_d_lf, ncol = ncol(design_train_lf))
    
    chains_lf[[i]] <- perform_sf_mcmc(
      output = output_lf_4d_trans_scaled[, , , -loo],
      design = design_train_lf,
      mcmc_setup_outputs = mcmc_prep_lf,
      n_iter = 4000,
      burn_in = 2000,
      r3 = r_d_lf,
      init_lambda_eta = init_lambda_eta_i,
      init_lambda_w_vec = init_lambda_w_vec_i,
      init_rho_mat = init_rho_mat_i,
      learning_rate = 75, # updating jump sizes during burn_in
      verbose = F
    )
    
    cat(paste0("\t  Chain ", i, "\t", format(Sys.time(), "%b %d %X"), "\n"))
  }
  
  end_mcmc <- Sys.time()
  
  timings[["mcmc"]] <- as.numeric(difftime(start_mcmc, end_mcmc, units = "secs"))
  
  
  # --------------------------------------
  # Check MCMC convergence
  # --------------------------------------  
  cat(paste0("\tCheck MCMC\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  # Geweke (within chain)
  retain_chains_lf <- check_geweke_multiple_chains(chains_lf, frac1 = 0.25, frac2 = 0.5)
  
  # Gelman-Reuben (between chain)
  convergence_results <- check_mcmc_convergence(retain_chains_lf)
  check_rhat(convergence_results$lambda_eta, name = "lambda_eta")
  
  check_rhat(convergence_results$lambda_w_vec, name = "lambda_w_vec")
  check_rhat(convergence_results$rho_mat, name = "rho_mat")
  
  check_rhat(convergence_results$log_lik, name = "log_lik")
  
  # Global PSRF
  check_mpsrf(convergence_results$lambda_w_vec, name = "lambda_w_vec", threshold = 1.2)  # Multivariate R̂
  check_mpsrf(convergence_results$rho_mat, name = "rho_mat", threshold = 1.2)  # Multivariate R̂
  
  # Rename for clarity
  chains_lf <- retain_chains_lf
  
  
  # --------------------------------------
  # Evaluate posterior predictive distribution
  # --------------------------------------
  cat(paste0("\tPost Preds\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_post_preds <- Sys.time()
  
  post_preds_lf_noise <- evaluate_sf_mcmc(
    untested_input = design_unif[loo, ],
    mcmc_setup_outputs = mcmc_prep_lf,
    mcmc_outputs = chains_lf,
    tucker_outputs = tuck_lf,
    n_post_draws = 250,
    design_train = design_train_lf,
    output_dim = dim(output_lf_4d_trans_scaled)[1:3],
    ranks = ranks_lf,
    hf = F,
    transform = "logistic-soft-clip",
    include_noise = TRUE,
    obs = output_hf_4d_unscaled[, , , loo_hf_idx]
  )
  
  end_post_preds <- Sys.time()
  
  timings[["post_preds"]] <- as.numeric(difftime(start_post_preds, end_post_preds, units = "secs"))
  
  # --------------------------------------
  # Save outputs to file
  # --------------------------------------
  cat(paste0("\tSaving\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  file_save <- paste0("save/mpas-lf-tensor/", loo, ".Rdata")
  save(ranks_lf, tuck_lf, mcmc_prep_lf, chains_lf, 
       post_preds_lf_noise, 
       timings,
       file = file_save)
}

rm(list = ls())
gc()