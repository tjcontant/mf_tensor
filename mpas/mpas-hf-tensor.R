setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved setup data
load("save/mpas-data.Rdata")


# Define prior hyperparameters for high-fidelity emulator
sigma_lambda_w <- 1
a_rho_w <- 0
b_rho_w <- 1
a_eta <- 1
b_eta <- 0.5

# Clear console
cat("\014")

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

set.seed(543)


# ------------------------------------------------------------------------------
# Begin Leave-One-Out Cross-Validation loop for high-fidelity runs
# ------------------------------------------------------------------------------

for (loo in runs_hf) {
  message <- paste0("HF: ", which(loo == runs_hf), " / ", length(runs_hf))
  pushover_quiet(message)
  
  loo_hf_idx <- which(loo == runs_hf)
  
  # Prepare training data by removing current LOO index
  d_train_hf <- d_hf - 1
  runs_train_hf <- runs_hf[-loo_hf_idx]
  design_train_hf <- design_unif[runs_train_hf, ]
  
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
  
  ranks_hf <- pmax(ranks_lf, ranks_hf)
  
  r_s_hf <- ranks_hf[1]
  r_m_hf <- ranks_hf[2]
  r_y_hf <- ranks_hf[3]
  r_d_hf <- ranks_hf[4]
  
  cat("Ranks:\t")
  cat(ranks_hf)
  cat("\t(Target 99% explained var)\n\n")
  
  end_rank_select <- Sys.time()
  
  timings[["rank_select"]] <- as.numeric(difftime(start_rank_select, end_rank_select, units = "secs"))
  
  
  # --------------------------------------
  # Tucker decomposition on training outputs
  # --------------------------------------
  cat(paste0("\tTucker\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_tucker <- Sys.time()
  
  tuck_hf <- apply_tucker_decomp(output = output_hf_4d_trans_scaled[, , , -loo_hf_idx], 
                                 ranks = ranks_hf, 
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
  
  mcmc_prep_hf <- prepare_sf_mcmc(output = output_hf_4d_trans_scaled[, , , -loo_hf_idx],
                                  tucker_outputs = tuck_hf,
                                  ranks = ranks_hf)
  
  end_mcmc_prep <- Sys.time()
  
  timings[["mcmc_prep"]] <- as.numeric(difftime(start_mcmc_prep, end_mcmc_prep, units = "secs"))
  
  
  # --------------------------------------
  # Run MCMC sampling
  # --------------------------------------
  cat(paste0("\tMCMC\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc <- Sys.time()
  
  n_chains <- 3 # as there are many chains, good for diagnosing convergence of multivariate chains (e.g. rho_mat)
  chains_hf <- vector("list", n_chains)
  for (i in 1:n_chains) {
    set.seed(1000 + i)  # different seed for each chain
    
    # Generate different initial values for each chain
    init_lambda_eta_i <- mcmc_prep_hf$a_eta_prime / mcmc_prep_hf$b_eta_prime
    init_lambda_w_vec_i <- runif(r_d_hf, 0.5, 5)
    init_rho_mat_i <- matrix(runif(r_d_hf * ncol(design_train_hf), 0, 1), nrow = r_d_hf, ncol = ncol(design_train_hf))
    
    chains_hf[[i]] <- perform_sf_mcmc(
      output = output_hf_4d_trans_scaled[, , , -loo_hf_idx],
      design = design_train_hf,
      mcmc_setup_outputs = mcmc_prep_hf,
      n_iter = 4000,
      burn_in = 2000,
      r3 = r_d_hf,
      init_lambda_eta = init_lambda_eta_i,
      init_lambda_w_vec = init_lambda_w_vec_i,
      init_rho_mat = init_rho_mat_i,
      learning_rate = 150, # updating jump sizes during burn_in
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
  retain_chains_hf <- check_geweke_multiple_chains(chains_hf, frac1 = 0.25, frac2 = 0.5)
  
  # Gelman-Reuben (between chain)
  convergence_results <- check_mcmc_convergence(retain_chains_hf)
  check_rhat(convergence_results$lambda_eta, name = "lambda_eta")
  
  check_rhat(convergence_results$lambda_w_vec, name = "lambda_w_vec")
  check_rhat(convergence_results$rho_mat, name = "rho_mat")
  
  check_rhat(convergence_results$log_lik, name = "log_lik")
  
  # Global PSRF
  check_mpsrf(convergence_results$lambda_w_vec, name = "lambda_w_vec", threshold = 1.2)  # Multivariate R̂
  check_mpsrf(convergence_results$rho_mat, name = "rho_mat", threshold = 1.2)  # Multivariate R̂
  
  # Rename for clarity
  chains_hf <- retain_chains_hf
  
  
  # --------------------------------------
  # Evaluate posterior predictive distribution
  # --------------------------------------
  cat(paste0("\tPost Preds\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_post_preds <- Sys.time()
  
  post_preds_hf_noise <- evaluate_sf_mcmc(
    untested_input = design_unif[loo, ],
    mcmc_setup_outputs = mcmc_prep_hf,
    mcmc_outputs = chains_hf,
    tucker_outputs = tuck_hf,
    n_post_draws = 250,
    design_train = design_train_hf,
    output_dim = dim(output_hf_4d_trans)[1:3],
    ranks = ranks_hf,
    hf = T,
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
  file_save <- paste0("save/maps-hf-tensor/", loo, ".Rdata")
  save(ranks_hf, tuck_hf, mcmc_prep_hf, chains_hf, 
       # post_preds_hf, cv_stats_hf,
       post_preds_hf_noise, 
       # cv_stats_hf_noise, 
       timings,
       file = file_save)
}

rm(list = ls())
gc()