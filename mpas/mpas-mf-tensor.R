setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved setup data
load("save/mpas-data.Rdata")


# Define prior hyperparameters for LF emulator
sigma_lambda_w <- 1
a_rho_w <- 0
b_rho_w <- 1
a_eta <- 1
b_eta <- 0.5

# Define prior hyperparameters for discrep emulator
sigma_lambda_v <- 1
a_rho_v <- 0
b_rho_v <- 1
a_delta <- 1
b_delta <- 0.5

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

discrep_4d_unscaled <- discrep_4d_unscaled[!no_ice_idx_hf, , , ]
discrep_4d_trans <- discrep_4d_trans[!no_ice_idx_hf, , , ]

output_lf_4d_trans_scaled <- scale_transformed(output_lf_4d_trans, max(output_hf_4d_trans))
output_hf_4d_trans_scaled <- scale_transformed(output_hf_4d_trans, max(output_hf_4d_trans))

output_interp_4d_trans_scaled <- scale_transformed(output_interp_4d_trans[!no_ice_idx_hf, , , ], max(output_hf_4d_trans))

discrep_4d_trans_scaled <- output_hf_4d_trans_scaled - output_interp_4d_trans_scaled

set.seed(543)

# Clear console
cat("\014")


# ------------------------------------------------------------------------------
# Begin Leave-One-Out Cross-Validation loop for high-fidelity runs
# ------------------------------------------------------------------------------

for (loo in runs_hf) {
  message <- paste0("MF: ", which(loo == runs_hf), " / ", length(runs_hf))
  pushover_quiet(message)
  
  loo_hf_idx <- which(loo == runs_hf)
  
  # Prepare training data by removing current LOO index
  d_train_lf <- d_lf - 1
  runs_train_lf <- runs_lf[-loo]
  design_train_lf <- design_unif[runs_train_lf, ]
  
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
  
  ranks_lf <- pmax(ranks_lf, ranks_hf)
  
  r_s_lf <- ranks_lf[1]
  r_m_lf <- ranks_lf[2]
  r_y_lf <- ranks_lf[3]
  r_d_lf <- ranks_lf[4]
  
  cat("\t\tRanks LF:\t")
  cat(ranks_lf)
  cat("\t(Target 98.5% explained var)\n")
  
  end_rank_select <- Sys.time()
  
  timings[["rank_select"]] <- as.numeric(difftime(start_rank_select, end_rank_select, units = "secs"))
  
  
  # --------------------------------------
  # LF Tucker decomposition
  # --------------------------------------
  cat(paste0("\tLF Tucker\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_tucker_lf <- Sys.time()
  
  tuck_lf_mf <- apply_tucker_decomp(output = output_lf_4d_trans_scaled[ , , , -loo], 
                                    ranks = ranks_lf, 
                                    threshold = 95,
                                    method = "hooi",
                                    calculate_bases = TRUE)
  
  end_tucker_lf <- Sys.time()
  
  timings[["tucker_lf"]] <- as.numeric(difftime(start_tucker_lf, end_tucker_lf, units = "secs"))
  gc()
  
  
  # --------------------------------------
  # Prepare structures for LF MCMC
  # --------------------------------------
  cat(paste0("\tLF MCMC Prep\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc_prep_lf <- Sys.time()
  
  mcmc_prep_lf_mf <- prepare_sf_mcmc(output = output_lf_4d_trans_scaled[ , , , -loo],
                                     tucker_outputs = tuck_lf_mf,
                                     ranks = ranks_lf)
  
  print(mcmc_prep_lf_mf$a_eta_prime / mcmc_prep_lf_mf$b_eta_prime)
  
  end_mcmc_prep_lf <- Sys.time()
  
  timings[["mcmc_prep_lf"]] <- as.numeric(difftime(start_mcmc_prep_lf, end_mcmc_prep_lf, units = "secs"))
  
  gc()
  
  # --------------------------------------
  # Discrep Tucker decomposition
  # --------------------------------------
  cat(paste0("\tDiscrep Tucker\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_tucker_discrep <- Sys.time()
  
  interp_reduced <- fast_knn_interpolate_first_dim(tuck_lf_mf$output_array_reduced, coord_lf, coord_hf, k = 1)
  discrep_interp <- output_hf_4d_unscaled[, , , -loo_hf_idx] - interp_reduced[, , , (runs_lf %in% runs_train_hf)[-loo]]
  
  ranks_discrep <- select_rank_variance(
    tnsr = as.tensor(discrep_interp),
    target_var = 0.80
  )
  
  r_s_discrep <- ranks_discrep[1]
  r_m_discrep <- ranks_discrep[2]
  r_y_discrep <- ranks_discrep[3]
  r_d_discrep <- ranks_discrep[4]
  
  cat("\t\tRanks Discrep:\t")
  cat(ranks_discrep)
  cat("\t(Target 80% explained var)\n\n")
  
  # note that the basis_matrices takes up A LOT of memory, avoid saving in RAM
  tuck_discrep <- apply_tucker_decomp(#output = discrep_4d_trans_scaled[ , , , -loo_hf_idx], 
    output = discrep_interp,
    ranks = ranks_discrep, 
    threshold = 70,
    method = "hosvd",
    calculate_bases = F)
  
  end_tucker_discrep <- Sys.time()
  
  timings[["tucker_discrep"]] <- as.numeric(difftime(start_tucker_discrep, end_tucker_discrep, units = "secs"))
  
  gc()
  
  
  # --------------------------------------
  # Interpolate LF Bases
  # --------------------------------------
  cat(paste0("\tBasis Interp\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_interp <- Sys.time()
  
  bases_interp <- interpolate_bases(tuck_lf_mf, k = 1)
  
  end_interp <- Sys.time()
  
  timings[["interp"]] <- as.numeric(difftime(start_interp, end_interp, units = "secs"))
  
  
  # --------------------------------------
  # Prepare MF MCMC
  # --------------------------------------
  cat(paste0("\tMF MCMC Prep\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc_prep_mf <- Sys.time()
  
  mcmc_prep_mf <- prepare_mf_mcmc(output_hf = output_hf_4d_trans_scaled[ , , , -loo_hf_idx],
                                  tucker_outputs_lf = tuck_lf_mf,
                                  tucker_outputs_discrep = tuck_discrep,
                                  interp_bases_outputs = bases_interp,
                                  ranks_lf = ranks_lf,
                                  ranks_discrep = ranks_discrep)
  
  print(mcmc_prep_mf$a_delta_prime / mcmc_prep_mf$b_delta_prime)
  
  end_mcmc_prep_mf <- Sys.time()
  
  timings[["mcmc_prep_mf"]] <- as.numeric(difftime(start_mcmc_prep_mf, end_mcmc_prep_mf, units = "secs"))
  
  
  # --------------------------------------
  # Run MCMC sampling
  # --------------------------------------
  cat(paste0("\tMF MCMC\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_mcmc <- Sys.time()
  
  n_chains <- 3 # as there are many chains, good for diagnosing convergence of multivariate chains (e.g. rho_mat)
  chains_mf <- vector("list", n_chains)
  
  for (i in 1:n_chains) {
    set.seed(1000 + i)  # different seed for each chain
    
    init_lambda_eta_i <- mcmc_prep_lf_mf$a_eta_prime / mcmc_prep_lf_mf$b_eta_prime
    init_lambda_w_vec_i <- runif(r_d_lf, 0.5, 5)
    init_rho_w_mat_i <- matrix(runif(r_d_lf * ncol(design_train_lf), 0.5, 5), nrow = r_d_lf, ncol = ncol(design_train_lf))
    
    init_lambda_delta_i <- mcmc_prep_mf$a_delta_prime / mcmc_prep_mf$b_delta_prime
    init_lambda_v_vec_i <- runif(r_d_discrep, 0.5, 5)
    init_rho_v_mat_i <- matrix(runif(r_d_discrep * ncol(design_train_hf), 0.5, 5), nrow = r_d_discrep, ncol = ncol(design_train_hf))
    
    chains_mf[[i]] <- perform_mf_mcmc(
      mcmc_lf_setup_outputs = mcmc_prep_lf_mf, 
      mcmc_mf_setup_outputs = mcmc_prep_mf,
      n_iter = 4000,
      burn_in = 2000,
      init_lambda_eta = init_lambda_eta_i,
      init_lambda_w_vec = init_lambda_w_vec_i,
      init_rho_w_mat = init_rho_w_mat_i,
      init_lambda_delta = init_lambda_delta_i,
      init_lambda_v_vec = init_lambda_v_vec_i,
      init_rho_v_mat = init_rho_v_mat_i,
      learning_rate = 150, # increase learning rate
      verbose = F
    )
    
    cat(paste0("\t  Chain ", i, " Complete ", format(Sys.time(), "%b %d %X"), "\n"))
  }
  
  end_mcmc <- Sys.time()
  
  timings[["mcmc"]] <- as.numeric(difftime(start_mcmc, end_mcmc, units = "secs"))
  
  
  # --------------------------------------
  # Check MCMC convergence
  # --------------------------------------  
  cat(paste0("\tCheck MCMC\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  # Geweke (within chain)
  retain_chains_mf <- check_geweke_multiple_chains(chains_mf, frac1 = 0.25, frac2 = 0.5)
  
  # # Gelman-Reuben (between chain)
  convergence_results <- check_mcmc_convergence_mf(retain_chains_mf)
  check_rhat(convergence_results$lambda_eta, name = "lambda_eta")
  check_rhat(convergence_results$lambda_w_vec, name = "lambda_w_vec")
  check_rhat(convergence_results$rho_w_mat, name = "rho_w_mat")
  
  check_rhat(convergence_results$lambda_delta, name = "lambda_delta")
  check_rhat(convergence_results$lambda_v_vec, name = "lambda_v_vec")
  check_rhat(convergence_results$rho_v_mat, name = "rho_v_mat")
  
  check_rhat(convergence_results$log_lik, name = "log_lik")
  
  # Global PSRF
  check_mpsrf(convergence_results$lambda_w_vec, name = "lambda_w_vec", threshold = 1.2)
  check_mpsrf(convergence_results$rho_w_mat, name = "rho_mat", threshold = 1.2)
  check_mpsrf(convergence_results$lambda_v_vec, name = "lambda_w_vec", threshold = 1.2)
  check_mpsrf(convergence_results$rho_v_mat, name = "rho_mat", threshold = 1.2)
  
  chains_mf <- retain_chains_mf
  
  
  # --------------------------------------
  # Evaluate posterior predictive distribution
  # --------------------------------------
  cat(paste0("\tMF Post Preds\t", format(Sys.time(), "%b %d %X"), "\n"))
  
  start_post_preds <- Sys.time()
  
  post_preds_mf_noise <- evaluate_mf_mcmc(design_unif[loo, ],
                                          tuck_lf_mf, tuck_discrep,
                                          bases_interp,
                                          mcmc_prep_lf_mf, mcmc_prep_mf,
                                          chains_mf,
                                          n_post_draws = 500,
                                          output_dim = dim(output_hf_4d_trans_scaled)[1:3],
                                          transform = "logistic-soft-clip",
                                          include_noise = TRUE,
                                          obs = output_hf_4d_unscaled[, , , loo_hf_idx])
  
  end_post_preds <- Sys.time()
  
  timings[["post_preds"]] <- as.numeric(difftime(start_post_preds, end_post_preds, units = "secs"))
  
  
  # --------------------------------------
  # Save outputs to file
  # --------------------------------------
  
  # Get rid of tuck basis arrays to save storage
  tuck_lf_mf$bases <- NULL
  gc()
  
  cat(paste0("\tSaving\t\t", format(Sys.time(), "%b %d %X"), "\n"))
  file_save <- paste0("save/mpas-mf-tensor/", loo, ".Rdata")
  save(ranks_lf, ranks_discrep, tuck_lf_mf, mcmc_prep_lf_mf, tuck_discrep, 
       bases_interp,
       mcmc_prep_mf, chains_mf, 
       post_preds_mf_noise, 
       timings,
       file = file_save)
  
  rm(ranks_lf, ranks_discrep, tuck_lf_mf, mcmc_prep_lf, tuck_discrep, bases_interp, mcmc_prep_mf, chains_mf, post_preds_mf, cv_stats_mf, timings)
  gc()
}

rm(list = ls())
gc()
