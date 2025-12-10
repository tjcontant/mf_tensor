setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved setup data
load("save/setup.Rdata")


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

set.seed(543)

# Clear console
cat("\014")

d_train_lf <- d_lf
runs_train_lf <- runs_lf
design_train_lf <- design[runs_train_lf, ]

d_train_hf <- d_hf
runs_train_hf <- runs_hf
design_train_hf <- design[runs_train_hf, ]

timings <- list()


# --------------------------------------
# Set Tucker ranks
# --------------------------------------  
cat(paste0("\tRank Select.\t", format(Sys.time(), "%b %d %X"), "\n"))

start_rank_select <- Sys.time()

ranks_lf <- select_rank_variance(
  tnsr = as.tensor(output_lf),
  target_var = 0.99
)

ranks_hf <- select_rank_variance(
  tnsr = as.tensor(output_hf),
  target_var = 0.99
)

ranks_lf <- pmax(ranks_lf, ranks_hf)

r_s_lf <- ranks_lf[1]
r_m_lf <- ranks_lf[2]
r_y_lf <- ranks_lf[3]
r_d_lf <- ranks_lf[4]

cat("\t\tRanks LF:\t")
cat(ranks_lf)
cat("\t(Target 99% explained var)\n")


end_rank_select <- Sys.time()

timings[["rank_select"]] <- as.numeric(difftime(start_rank_select, end_rank_select, units = "secs"))


# --------------------------------------
# LF Tucker decomposition
# --------------------------------------
cat(paste0("\tLF Tucker\t", format(Sys.time(), "%b %d %X"), "\n"))

start_tucker_lf <- Sys.time()

tuck_lf_mf <- apply_tucker_decomp(output = output_lf, 
                                  ranks = ranks_lf, 
                                  threshold = 95,
                                  method = "hooi",
                                  calculate_bases = TRUE)

end_tucker_lf <- Sys.time()

timings[["tucker_lf"]] <- as.numeric(difftime(start_tucker_lf, end_tucker_lf, units = "secs"))


# --------------------------------------
# Prepare structures for LF MCMC
# --------------------------------------
cat(paste0("\tLF MCMC Prep\t", format(Sys.time(), "%b %d %X"), "\n"))

start_mcmc_prep_lf <- Sys.time()

mcmc_prep_lf_mf <- prepare_sf_mcmc(output = output_lf,
                                   tucker_outputs = tuck_lf_mf,
                                   ranks = ranks_lf)

end_mcmc_prep_lf <- Sys.time()

timings[["mcmc_prep_lf"]] <- as.numeric(difftime(start_mcmc_prep_lf, end_mcmc_prep_lf, units = "secs"))


# --------------------------------------
# Interpolate LF Bases
# --------------------------------------
cat(paste0("\tBasis Interp\t", format(Sys.time(), "%b %d %X"), "\n"))

start_interp <- Sys.time()

bases_interp <- interpolate_bases(tuck_lf_mf, k = 10)

end_interp <- Sys.time()

timings[["interp"]] <- as.numeric(difftime(start_interp, end_interp, units = "secs"))


# --------------------------------------
# Discrep Tucker decomposition
# --------------------------------------
cat(paste0("\tDiscrep Tucker\t", format(Sys.time(), "%b %d %X"), "\n"))

start_tucker_discrep <- Sys.time()

low_interp <- ttl(as.tensor(tuck_lf_mf$core_tensor), bases_interp$factor_matrices_interp, 1:4)
low_interp <- low_interp@data
discrep <- output_hf - low_interp[, , , runs_hf]

# discrep <- output_hf - knn_interpolate_first_dim(tuck_lf_mf$output_array_reduced, coord_lf, coord_hf, k = 10)[, , , runs_hf]

ranks_discrep <- select_rank_variance(
  tnsr = as.tensor(discrep),
  target_var = 0.99
)

r_s_discrep <- ranks_discrep[1]
r_m_discrep <- ranks_discrep[2]
r_y_discrep <- ranks_discrep[3]
r_d_discrep <- ranks_discrep[4]

cat("\t\tRanks Discrep:\t")
cat(ranks_discrep)
cat("\t(Target 99% explained var)\n\n")

# note that the basis_matrices takes up A LOT of memory, avoid saving in RAM
tuck_discrep <- apply_tucker_decomp(output = discrep, 
                                    ranks = ranks_discrep, 
                                    threshold = 70,
                                    method = "hooi",
                                    calculate_bases = F)

end_tucker_discrep <- Sys.time()

timings[["tucker_discrep"]] <- as.numeric(difftime(start_tucker_discrep, end_tucker_discrep, units = "secs"))


# --------------------------------------
# Prepare MF MCMC
# --------------------------------------
cat(paste0("\tMF MCMC Prep\t", format(Sys.time(), "%b %d %X"), "\n"))

start_mcmc_prep_mf <- Sys.time()

mcmc_prep_mf <- prepare_mf_mcmc(output_hf = output_hf,
                                tucker_outputs_lf = tuck_lf_mf,
                                tucker_outputs_discrep = tuck_discrep,
                                interp_bases_outputs = bases_interp,
                                ranks_lf = ranks_lf,
                                ranks_discrep = ranks_discrep)

end_mcmc_prep_mf <- Sys.time()

timings[["mcmc_prep_mf"]] <- as.numeric(difftime(start_mcmc_prep_mf, end_mcmc_prep_mf, units = "secs"))


# --------------------------------------
# Run MCMC sampling
# --------------------------------------
cat(paste0("\tMF MCMC\t\t", format(Sys.time(), "%b %d %X"), "\n"))

start_mcmc <- Sys.time()

n_chains <- 2
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
    learning_rate = 125,
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

# Gelman-Reuben (between chain)
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

test_input <- matrix(c(-1.7, -1.7, -1.7), nrow = 1)
true_output <- simulate_for_design(test_input, 
                                   coord_lf, coord_hf, 
                                   year_ids, month_ids,
                                   start_lf = 1, end_lf = 5)

post_preds_mf_noise <- evaluate_mf_mcmc(test_input,
                                        tuck_lf_mf, tuck_discrep,
                                        bases_interp,
                                        mcmc_prep_lf_mf, mcmc_prep_mf,
                                        chains_mf,
                                        n_post_draws = 200,
                                        output_dim = dim(output_hf)[1:3],
                                        transform = "none",
                                        include_noise = TRUE,
                                        obs = true_output$high_fi[, , , , drop = T])

end_post_preds <- Sys.time()

timings[["post_preds"]] <- as.numeric(difftime(start_post_preds, end_post_preds, units = "secs"))


# --------------------------------------
# Save outputs to file
# --------------------------------------
cat(paste0("\tSaving\t\t", format(Sys.time(), "%b %d %X"), "\n"))
file_save <- paste0("save/mf-tensor.Rdata")
save(ranks_lf, ranks_discrep, tuck_lf_mf, mcmc_prep_lf_mf, tuck_discrep,
     bases_interp, # remove to save computer storage
     mcmc_prep_mf, chains_mf,
     post_preds_mf_noise, 
     timings,
     file = file_save)

rm(ranks_lf, ranks_discrep, tuck_lf_mf, mcmc_prep_lf_mf, tuck_discrep, bases_interp, mcmc_prep_mf, chains_mf, post_preds_mf, cv_stats_mf, timings)
gc()
