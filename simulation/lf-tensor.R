setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved setup data
load("save/setup.Rdata")

# Define prior hyperparameters for high-fidelity emulator
sigma_lambda_w <- 1
a_rho_w <- 0
b_rho_w <- 1
a_eta <- 1
b_eta <- 0.5

# Clear console
cat("\014")

set.seed(543)

d_train_lf <- d_lf
runs_train_lf <- runs_lf
design_train_lf <- design[runs_train_lf, ]

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

tuck_lf <- apply_tucker_decomp(output = output_lf, 
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

mcmc_prep_lf <- prepare_sf_mcmc(output = output_lf,
                                tucker_outputs = tuck_lf,
                                ranks = ranks_lf)

end_mcmc_prep <- Sys.time()

timings[["mcmc_prep"]] <- as.numeric(difftime(start_mcmc_prep, end_mcmc_prep, units = "secs"))


# --------------------------------------
# Run MCMC sampling
# --------------------------------------
cat(paste0("\tMCMC\t\t", format(Sys.time(), "%b %d %X"), "\n"))

start_mcmc <- Sys.time()

n_chains <- 2 # as there are many chains, good for diagnosing convergence of multivariate chains (e.g. rho_mat)
chains_lf <- vector("list", n_chains)
for (i in 1:n_chains) {
  set.seed(1000 + i)  # different seed for each chain
  
  # Generate different initial values for each chain
  init_lambda_eta_i <- mcmc_prep_lf$a_eta_prime / mcmc_prep_lf$b_eta_prime
  init_lambda_w_vec_i <- runif(r_d_lf, 0.5, 5)
  init_rho_mat_i <- matrix(runif(r_d_lf * ncol(design_train_lf), 0.5, 5), nrow = r_d_lf, ncol = ncol(design_train_lf))
  
  chains_lf[[i]] <- perform_sf_mcmc(
    output = output_lf,
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

test_input <- matrix(c(-1.7, -1.7, -1.7), nrow = 1)
true_output <- simulate_for_design(test_input, 
                                   coord_lf, coord_hf, 
                                   year_ids, month_ids,
                                   start_lf = 1, end_lf = 5)


post_preds_lf_noise <- evaluate_sf_mcmc(
  untested_input = test_input,
  mcmc_setup_outputs = mcmc_prep_lf,
  mcmc_outputs = chains_lf,
  tucker_outputs = tuck_lf,
  n_post_draws = 200,
  design_train = design_train_lf,
  output_dim = dim(output_lf)[1:3],
  ranks = ranks_lf,
  hf = F,
  transform = "none",
  include_noise = TRUE,
  obs = true_output$high_fi[, , , , drop = T]
)

end_post_preds <- Sys.time()

timings[["post_preds"]] <- as.numeric(difftime(start_post_preds, end_post_preds, units = "secs"))


# --------------------------------------
# Save outputs to file
# --------------------------------------
cat(paste0("\tSaving\t\t", format(Sys.time(), "%b %d %X"), "\n"))
file_save <- paste0("save/lf-tensor.Rdata")
save(ranks_lf, tuck_lf, mcmc_prep_lf, chains_lf, 
     post_preds_lf_noise, 
     timings,
     file = file_save)

# Clean up environment
rm(list = ls())
gc()
