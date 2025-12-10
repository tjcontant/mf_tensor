###
### SETUP
###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(DiceKriging)

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


# -------------------------------
# LOO GP prediction function
# -------------------------------

max_hf_trans <- max(output_hf_4d_trans)

loo_cv_naive <- function(X, y, obs, covtype = "gauss", nugget = 1e-8, type = "SK") {
  
  y <- round(y, 6)
  
  # Check if y is constant
  if (sd(y) < 1e-12) {
    preds <- rep(mean(y), length(y))
    sds   <- rep(1e-12, length(y))
    
    draws <- rmvnorm(50, mu = preds, Sigma = diag(sds))
    
    # truncate
    # draws <- apply(draws, c(1, 2), function(x) pmin(pmax(x, 0), 1))
    draws <- apply(draws, c(1, 2), function(x) logistic_softclip(unscaled_transformed(x, max_hf_trans), epsilon = 0.0015)) 
    
    # stats
    preds <- colMeans(draws)
    sds   <- colSds(draws)
    
    lower <- apply(draws, 2, function(x) quantile(x, probs = 0.025))
    upper <- apply(draws, 2, function(x) quantile(x, probs = 0.975))
    
  } else {
    # Fit GP once
    fit <- DiceKriging::km(
      design       = X,
      response     = data.frame(y = y),
      covtype      = covtype,
      nugget       = nugget,
      nugget.estim = FALSE,
      control      = list(trace = FALSE)
    )
    
    # Compute LOO predictions
    loo <- DiceKriging::leaveOneOut.km(fit, type = type)
    
    # Transform back to original scale
    preds <- logistic_softclip(unscaled_transformed(loo$mean, max_hf_trans), epsilon = 0.0015)
    
    draws <- rmvnorm(200, mu=loo$mean, Sigma=diag(loo$sd))
    # draws_trunc <- pmin(pmax(draws, 0), 1)  # truncate only for coverage calculation
    draws_unscaled <- apply(draws, c(1, 2), function (x) logistic_softclip(unscaled_transformed(x, max_hf_trans), epsilon = 0.0015))
    
    # stats
    preds <- colMeans(draws_unscaled)
    sds   <- colSds(draws_unscaled)
    
    lower <- apply(draws_unscaled, 2, quantile, probs=0.025)
    upper <- apply(draws_unscaled, 2, quantile, probs=0.975)
    cover <- mean(obs >= lower & obs <= upper)
  }
  
  # Compute metrics
  mse <- (obs - preds)^2
  cover <- lower <= obs & obs <= upper
  
  return(list(mean = preds, sd = sds, mse = mse, cover = cover))
}


###
### Subset of spatial locations
###

set.seed(952)
spatial_reduction <- 0.01
s_hf_reduced_idx <- sample(1:s_hf, size = round(s_hf * spatial_reduction), replace = F)
s_hf_reduced <- length(s_hf_reduced_idx)

spatial_keep <- 1:s_hf %in% s_hf_reduced_idx

output_hf_4d_unscaled_reduced <- output_hf_4d_unscaled[spatial_keep, , , ]
output_hf_4d_trans_scaled_reduced <- output_hf_4d_trans_scaled[spatial_keep, , , ]


par(mfrow = c(1, 2))
plot(coord_hf$x[spatial_keep], coord_hf$y[spatial_keep], pch = 16, cex = 0.2)
plot(coord_hf, pch = 16, cex = 0.2)


t_total * s_hf_reduced *t_hf / 60 / 60 / 24


###
### INDIVIDUAL GP
###

naive_gp_mean <- array(NA, dim(output_hf_4d_unscaled_reduced))
naive_gp_sd <- array(NA, dim(output_hf_4d_unscaled_reduced))
naive_gp_mse <- array(NA, dim(output_hf_4d_unscaled_reduced))
naive_gp_cover <- array(NA, dim(output_hf_4d_unscaled_reduced))

pb <- progress_bar$new(
  format = " naive gp [:bar] :current/:total (:percent) Rate: :rate ETA: :eta",
  total = s_hf_reduced * t_hf, clear = FALSE
)

for (s in 1:s_hf_reduced) {
  for (m in 9) {
    for (year in 1) {
      s_idx <- s_hf_reduced_idx[s]
      gpr <- loo_cv_naive(
        X = design_unif[runs_hf, ], 
        y = output_hf_4d_trans_scaled[s_idx, m, year, ],
        obs = output_hf_4d_unscaled[s_idx, m, year, ]
      )
      
      naive_gp_mean[s, m, year, ] <- gpr$mean
      naive_gp_sd[s, m, year, ] <- gpr$sd
      naive_gp_mse[s, m, year, ] <- gpr$mse
      naive_gp_cover[s, m, year, ] <- gpr$cover
      
      pb$tick()
    }
  }
}

save(
  spatial_reduction,
  s_hf_reduced_idx,
  spatial_keep,
  output_hf_4d_unscaled_reduced,
  output_hf_4d_trans_scaled_reduced,
  naive_gp_mean,
  naive_gp_sd,
  naive_gp_mse,
  naive_gp_cover,
  file = "sav/mpas-naive-gp/mpas-naive_gp.Rdata"
)

# rm(list = ls())