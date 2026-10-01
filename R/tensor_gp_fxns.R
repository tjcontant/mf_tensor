# Bayesian tensor-GP model: Tucker-reduced design-mode coefficients are each an
# independent GP over the design space, fit via adaptive Metropolis-Hastings
# MCMC

library(progress)
library(rTensor)
library(Matrix)
library(mvtnorm)
library(matrixStats)
library(truncnorm)
library(pbapply)

# --- Kernel / distance-cache / density helpers ---

R_mat <- function(design_mat_1, design_mat_2, lambda_w, rho_vec) {
  # Scale each column of both design matrices by the inverse of rho_vec
  scaled_1 <- sweep(design_mat_1, 2, rho_vec, FUN = "/")
  scaled_2 <- sweep(design_mat_2, 2, rho_vec, FUN = "/")

  # Compute pairwise Euclidean distances between all rows
  dists <- fields::rdist(scaled_1, scaled_2)

  # Apply Gaussian RBF kernel: K_ij = exp(-0.5 * squared distance)
  K <- exp(-0.5 * dists^2)

  # Scale kernel matrix by 1 / lambda_w
  out <- (1 / lambda_w) * K

  return(out)
}

Sigma_general_fast <- function(dists_sq_list, lambda_vec) {
  r3 <- length(lambda_vec) # Number of groups
  d1 <- nrow(dists_sq_list[[1]]) # Number of points in design_mat_1 per group
  d2 <- ncol(dists_sq_list[[1]]) # Number of points in design_mat_2 per group

  # Pre-allocate output covariance matrix
  out <- matrix(0, d1 * r3, d2 * r3)

  # Fill each block with group-specific Gaussian kernels
  for (i in 1:r3) {
    row_idx_1 <- ((i - 1) * d1 + 1):(i * d1)
    row_idx_2 <- ((i - 1) * d2 + 1):(i * d2)

    # Compute Gaussian kernel from squared distances
    K <- exp(-0.5 * dists_sq_list[[i]])

    # Scale kernel by inverse lambda and assign to block
    out[row_idx_1, row_idx_2] <- (1 / lambda_vec[i]) * K
  }

  return(out)
}

# Precomputes RAW (rho-independent) squared per-dimension distances between two
# design matrices
initialize_raw_dist_cache <- function(design_mat_1, design_mat_2) {
  p <- ncol(design_mat_1)
  lapply(1:p, function(d) outer(design_mat_1[, d], design_mat_2[, d], "-")^2)
}

# Builds the r3-component covariance matrix from a raw per-dimension squared
# distance cache (initialize_raw_dist_cache()) for a given lambda_vec/rho_mat
Sigma_general_from_raw <- function(raw_dist_cache, lambda_vec, rho_mat, r3 = NULL) {
  if (is.null(r3)) r3 <- nrow(rho_mat)
  d1 <- nrow(raw_dist_cache[[1]])
  d2 <- ncol(raw_dist_cache[[1]])
  out <- matrix(0, d1 * r3, d2 * r3)
  for (i in 1:r3) {
    rho_i <- rho_mat[i, ]
    scaled_sq_dist <- Reduce(`+`, Map(function(rd, rho) rd / rho^2, raw_dist_cache, rho_i))
    K <- exp(-0.5 * scaled_sq_dist)
    row_idx <- ((i - 1) * d1 + 1):(i * d1)
    col_idx <- ((i - 1) * d2 + 1):(i * d2)
    out[row_idx, col_idx] <- (1 / lambda_vec[i]) * K
  }
  out
}

initialize_distance_cache_general <- function(design_1, design_2, rho_mat) {
  r3 <- nrow(rho_mat) # Number of components/groups
  dists_sq_list <- vector("list", r3) # Preallocate list for squared distance matrices

  for (i in 1:r3) {
    rho_i <- rho_mat[i, ] # Range parameters for component i

    # Scale columns of design matrices by 1 / rho_i (element-wise division)
    scaled_1 <- sweep(design_1, 2, rho_i, "/") # (n1 × d)
    scaled_2 <- sweep(design_2, 2, rho_i, "/") # (n2 × d)

    # Compute pairwise squared Euclidean distances between scaled points
    dists_sq_list[[i]] <- fields::rdist(scaled_1, scaled_2)^2
  }

  return(dists_sq_list)
}

# Updates the squared distance matrix cache for a single component
update_distance_cache_element_general <- function(design_1, design_2, rho_mat, dists_sq_list, i) {
  # Extract range parameters for component i
  rho_i <- rho_mat[i, ]

  # Scale design matrices by 1 / rho_i
  scaled_1 <- sweep(design_1, 2, rho_i, "/")
  scaled_2 <- sweep(design_2, 2, rho_i, "/")

  # Recompute squared distances and update cache element i
  dists_sq_list[[i]] <- fields::rdist(scaled_1, scaled_2)^2

  return(dists_sq_list)
}

# Computes the density of the half-normal distribution at given points
dhalfnorm <- function(x, sigma = 1, log = FALSE) {
  # Density is zero for negative values (half-normal is zero for x < 0)
  if (any(x < 0)) {
    return(0)
  }

  # Compute log density using normal density and adding log(2)
  log_density <- log(2) + dnorm(x, mean = 0, sd = sigma, log = TRUE)

  # Return log density or density according to user request
  if (log) {
    return(log_density)
  }
  return(exp(log_density))
}

# Computes the log-density of a multivariate normal distribution using Cholesky
# decomposition
log_mvnorm_cholesky <- function(x, mu, Sigma = NULL, L = NULL) {
  # Cholesky decomposition: Sigma = L * L^T. Pass a precomputed L (e.g. one
  # already produced by a PD check) to skip recomputing it here.
  if (is.null(L)) L <- base::chol(Sigma)

  # Solve for z = L^{-T} * (x - mu)
  # Equivalent to solving L * z = (x - mu), with L upper-triangular
  z <- backsolve(L, x - mu, transpose = TRUE)

  # Compute log determinant: log(det(Sigma)) = 2 * sum(log(diag(L)))
  log_det_Sigma <- 2 * sum(log(diag(L)))

  # Compute quadratic form: (x - mu)^T * Sigma^{-1} * (x - mu) = z^T * z
  quadratic_form <- sum(z^2)

  # Final log-density formula:
  # -0.5 * (log det + quadratic term + d * log(2π))
  log_density <- -0.5 * (log_det_Sigma + quadratic_form + length(x) * log(2 * pi))

  return(log_density)
}

matrix_calc <- function(vectors, A, x = NULL, d) {
  r <- length(vectors) # Number of vectors
  n <- length(vectors[[1]]) # Length of each vector (s x t)
  r3 <- ncol(A) # Number of output columns per block

  # Progress bar for weighted vector combinations
  pb1 <- progress_bar$new(
    format = "  Combining vectors [:bar] :percent eta: :eta",
    total = r3, clear = TRUE, width = 70
  )
  results <- list()
  for (j in 1:r3) {
    # Weighted combination of all k vectors using column j of A
    weights <- A[, j]
    vec <- Reduce(`+`, Map(`*`, weights, vectors))
    results[[j]] <- vec
    pb1$tick()
  }

  # Progress bar for coef_matrix calculation (r3 x r3 iterations)
  pb2 <- progress_bar$new(
    format = "  Computing coef_matrix [:bar] :percent eta: :eta",
    total = r3 * r3, clear = TRUE, width = 70
  )

  coef_matrix <- matrix(0, nrow = r3, ncol = r3)
  # Compute dot product for each pair
  for (i in 1:r3) {
    for (j in 1:r3) {
      coef_matrix[i, j] <- sum(results[[i]] * results[[j]])
      pb2$tick()
    }
  }

  solve_Ct_C <- kronecker(solve(coef_matrix), Diagonal(d))

  if (!is.null(x)) {
    total_iters <- r3 * d
    pb3 <- progress_bar$new(
      format = "  Computing Ct_x [:bar] :percent eta: :eta",
      total = total_iters, clear = TRUE, width = 70
    )
    results_2 <- vector("numeric", total_iters)
    counter <- 1
    for (i in 1:r3) {
      for (j in 1:d) {
        idx <- ((j - 1) * n + 1):(j * n)
        dot_product <- sum(results[[i]] * x[idx]) # scalar
        results_2[counter] <- dot_product
        counter <- counter + 1
        pb3$tick()
      }
    }

    Ct_x <- unlist(results_2)

    return(list(
      solve_Ct_C = solve_Ct_C,
      Ct_x       = Ct_x
    ))
  }

  return(list(
    solve_Ct_C = solve_Ct_C,
  ))
}

flatten_last_dim <- function(arr) {
  d <- dim(arr)
  dim(arr) <- c(prod(d[-length(d)]), d[length(d)])
  return(arr)
}

drop_first_dim <- function(x) {
  dim_x <- dim(x)
  if (!is.null(dim_x) && dim_x[1] == 1) {
    array(x[1, , , drop = FALSE], dim = dim_x[-1])
  } else {
    x
  }
}

# --- LF basis interpolation (used by the MF driver to bring the LF
# Tucker reconstruction onto the HF grid) ---

interpolate_bases <- function(tucker_output_lf, coord_lf, coord_hf, k) {
  G <- tucker_output_lf$core_tensor
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)

  factor_matrices_lf <- tucker_output_lf$factor_matrices
  factor_matrices_lf_interp <- factor_matrices_lf
  factor_matrices_lf_interp[[1]] <- fast_knn_interpolate_first_dim(factor_matrices_lf[[1]], coord_lf, coord_hf, k)

  low_interp <- ttl(as.tensor(G), factor_matrices_lf_interp, 1:4)

  B_eta_tile_x_G_unfold <- Reduce(
    function(x, y) kronecker(y, x),
    rev(factor_matrices_lf_interp[3:1])
  ) %*% G_unfold

  return(list(
    B_eta_tile_x_G_unfold = B_eta_tile_x_G_unfold,
    low_interp = low_interp@data,
    factor_matrices_interp = factor_matrices_lf_interp
  ))
}

# --- Single-fidelity model (used identically for LF-only and HF-only) ---

prepare_sf_mcmc <- function(output, tucker_outputs, ranks, a_eta, b_eta) {
  # Number of simulations (mode-3 dimension)
  d_train <- dim(output)[length(dim(output))]

  # Extract core tensor and vectorized basis matrix
  G <- tucker_outputs$core_tensor # Core tensor (3D array)
  bases <- tucker_outputs$bases

  # Unfold the core tensor along last mode (simulations) and transpose
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)

  # Flatten the 3D output array into a long vector
  output_vec <- as.vector(output)

  r_d <- ranks[length(ranks)]

  n_out <- prod(dim(output)[-length(dim(output))])

  # if bases already calculated
  if (!is.null(tucker_outputs$bases)) {
    matrix_calc_results <- matrix_calc(
      pblapply(bases, as.vector),
      G_unfold,
      output_vec,
      d_train
    )

    solve_Ct_C <- matrix_calc_results$solve_Ct_C
    Ct_x <- matrix_calc_results$Ct_x

    # bases not already calculated
  } else {
    index_grid <- expand.grid(lapply(ranks[-length(ranks)], function(n) 1:n))

    pb <- progress_bar$new(
      total = r_d * nrow(index_grid) + r_d * d_train,
      format = " Calculating Matrices [:bar] :percent eta: :eta",
      clear = TRUE,
      width = 70
    )

    # calculate c_k
    c_vecs <- list()

    for (k in 1:r_d) {
      out <- rep(0, n_out)
      for (i in 1:nrow(index_grid)) {
        b <- as.vector(calculate_basis(
          U_list = tucker_outputs$factor_matrices,
          idx = index_grid[i, ]
        ))

        out <- out + G_unfold[i, k] * b

        pb$tick()
      }
      c_vecs[[k]] <- out
    }

    n <- length(c_vecs)
    C_coef <- matrix(0, n, n)

    for (i in 1:n) {
      for (j in 1:n) {
        C_coef[i, j] <- sum(c_vecs[[j]] * c_vecs[[i]])
      }
    }

    solve_Ct_C <- kronecker(solve(C_coef), diag(d_train))

    Ct_x <- c()
    for (k in 1:r_d) {
      for (d in 1:d_train) {
        nd <- length(dim(output))
        z_d <- as.vector(do.call(`[`, c(list(output), rep(list(TRUE), nd - 1), list(d))))
        Ct_x <- c(
          Ct_x,
          t(c_vecs[[k]]) %*% z_d
        )

        pb$tick()
      }
    }
  }

  # Estimate gamma_hat (posterior mode/mean for regression coefficients)
  gamma_hat <- solve_Ct_C %*% Ct_x

  a_eta_prime <- a_eta + (d_train * (n_out - r_d)) / 2

  b_eta_prime <- b_eta + 0.5 * (
    sum(output_vec^2) -
      t(Ct_x) %*% solve_Ct_C %*% Ct_x
  )
  b_eta_prime <- as.numeric(b_eta_prime) # Ensure scalar output

  return(list(
    solve_C_t_C  = solve_Ct_C,
    gamma_hat    = gamma_hat,
    a_eta_prime  = a_eta_prime,
    b_eta_prime  = b_eta_prime
  ))
}

# Runs adaptive Metropolis-Hastings MCMC for single-fidelity model parameters
perform_sf_mcmc <- function(
  output, design, mcmc_setup_outputs, n_iter, burn_in,
  r3,
  init_lambda_eta = NULL,
  init_lambda_w_vec = NULL,
  init_rho_mat = NULL,
  omega_init = list(
    lambda_eta = 0.01,
    lambda_w = 0.5,
    rho = 0.3
  ),
  learning_rate,
  prior,
  progress_file = NULL
) {
  d_train <- dim(output)[length(dim(output))]

  # Extract fixed components
  solve_C_t_C <- mcmc_setup_outputs$solve_C_t_C
  gamma_hat <- mcmc_setup_outputs$gamma_hat
  a_eta_prime <- mcmc_setup_outputs$a_eta_prime
  b_eta_prime <- mcmc_setup_outputs$b_eta_prime

  # Initialize parameters
  lambda_eta <- init_lambda_eta %||% (a_eta_prime / b_eta_prime)
  lambda_w_vec <- init_lambda_w_vec %||% rep(1, r3)
  rho_mat <- init_rho_mat %||% matrix(2, nrow = r3, ncol = ncol(design))

  # Initialize squared-distance cache
  dists_sq_list <- initialize_distance_cache_general(design, design, rho_mat)

  # Log-posterior function
  log_posterior <- function(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat) {
    Sigma <- lambda_eta^(-1) * solve_C_t_C + Sigma_general_fast(dists_sq_list, lambda_w_vec)

    # chol() doubles as the positive-definiteness check
    L <- tryCatch(base::chol(Sigma), error = function(e) NULL)
    if (is.null(L)) {

      Sigma <- as.matrix(forceSymmetric(
        nearPD(Sigma, ensureSymmetry = TRUE)$mat
      ))

      # Add small jitter to diagonal to ensure strict PD
      L <- base::chol(Sigma)
    }

    ll <- log_mvnorm_cholesky(
      as.vector(gamma_hat),
      mu = rep(0, d_train * r3),
      L = L
    )
    # Gamma prior times the likelihood of the Tucker residual (a_eta_prime, b_eta_prime)
    prior_lambda_eta <- dgamma(lambda_eta, a_eta_prime, b_eta_prime, log = TRUE)
    prior_lambda_w <- sum(dhalfnorm(lambda_w_vec, prior$sigma_lambda_w, log = TRUE))
    prior_rho <- sum(dlnorm(rho_mat, prior$a_rho_w, prior$b_rho_w, log = TRUE))

    return(as.numeric(ll + prior_lambda_eta + prior_lambda_w + prior_rho))
  }

  # --- Metropolis updates ---
  # log_post_current is the already-known log-posterior at the current state
  update_lambda_eta_MH <- function(lambda_eta, ..., log_post_current) {
    # Propose from truncated normal (truncated at 0 from below)
    lambda_eta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_eta, sd = omega_lambda_eta)

    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_eta_prop, a = 0, b = Inf, mean = lambda_eta, sd = omega_lambda_eta) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_eta, a = 0, b = Inf, mean = lambda_eta_prop, sd = omega_lambda_eta) |> log()

    log_post_prop <- log_posterior(lambda_eta_prop, ...)

    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = lambda_eta_prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_eta, accept = FALSE, log_post = log_post_current))
    }
  }

  update_lambda_w_vec_i_MH <- function(i, lambda_eta, lambda_w_vec, ..., log_post_current) {
    prop <- lambda_w_vec

    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega_lambda_w[i])
    prop[i] <- prop_i_new

    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega_lambda_w[i]) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_w_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega_lambda_w[i]) |> log()

    log_post_prop <- log_posterior(lambda_eta, prop, ...)

    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_w_vec, accept = FALSE, log_post = log_post_current))
    }
  }

  update_rho_mat_ij_MH <- function(i, j, lambda_eta, lambda_w_vec, rho_mat, dists_sq_list, log_post_current) {
    rho_pro <- rho_mat

    # Propose new rho[i,j] from truncated normal > 0
    rho_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_mat[i, j], sd = omega_rho[i, j])
    rho_pro[i, j] <- rho_prop_ij

    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_prop_ij, a = 0, b = Inf, mean = rho_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_mat[i, j], a = 0, b = Inf, mean = rho_prop_ij, sd = omega_rho[i, j]) |> log()

    # Update the distance cache with the proposed rho
    dists_sq_list_pro <- update_distance_cache_element_general(design, design, rho_pro, dists_sq_list, i)

    log_post_prop <- log_posterior(lambda_eta, lambda_w_vec, dists_sq_list_pro, rho_pro)

    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = rho_pro, dists_sq_list = dists_sq_list_pro, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = rho_mat, dists_sq_list = dists_sq_list, accept = FALSE, log_post = log_post_current))
    }
  }

  # Storage
  n_save <- n_iter - burn_in
  lambda_eta_chain <- numeric(n_save)
  lambda_w_vec_chain <- matrix(NA, nrow = n_save, ncol = r3)
  rho_mat_chain <- array(NA, dim = c(n_save, r3, ncol(design)))
  log_lik_chain <- numeric(n_save)

  # Initialize step sizes
  omega_lambda_eta <- omega_init$lambda_eta
  omega_lambda_w <- rep(omega_init$lambda_w, r3)
  omega_rho <- matrix(omega_init$rho, r3, ncol(design))

  accept_lambda_eta <- 0
  accept_lambda_w <- rep(0, r3)
  accept_rho <- matrix(0, r3, ncol(design))
  update_omega_iter <- 50

  target_accept <- 0.45

  pb <- progress_bar$new(
    format = "  SF MCMC [:bar] :current/:total (:percent) ETA: :eta",
    total = n_iter, width = 70
  )

  # Current log-posterior, carried forward across MH updates
  log_post_curr <- log_posterior(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat)

  for (iter in 1:n_iter) {
    # lambda_eta
    res <- update_lambda_eta_MH(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat, log_post_current = log_post_curr)
    lambda_eta <- res$value
    log_post_curr <- res$log_post
    if (iter <= burn_in) {
      accept_lambda_eta <- accept_lambda_eta + res$accept
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_eta / update_omega_iter
        old_omega <- omega_lambda_eta
        omega_lambda_eta <- omega_lambda_eta * exp((rate - target_accept) / sqrt(iter) * learning_rate)
        accept_lambda_eta <- 0
      }
    }

    # lambda_w
    for (j in 1:r3) {
      res <- update_lambda_w_vec_i_MH(j, lambda_eta, lambda_w_vec, dists_sq_list, rho_mat, log_post_current = log_post_curr)
      lambda_w_vec <- res$value
      log_post_curr <- res$log_post
      if (iter <= burn_in) {
        accept_lambda_w[j] <- accept_lambda_w[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_w[j] / update_omega_iter
          old_omega <- omega_lambda_w[j]
          omega_lambda_w[j] <- omega_lambda_w[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_w[j] <- 0
        }
      }
    }

    # rho
    for (i_rho in 1:r3) {
      for (j_rho in 1:ncol(design)) {
        res <- update_rho_mat_ij_MH(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_mat, dists_sq_list, log_post_current = log_post_curr)
        rho_mat <- res$value
        dists_sq_list <- res$dists_sq_list
        log_post_curr <- res$log_post
        if (iter <= burn_in) {
          accept_rho[i_rho, j_rho] <- accept_rho[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho[i_rho, j_rho] / update_omega_iter
            old_omega <- omega_rho[i_rho, j_rho]
            omega_rho[i_rho, j_rho] <- omega_rho[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho[i_rho, j_rho] <- 0
          }
        }
      }
    }

    # Store samples
    if (iter > burn_in) {
      idx <- iter - burn_in
      lambda_eta_chain[idx] <- lambda_eta
      lambda_w_vec_chain[idx, ] <- lambda_w_vec
      rho_mat_chain[idx, , ] <- rho_mat
      log_lik_chain[idx] <- log_post_curr
    }

    if (!is.null(progress_file)) cat(iter, file = progress_file)

    pb$tick()
  }

  return(list(
    n_iter = n_iter,
    burn_in = burn_in,
    init_params = list(
      lambda_eta_init = init_lambda_eta,
      lambda_w_vec_init = init_lambda_w_vec,
      rho_mat_init = init_rho_mat
    ),
    post_chains = list(
      lambda_eta_chain = lambda_eta_chain,
      lambda_w_vec_chain = lambda_w_vec_chain,
      rho_mat_chain = rho_mat_chain,
      log_lik_chain = log_lik_chain
    ),
    omega = list(
      omega_lambda_eta,
      omega_lambda_w,
      omega_rho
    )
  ))
}

# Predictions at the rows of untested_input; V_11's Cholesky is computed once
# per draw and shared across all test points
evaluate_sf_mcmc <- function(untested_input,
                             mcmc_setup_outputs,
                             mcmc_outputs,
                             tucker_outputs,
                             n_post_draws,
                             design_train,
                             output_dim,
                             ranks,
                             hf = FALSE,
                             coord_lf = NULL,
                             coord_hf = NULL,
                             include_noise = FALSE,
                             obs = NULL,
                             verbose = TRUE) {
  n_chains <- length(mcmc_outputs)

  r_d <- ranks[length(ranks)]

  n_iter <- mcmc_outputs[[1]]$n_iter
  burn_in <- mcmc_outputs[[1]]$burn_in

  bases <- tucker_outputs$bases

  if (!is.null(bases)) {
    B_eta <- do.call(cbind, lapply(bases, as.vector))
  } else {
    n <- prod(output_dim)
    n_r <- prod(tucker_outputs$ranks[-length(tucker_outputs$ranks)])
    B_eta <- matrix(NA, n, n_r)

    index_grid <- expand.grid(lapply(ranks[-length(ranks)], function(n) 1:n))

    pb <- progress_bar$new(
      format = "  Calculating Bases [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
      total = nrow(index_grid),
      width = 70
    )

    for (i in 1:nrow(index_grid)) {
      B_eta[, i] <- as.vector(calculate_basis(
        U_list <- tucker_outputs$factor_matrices,
        idx = index_grid[i, ]
      ))

      pb$tick()
    }
  }

  solve_C_t_C <- mcmc_setup_outputs$solve_C_t_C

  gamma_hat_1 <- mcmc_setup_outputs$gamma_hat

  G <- tucker_outputs$core_tensor
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)

  post_pred_sf <- function(untested_input, n_post_draws) {
    untested_input <- as.matrix(untested_input)
    n_new <- nrow(untested_input)

    # Raw (rho-independent) squared distances, computed ONCE and reused for
    # every posterior draw's covariance construction below (see
    # Sigma_general_from_raw()'s comment)
    raw_11 <- initialize_raw_dist_cache(design_train, design_train)
    raw_12 <- initialize_raw_dist_cache(design_train, untested_input)
    raw_22 <- initialize_raw_dist_cache(untested_input, untested_input)

    gamma_output_list <- vector("list", n_new)
    pred_output_list <- vector("list", n_new)
    for (k in 1:n_new) {
      gamma_output_list[[k]] <- list()
      pred_output_list[[k]] <- list()
    }

    count <- 0
    attempts <- 0
    max_attempts <- 3 * n_post_draws # To prevent infinite loops

    if (verbose) {
      pb <- progress_bar$new(
        format = ifelse(hf,
          "  HF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
          "  LF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate"
        ),
        total = n_post_draws,
        width = 70
      )
    }

    while (count < n_post_draws && attempts < max_attempts) {
      attempts <- attempts + 1

      # Sample chain and draw index
      chain_i <- sample(1:n_chains, 1)
      chain_length <- length(mcmc_outputs[[chain_i]]$post_chains$lambda_eta_chain)
      draw_i <- sample(1:chain_length, 1)

      chain <- mcmc_outputs[[chain_i]]$post_chains

      lambda_eta <- chain$lambda_eta_chain[draw_i]
      lambda_w_vec <- chain$lambda_w_vec_chain[draw_i, ]
      rho_w_mat <- chain$rho_mat_chain[draw_i, , ]

      # Ensure rho_w_mat is a proper matrix
      if (is.null(dim(rho_w_mat))) {
        rho_w_mat <- matrix(rho_w_mat,
          nrow = dim(chain$rho_mat_chain)[2],
          ncol = dim(chain$rho_mat_chain)[3]
        )
      }

      # Construct GP components
      A <- Sigma_general_from_raw(raw_11, lambda_w_vec, rho_w_mat, r3 = r_d)
      B <- (1 / lambda_eta * solve_C_t_C)

      V_11 <- A + B

      # Predictive mean and covariance
      eigvals <- tryCatch(
        eigen(forceSymmetric(V_11), symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA
      )
      if (any(is.na(eigvals)) || any(eigvals <= 0)) {
        next # skip non-PD draw
      }

      L <- chol(V_11)
      y <- backsolve(L, forwardsolve(t(L), gamma_hat_1))

      # All n_new test points' cross-covariance built and back-solved in ONE
      # shot (V_12_all has r_d columns per test point).
      V_12_all <- Sigma_general_from_raw(raw_12, lambda_w_vec, rho_w_mat, r3 = r_d)
      V_22_all <- Sigma_general_from_raw(raw_22, lambda_w_vec, rho_w_mat, r3 = r_d)
      Y_all <- backsolve(L, forwardsolve(t(L), V_12_all))

      # A failed draw at any test point skips this posterior draw
      draws_this_iter <- vector("list", n_new)
      gammas_this_iter <- vector("list", n_new)
      draw_failed <- FALSE

      for (k in 1:n_new) {
        # Columns are grouped by rank component, so test point k's columns
        # are strided
        idx_k <- k + (0:(r_d - 1)) * n_new
        V_12_k <- V_12_all[, idx_k, drop = FALSE]
        V_21_k <- t(V_12_k)
        V_22_k <- V_22_all[idx_k, idx_k, drop = FALSE]
        Y_k <- Y_all[, idx_k, drop = FALSE]

        mu <- V_21_k %*% y
        sigma <- as.matrix(V_22_k - V_21_k %*% Y_k)

        # in case of negative diagonals due to numerical imprecision
        if (any(diag(sigma) < 0)) {
          diag(sigma)[diag(sigma) < 0] <- 1e-12
        }

        # Handle scalar sigma case
        if (length(sigma) == 1 && isTRUE(as.numeric(sigma) < 0)) {
          draw_failed <- TRUE
          break
        }

        # Check if sigma is positive definite
        eigvals <- tryCatch(
          eigen(forceSymmetric(sigma), symmetric = TRUE, only.values = TRUE)$values,
          error = function(e) NA
        )
        if (any(is.na(eigvals)) || any(eigvals <= 0)) {
          draw_failed <- TRUE
          break
        }

        # Draw from predictive distribution
        gamma_star <- t(rmvnorm(1, mu, sigma))

        gamma_part <- B_eta %*% G_unfold %*% gamma_star
        noise_part <- rnorm(n, 0, sd = sqrt(1 / lambda_eta))

        if (include_noise) {
          pred <- gamma_part + noise_part
        } else {
          pred <- gamma_part
        }

        gammas_this_iter[[k]] <- gamma_star
        draws_this_iter[[k]] <- pred
      }

      if (draw_failed) next

      count <- count + 1
      for (k in 1:n_new) {
        gamma_output_list[[k]][[count]] <- gammas_this_iter[[k]]
        pred_output_list[[k]][[count]] <- draws_this_iter[[k]]
      }
      if (verbose) pb$tick()
    }

    if (count == 0) {
      stop("❌ No valid posterior draws were collected.")
    }

    list(
      pred_output_mats  = lapply(pred_output_list, function(l) do.call(cbind, l)),
      gamma_output_mats = lapply(gamma_output_list, function(l) do.call(cbind, l))
    )
  }

  crps_from_samples <- function(y, draws) {
    # y: vector of observations (length N)
    # draws: matrix [N x S] (S posterior samples per location)

    S <- ncol(draws)

    term1 <- rowMeans(abs(draws - y))

    term2 <- numeric(nrow(draws))
    for (i in seq_len(nrow(draws))) {
      di <- draws[i, ]
      term2[i] <- mean(abs(outer(di, di, "-")))
    }

    term1 - 0.5 * term2
  }

  out <- post_pred_sf(untested_input, n_post_draws)

  n_new <- length(out$pred_output_mats)

  # Per-test-input post-processing
  draws_list <- vector("list", n_new)
  mean_list <- vector("list", n_new)
  sd_list <- vector("list", n_new)

  for (idx in 1:n_new) {
    post_pred_draws_mat <- out$pred_output_mats[[idx]]
    post_pred_draws_arr <- array(post_pred_draws_mat, c(output_dim, n_post_draws))

    # Interpolate posterior preds
    if (!hf) {
      post_pred_draws_arr <- fast_knn_interpolate_first_dim(post_pred_draws_arr, coord_lf, coord_hf, k = 4)
      post_pred_draws_mat <- flatten_last_dim(post_pred_draws_arr)
    }

    spatiotemp_modes <- 1:(length(dim(post_pred_draws_arr)) - 1)

    draws_list[[idx]] <- post_pred_draws_arr
    mean_list[[idx]] <- array(rowMeans(post_pred_draws_mat), dim = dim(post_pred_draws_arr)[spatiotemp_modes])
    sd_list[[idx]] <- array(rowSds(post_pred_draws_mat), dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  }

  spatiotemp_dim <- dim(mean_list[[1]])

  post_pred_draws_arr <- array(unlist(draws_list), dim = c(dim(draws_list[[1]]), n_new))
  post_pred_mean <- array(unlist(mean_list), dim = c(spatiotemp_dim, n_new))
  post_pred_sd <- array(unlist(sd_list), dim = c(spatiotemp_dim, n_new))

  return(list(
    draws = post_pred_draws_arr,
    mean = post_pred_mean,
    sd = post_pred_sd
  ))
}

# --- Multi-fidelity model (LF + discrepancy, jointly sampled) ---

# Prepares matrices and posterior parameters for multifidelity MCMC sampling
prepare_mf_mcmc <- function(output_hf, tucker_outputs_lf, tucker_outputs_discrep, interp_bases_outputs, ranks_lf, ranks_discrep,
                            a_delta, b_delta) {
  # Number of high-fidelity simulations
  d_train_hf <- dim(output_hf)[length(dim(output_hf))]

  # Extract core tensors
  G <- tucker_outputs_lf$core_tensor
  G_discrep <- tucker_outputs_discrep$core_tensor

  # Unfold core tensors along last mode (design) and transpose
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  G_unfold_discrep <- t(k_unfold(as.tensor(G_discrep), m = length(dim(G_discrep)))@data)

  # Extract interpolated LF bases and discrepancy bases
  bases_discrep <- tucker_outputs_discrep$bases

  output_hf_vec <- as.vector(output_hf)

  r_d_lf <- ranks_lf[length(ranks_lf)]
  r_d_discrep <- ranks_discrep[length(ranks_discrep)]

  n_out <- prod(dim(output_hf)[-length(dim(output_hf))])

  if (!is.null(bases_discrep)) {
    # ...
  } else {
    # calculate c_vecs
    c_vecs <- list()

    index_grid_lf <- expand.grid(lapply(ranks_lf[-length(ranks_lf)], function(n) 1:n))
    index_grid_discrep <- expand.grid(lapply(ranks_discrep[-length(ranks_discrep)], function(n) 1:n))

    pb <- progress_bar$new(
      total = r_d_lf * nrow(index_grid_lf) + r_d_discrep * nrow(index_grid_discrep) + (r_d_lf + r_d_discrep) * d_train_hf,
      format = " Calculating Matrices [:bar] :current/:total :percent ETA: :eta",
      clear = TRUE,
      width = 70
    )

    # lf_interp vecs (going to calculate anyways)
    for (k in 1:r_d_lf) {
      out <- rep(0, n_out)
      for (i in 1:nrow(index_grid_lf)) {
        b <- as.vector(calculate_basis(
          U_list = interp_bases_outputs$factor_matrices_interp,
          idx = index_grid_lf[i, ]
        ))

        out <- out + G_unfold[i, k] * b

        pb$tick()

        if (i %% 100 == 0) {
          gc()
        }
      }
      c_vecs[[k]] <- out
    }

    # discrep vecs
    for (k in 1:r_d_discrep) {
      out <- rep(0, n_out)
      for (i in 1:nrow(index_grid_discrep)) {
        b <- as.vector(calculate_basis(
          U_list = tucker_outputs_discrep$factor_matrices,
          idx = index_grid_discrep[i, ]
        ))

        out <- out + G_unfold_discrep[i, k] * b

        pb$tick()

        if (i %% 100 == 0) {
          gc()
        }
      }
      c_vecs[[r_d_lf + k]] <- out
    }

    # Calculate C_coef
    n <- length(c_vecs)
    C_coef <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        C_coef[i, j] <- sum(c_vecs[[j]] * c_vecs[[i]])
      }
    }

    solve_Kt_K <- kronecker(solve(C_coef), diag(d_train_hf))

    Kt_x <- c()
    for (k in 1:length(c_vecs)) {
      for (d in 1:d_train_hf) {
        nd <- length(dim(output_hf))
        z_d <- as.vector(do.call(`[`, c(list(output_hf), rep(list(TRUE), nd - 1), list(d))))
        Kt_x <- c(
          Kt_x,
          t(c_vecs[[k]]) %*% z_d
        )

        pb$tick()
      }
    }
  }

  # Estimate joint coefficient vector (u_hat = [gamma_hat_hf; zeta_hat])
  u_hat <- solve_Kt_K %*% Kt_x

  # Extract LF and discrepancy coefficients from joint estimate
  gamma_hat_hf <- u_hat[1:(r_d_lf * d_train_hf)]
  zeta_hat <- u_hat[-(1:(r_d_lf * d_train_hf))]

  # Update inverse-gamma hyperparameters for variance (discrepancy noise)
  a_delta_prime <- a_delta + (d_train_hf * (n_out - r_d_lf - r_d_discrep)) / 2
  b_delta_prime <- b_delta + 0.5 * (
    sum(output_hf_vec^2) -
      t(Kt_x) %*% solve_Kt_K %*% Kt_x
  )

  b_delta_prime <- as.numeric(b_delta_prime)

  # Return results for MCMC sampling
  return(list(
    solve_K_t_K   = solve_Kt_K,
    gamma_hat_hf  = gamma_hat_hf,
    zeta_hat      = zeta_hat,
    a_delta_prime = a_delta_prime,
    b_delta_prime = b_delta_prime,
    r_d_lf        = r_d_lf,
    r_d_discrep   = r_d_discrep
  ))
}

# Performs multifidelity MCMC sampling for joint posterior estimation
perform_mf_mcmc <- function(
  mcmc_lf_setup_outputs, mcmc_mf_setup_outputs, n_iter, burn_in,
  design_train_lf, design_train_hf,
  init_lambda_eta = NULL,
  init_lambda_w_vec = NULL,
  init_rho_w_mat = NULL,
  init_lambda_delta = NULL,
  init_lambda_v_vec = NULL,
  init_rho_v_mat = NULL,
  omega_init = list(
    lambda_eta = 0.5,
    lambda_w = 0.5,
    rho_w = 1.5,
    lambda_delta = 0.5,
    lambda_v = 1.0,
    rho_v = 0.175
  ),
  learning_rate,
  prior,
  progress_file = NULL
) {
  r_d_lf <- mcmc_mf_setup_outputs$r_d_lf
  r_d_discrep <- mcmc_mf_setup_outputs$r_d_discrep
  solve_C_t_C <- mcmc_lf_setup_outputs$solve_C_t_C
  solve_K_t_K <- mcmc_mf_setup_outputs$solve_K_t_K

  gamma_hat_lf <- mcmc_lf_setup_outputs$gamma_hat
  gamma_hat_hf <- mcmc_mf_setup_outputs$gamma_hat_hf
  zeta_hat <- mcmc_mf_setup_outputs$zeta_hat
  hat_vals <- c(as.vector(gamma_hat_lf), as.vector(gamma_hat_hf), as.vector(zeta_hat))

  a_eta_prime <- mcmc_lf_setup_outputs$a_eta_prime
  b_eta_prime <- mcmc_lf_setup_outputs$b_eta_prime
  a_delta_prime <- mcmc_mf_setup_outputs$a_delta_prime
  b_delta_prime <- mcmc_mf_setup_outputs$b_delta_prime

  # Builds the block-diagonal noise/nugget term
  build_sigma1 <- function(lambda_eta, lambda_delta) {
    Sigma_C <- lambda_eta^(-1) * solve_C_t_C
    Sigma_K <- lambda_delta^(-1) * solve_K_t_K
    Sigma.1 <- matrix(0, nrow = nrow(Sigma_C) + nrow(Sigma_K), ncol = ncol(Sigma_C) + ncol(Sigma_K))
    Sigma.1[1:nrow(Sigma_C), 1:ncol(Sigma_C)] <- as.matrix(Sigma_C)
    Sigma.1[
      (nrow(Sigma_C) + 1):(nrow(Sigma_C) + nrow(Sigma_K)),
      (ncol(Sigma_C) + 1):(ncol(Sigma_C) + ncol(Sigma_K))
    ] <- as.matrix(Sigma_K)
    Sigma.1
  }

  log_posterior_mfm <- function(lambda_eta, lambda_w_vec, rho_w_mat,
                                lambda_delta, lambda_v_vec, rho_v_mat,
                                dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33,
                                Sigma.1) {
    # Build covariances from cached rho-scaled squared distances
    Sigma11 <- Sigma_general_fast(dists_sq_11, lambda_w_vec)
    Sigma22 <- Sigma_general_fast(dists_sq_22, lambda_w_vec)
    Sigma33 <- Sigma_general_fast(dists_sq_33, lambda_v_vec)

    Sigma12 <- Sigma_general_fast(dists_sq_12, lambda_w_vec)
    Sigma21 <- t(Sigma12)

    n1 <- nrow(Sigma11)
    n2 <- nrow(Sigma22)
    n3 <- nrow(Sigma33)

    Sigma.2 <- matrix(0, n1 + n2 + n3, n1 + n2 + n3)

    # Fill blocks
    Sigma.2[1:n1, 1:n1] <- Sigma11
    Sigma.2[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)] <- Sigma22
    Sigma.2[(n1 + n2 + 1):(n1 + n2 + n3), (n1 + n2 + 1):(n1 + n2 + n3)] <- Sigma33

    Sigma.2[1:n1, (n1 + 1):(n1 + n2)] <- Sigma12
    Sigma.2[(n1 + 1):(n1 + n2), 1:n1] <- Sigma21

    Sigma <- Sigma.1 + Sigma.2

    # chol() doubles as the positive-definiteness check
    L <- tryCatch(base::chol(Sigma), error = function(e) NULL)
    if (is.null(L)) {

      Sigma <- as.matrix(forceSymmetric(
        nearPD(Sigma, ensureSymmetry = TRUE)$mat
      ))

      # Add small jitter to diagonal to ensure strict PD
      L <- base::chol(Sigma)
    }

    out.1 <- log_mvnorm_cholesky(hat_vals, rep(0, length(hat_vals)), L = L)

    # Gamma priors times the Tucker-residual likelihoods (the *_prime parameters)
    out.2 <- dgamma(lambda_eta, a_eta_prime, b_eta_prime, log = TRUE)

    out.3 <- sum(dhalfnorm(lambda_w_vec, prior$sigma_lambda_w, log = TRUE))

    out.4 <- sum(dlnorm(rho_w_mat, prior$a_rho_w, prior$b_rho_w, log = TRUE))

    out.5 <- dgamma(lambda_delta, a_delta_prime, b_delta_prime, log = TRUE)

    out.6 <- sum(dhalfnorm(lambda_v_vec, prior$sigma_lambda_v, log = TRUE))

    out.7 <- sum(dlnorm(rho_v_mat, prior$a_rho_v, prior$b_rho_v, log = TRUE))

    return(as.numeric(out.1 + out.2 + out.3 + out.4 + out.5 + out.6 + out.7))
  }

  # log_post_current is the log-posterior at the current state
  update_lambda_eta_MH_MFM <- function(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                       dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33,
                                       omega, log_post_current) {
    # Propose from truncated normal (truncated at 0 from below)
    lambda_eta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_eta, sd = omega)

    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_eta_prop, a = 0, b = Inf, mean = lambda_eta, sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_eta, a = 0, b = Inf, mean = lambda_eta_prop, sd = omega) |> log()

    # Sigma.1 depends on lambda_eta, so it must be rebuilt for this proposal
    # specifically (can't reuse a cached value here) -- the main loop rebuilds
    # its own carried-forward copy once more after this returns
    Sigma1_prop <- build_sigma1(lambda_eta_prop, lambda_delta)
    log_post_prop <- log_posterior_mfm(
      lambda_eta_prop, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_prop
    )

    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = lambda_eta_prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_eta, accept = FALSE, log_post = log_post_current))
    }
  }

  update_lambda_w_vec_i_MH_MFM <- function(i, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                           dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
                                           omega, log_post_current) {
    prop <- lambda_w_vec

    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega)
    prop[i] <- prop_i_new

    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_w_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega) |> log()

    log_post_prop <- log_posterior_mfm(
      lambda_eta, prop, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr
    )

    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_w_vec, accept = FALSE, log_post = log_post_current))
    }
  }

  update_rho_w_mat_ij_MH_MFM <- function(i, j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                         dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
                                         omega_rho, log_post_current) {
    rho_w_pro <- rho_w_mat

    # Propose new rho[i,j] from truncated normal > 0
    rho_w_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_w_mat[i, j], sd = omega_rho[i, j])
    rho_w_pro[i, j] <- rho_w_prop_ij

    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_w_prop_ij, a = 0, b = Inf, mean = rho_w_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_w_mat[i, j], a = 0, b = Inf, mean = rho_w_prop_ij, sd = omega_rho[i, j]) |> log()

    # Only component i's distances changed, so update just that cache element
    dists_sq_11_prop <- update_distance_cache_element_general(design_train_lf, design_train_lf, rho_w_pro, dists_sq_11, i)
    dists_sq_22_prop <- update_distance_cache_element_general(design_train_hf, design_train_hf, rho_w_pro, dists_sq_22, i)
    dists_sq_12_prop <- update_distance_cache_element_general(design_train_lf, design_train_hf, rho_w_pro, dists_sq_12, i)

    log_post_prop <- log_posterior_mfm(
      lambda_eta, lambda_w_vec, rho_w_pro, lambda_delta, lambda_v_vec, rho_v_mat,
      dists_sq_11_prop, dists_sq_22_prop, dists_sq_12_prop, dists_sq_33, Sigma1_curr
    )

    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(
        value = rho_w_pro, dists_sq_11 = dists_sq_11_prop, dists_sq_22 = dists_sq_22_prop,
        dists_sq_12 = dists_sq_12_prop, accept = TRUE, log_post = log_post_prop
      ))
    } else {
      return(list(
        value = rho_w_mat, dists_sq_11 = dists_sq_11, dists_sq_22 = dists_sq_22,
        dists_sq_12 = dists_sq_12, accept = FALSE, log_post = log_post_current
      ))
    }
  }

  update_lambda_delta_MH_MFM <- function(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                         dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33,
                                         omega, log_post_current) {
    # Propose from truncated normal (truncated at 0 from below)
    lambda_delta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_delta, sd = omega)

    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_delta_prop, a = 0, b = Inf, mean = lambda_delta, sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_delta, a = 0, b = Inf, mean = lambda_delta_prop, sd = omega) |> log()

    # Sigma.1 depends on lambda_delta too -- rebuild for this proposal (see
    # update_lambda_eta_MH_MFM's comment above for why).
    Sigma1_prop <- build_sigma1(lambda_eta, lambda_delta_prop)
    log_post_prop <- log_posterior_mfm(
      lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta_prop, lambda_v_vec, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_prop
    )

    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = lambda_delta_prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_delta, accept = FALSE, log_post = log_post_current))
    }
  }

  update_lambda_v_vec_i_MH_MFM <- function(i, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                           dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
                                           omega, log_post_current) {
    prop <- lambda_v_vec

    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_v_vec[i], sd = omega)
    prop[i] <- prop_i_new

    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_v_vec[i], sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_v_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega) |> log()

    log_post_prop <- log_posterior_mfm(
      lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, prop, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr
    )

    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = lambda_v_vec, accept = FALSE, log_post = log_post_current))
    }
  }

  update_rho_v_mat_ij_MH_MFM <- function(i, j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
                                         dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
                                         omega_rho, log_post_current) {
    rho_v_pro <- rho_v_mat

    # Propose new rho[i,j] from truncated normal > 0
    rho_v_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_v_mat[i, j], sd = omega_rho[i, j])
    rho_v_pro[i, j] <- rho_v_prop_ij

    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_v_prop_ij, a = 0, b = Inf, mean = rho_v_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_v_mat[i, j], a = 0, b = Inf, mean = rho_v_prop_ij, sd = omega_rho[i, j]) |> log()

    # Only component i's distances actually changed -- update just that
    # element (Sigma33 is the only cache rho_v_mat feeds).
    dists_sq_33_prop <- update_distance_cache_element_general(design_train_hf, design_train_hf, rho_v_pro, dists_sq_33, i)

    log_post_prop <- log_posterior_mfm(
      lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_pro,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33_prop, Sigma1_curr
    )

    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)

    if (log(runif(1)) < lpr) {
      return(list(value = rho_v_pro, dists_sq_33 = dists_sq_33_prop, accept = TRUE, log_post = log_post_prop))
    } else {
      return(list(value = rho_v_mat, dists_sq_33 = dists_sq_33, accept = FALSE, log_post = log_post_current))
    }
  }

  # Initialize parameters
  lambda_eta <- init_lambda_eta %||% (a_eta_prime / b_eta_prime)
  lambda_w_vec <- init_lambda_w_vec %||% rep(1, r_d_lf)
  rho_w_mat <- init_rho_w_mat %||% matrix(1, r_d_lf, ncol(design_train_lf))

  lambda_delta <- init_lambda_delta %||% (a_delta_prime / b_delta_prime)
  lambda_v_vec <- init_lambda_v_vec %||% rep(1, r_d_discrep)
  rho_v_mat <- init_rho_v_mat %||% matrix(1, r_d_discrep, ncol(design_train_hf))

  # Squared-distance caches (see log_posterior_mfm's comment) -- rho_w_mat feeds
  # THREE of them (design_train_lf x design_train_lf, design_train_hf x
  # design_train_hf, design_train_lf x design_train_hf); rho_v_mat feeds its own
  # (design_train_hf x design_train_hf)
  dists_sq_11 <- initialize_distance_cache_general(design_train_lf, design_train_lf, rho_w_mat)
  dists_sq_22 <- initialize_distance_cache_general(design_train_hf, design_train_hf, rho_w_mat)
  dists_sq_12 <- initialize_distance_cache_general(design_train_lf, design_train_hf, rho_w_mat)
  dists_sq_33 <- initialize_distance_cache_general(design_train_hf, design_train_hf, rho_v_mat)

  # Storage for posterior samples
  n_save <- n_iter - burn_in

  lambda_eta_chain <- numeric(n_save)
  lambda_w_vec_chain <- matrix(NA, nrow = n_save, ncol = r_d_lf)
  rho_w_mat_chain <- array(NA, dim = c(n_save, r_d_lf, ncol(design_train_lf)))

  lambda_delta_chain <- numeric(n_save)
  lambda_v_vec_chain <- matrix(NA, nrow = n_save, ncol = r_d_discrep)
  rho_v_mat_chain <- array(NA, dim = c(n_save, r_d_discrep, ncol(design_train_hf)))

  log_lik_chain <- numeric(n_save)

  # Track acceptance counts and step sizes
  target_accept <- 0.44 # https://support.sas.com/documentation/cdl/en/statug/68162/HTML/default/viewer.htm#statug_bchoice_details17.htm

  accept_lambda_eta <- 0
  omega_lambda_eta <- omega_init$lambda_eta

  accept_lambda_w <- rep(0, r_d_lf)
  omega_lambda_w <- rep(omega_init$lambda_w, r_d_lf)

  accept_rho_w <- matrix(0, r_d_lf, ncol(design_train_lf))
  omega_rho_w <- matrix(omega_init$rho_w, r_d_lf, ncol(design_train_lf))

  accept_lambda_delta <- 0
  omega_lambda_delta <- omega_init$lambda_delta

  accept_lambda_v <- rep(0, r_d_discrep)
  omega_lambda_v <- rep(omega_init$lambda_v, r_d_discrep)

  accept_rho_v <- matrix(0, r_d_discrep, ncol(design_train_hf))
  omega_rho_v <- matrix(omega_init$rho_v, r_d_discrep, ncol(design_train_hf))

  update_omega_iter <- 50

  # Progress bar
  pb <- progress_bar$new(
    format = "  MF MCMC [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
    total = n_iter,
    width = 70
  )

  # Current log-posterior and Sigma.1, carried forward across MH updates
  Sigma1_curr <- build_sigma1(lambda_eta, lambda_delta)
  log_post_curr <- log_posterior_mfm(
    lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
    dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr
  )

  for (iter in 1:n_iter) {
    # Update lambda_eta
    res <- update_lambda_eta_MH_MFM(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33,
      omega_lambda_eta,
      log_post_current = log_post_curr
    )
    lambda_eta <- res$value
    log_post_curr <- res$log_post
    # Sigma.1 depends on lambda_eta, so rebuild it after each lambda_eta update
    Sigma1_curr <- build_sigma1(lambda_eta, lambda_delta)
    # Update omega during burn-in
    if (iter <= burn_in) {
      accept_lambda_eta <- accept_lambda_eta + res$accept
      # Adapt every few steps
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_eta / update_omega_iter # only look at local acceptance
        omega_lambda_eta <- omega_lambda_eta * exp((rate - target_accept) / sqrt(iter) * learning_rate) # https://mvihola.github.io/docs/AdaptiveMCMC.jl/adapt/?utm_source=chatgpt.com
        accept_lambda_eta <- 0
      }
    }

    # Update lambda_w_vec
    for (j in 1:r_d_lf) {
      res <- update_lambda_w_vec_i_MH_MFM(j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
        dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
        omega_lambda_w[j],
        log_post_current = log_post_curr
      )
      lambda_w_vec <- res$value
      log_post_curr <- res$log_post
      if (iter <= burn_in) {
        accept_lambda_w[j] <- accept_lambda_w[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_w[j] / update_omega_iter
          omega_lambda_w[j] <- omega_lambda_w[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_w[j] <- 0
        }
      }
    }

    # Update rho_w_mat
    for (i_rho in 1:r_d_lf) {
      for (j_rho in 1:ncol(design_train_lf)) {
        res <- update_rho_w_mat_ij_MH_MFM(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
          dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
          omega_rho_w,
          log_post_current = log_post_curr
        )
        rho_w_mat <- res$value
        dists_sq_11 <- res$dists_sq_11
        dists_sq_22 <- res$dists_sq_22
        dists_sq_12 <- res$dists_sq_12
        log_post_curr <- res$log_post
        if (iter <= burn_in) {
          accept_rho_w[i_rho, j_rho] <- accept_rho_w[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho_w[i_rho, j_rho] / update_omega_iter
            omega_rho_w[i_rho, j_rho] <- omega_rho_w[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho_w[i_rho, j_rho] <- 0
          }
        }
      }
    }

    # Update lambda_delta
    res <- update_lambda_delta_MH_MFM(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
      dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33,
      omega_lambda_delta,
      log_post_current = log_post_curr
    )
    lambda_delta <- res$value
    log_post_curr <- res$log_post
    # Sigma.1 also depends on lambda_delta -- rebuild once here for the
    # lambda_v_vec/rho_v_mat sweep below, same reasoning as after lambda_eta.
    Sigma1_curr <- build_sigma1(lambda_eta, lambda_delta)
    # Update omega during burn-in
    if (iter <= burn_in) {
      accept_lambda_delta <- accept_lambda_delta + res$accept
      # Adapt every few steps
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_delta / update_omega_iter # only look at local acceptance
        omega_lambda_delta <- omega_lambda_delta * exp((rate - target_accept) / sqrt(iter) * learning_rate) # https://mvihola.github.io/docs/AdaptiveMCMC.jl/adapt/?utm_source=chatgpt.com
        accept_lambda_delta <- 0
      }
    }

    # Update lambda_v_vec
    for (j in 1:r_d_discrep) {
      res <- update_lambda_v_vec_i_MH_MFM(j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
        dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
        omega_lambda_v[j],
        log_post_current = log_post_curr
      )
      lambda_v_vec <- res$value
      log_post_curr <- res$log_post
      if (iter <= burn_in) {
        accept_lambda_v[j] <- accept_lambda_v[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_v[j] / update_omega_iter
          omega_lambda_v[j] <- omega_lambda_v[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_v[j] <- 0
        }
      }
    }

    # Update rho_v_mat
    for (i_rho in 1:r_d_discrep) {
      for (j_rho in 1:ncol(design_train_lf)) {
        res <- update_rho_v_mat_ij_MH_MFM(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat,
          dists_sq_11, dists_sq_22, dists_sq_12, dists_sq_33, Sigma1_curr,
          omega_rho_v,
          log_post_current = log_post_curr
        )
        rho_v_mat <- res$value
        dists_sq_33 <- res$dists_sq_33
        log_post_curr <- res$log_post
        if (iter <= burn_in) {
          accept_rho_v[i_rho, j_rho] <- accept_rho_v[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho_v[i_rho, j_rho] / update_omega_iter
            omega_rho_v[i_rho, j_rho] <- omega_rho_v[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho_v[i_rho, j_rho] <- 0
          }
        }
      }
    }

    # Store the results after the burn-in period
    if (iter > burn_in) {
      lambda_eta_chain[iter - burn_in] <- lambda_eta
      lambda_w_vec_chain[iter - burn_in, ] <- lambda_w_vec
      rho_w_mat_chain[iter - burn_in, , ] <- rho_w_mat

      lambda_delta_chain[iter - burn_in] <- lambda_delta
      lambda_v_vec_chain[iter - burn_in, ] <- lambda_v_vec
      rho_v_mat_chain[iter - burn_in, , ] <- rho_v_mat

      log_lik_chain[iter - burn_in] <- log_post_curr
    }

    if (!is.null(progress_file)) cat(iter, file = progress_file)

    pb$tick()
  }

  return(list(
    n_iter = n_iter,
    burn_in = burn_in,
    init_params = list(
      lambda_eta_init = init_lambda_eta,
      lambda_w_vec_init = init_lambda_w_vec,
      rho_w_mat_init = init_rho_w_mat,
      lambda_delta_init = init_lambda_delta,
      lambda_v_vec_init = init_lambda_v_vec,
      rho_v_mat_init = init_rho_v_mat
    ),
    post_chains = list(
      lambda_eta_chain = lambda_eta_chain,
      lambda_w_vec_chain = lambda_w_vec_chain,
      rho_w_mat_chain = rho_w_mat_chain,
      lambda_delta_chain = lambda_delta_chain,
      lambda_v_vec_chain = lambda_v_vec_chain,
      rho_v_mat_chain = rho_v_mat_chain,
      log_lik_chain = log_lik_chain
    )
  ))
}

# Posterior predictive mean/sd (and draws) of the MF model at new inputs
evaluate_mf_mcmc <- function(untested_input,
                             tucker_outputs_lf, tucker_outputs_discrep,
                             interpolate_bases_outputs,
                             mcmc_lf_setup_outputs, mcmc_mf_setup_outputs,
                             mf_mcmc_outputs,
                             n_post_draws,
                             output_dim,
                             design_train_lf, design_train_hf,
                             include_noise = FALSE,
                             obs = NULL,
                             verbose = TRUE,
                             r_d_lf,
                             r_d_discrep,
                             D_G = NULL) {
  tuck_discrep <- tucker_outputs_discrep

  n_chains <- length(mf_mcmc_outputs)

  solve_C_t_C <- mcmc_lf_setup_outputs$solve_C_t_C
  solve_K_t_K <- mcmc_mf_setup_outputs$solve_K_t_K

  gamma_hat_lf <- mcmc_lf_setup_outputs$gamma_hat
  gamma_hat_hf <- mcmc_mf_setup_outputs$gamma_hat_hf
  zeta_hat <- mcmc_mf_setup_outputs$zeta_hat
  hat_vals <- c(as.vector(gamma_hat_lf), as.vector(gamma_hat_hf), as.vector(zeta_hat))

  G <- tucker_outputs_lf$core_tensor
  G_discrep <- tucker_outputs_discrep$core_tensor

  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  G_unfold_discrep <- t(k_unfold(as.tensor(G_discrep), m = length(dim(G_discrep)))@data)

  # Suppose tuck_discrep$factor_matrices has N modes
  ranks <- sapply(tucker_outputs_discrep$factor_matrices, ncol) # number of ranks per mode

  # Generate all combinations of indices
  idx_combinations <- expand.grid(lapply(ranks[-length(ranks)], seq_len))

  if (is.null(D_G)) {
    D_G <- matrix(0, nrow = prod(output_dim), ncol = ncol(G_unfold_discrep))

    if (verbose) {
      pb <- progress_bar$new(
        format = "  Matrix Calc [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
        total = prod(ranks),
        width = 70
      )
    }

    for (counter in seq_len(nrow(idx_combinations))) {
      idx <- as.numeric(idx_combinations[counter, ])
      basis_matrix <- calculate_basis(U_list = tuck_discrep$factor_matrices, idx = idx)

      D_G <- D_G + as.vector(basis_matrix) %*% G_unfold_discrep[counter, , drop = FALSE]

      if (verbose) pb$tick()
    }

  }

  # Prepare chain info for sampling
  chain_lengths <- sapply(mf_mcmc_outputs, function(chain) {
    length(chain$post_chains$lambda_eta_chain)
  })

  # Prediction function; A+B's Cholesky is computed once per draw for all test points
  post_pred_mf <- function(theta_star, n_post_draws) {
    theta_star <- as.matrix(theta_star)
    n_new <- nrow(theta_star)

    # Raw (rho-independent) squared distances, computed ONCE and reused for
    # every posterior draw's covariance construction below (see
    # Sigma_general_from_raw()'s comment)
    raw_lf_lf <- initialize_raw_dist_cache(design_train_lf, design_train_lf)
    raw_lf_hf <- initialize_raw_dist_cache(design_train_lf, design_train_hf)
    raw_hf_hf <- initialize_raw_dist_cache(design_train_hf, design_train_hf)
    raw_lf_theta <- initialize_raw_dist_cache(design_train_lf, theta_star)
    raw_hf_theta <- initialize_raw_dist_cache(design_train_hf, theta_star)
    raw_theta_theta <- initialize_raw_dist_cache(theta_star, theta_star)

    col_width <- r_d_lf + r_d_discrep
    pred_output_list <- vector("list", n_new)
    for (k in 1:n_new) pred_output_list[[k]] <- list()

    # Parallel per-draw storage for the LF-only (gamma_part) and discrepancy-
    # only (zeta_part [+ noise_part]) components
    pred_output_list_lf <- vector("list", n_new)
    pred_output_list_discrep <- vector("list", n_new)
    for (k in 1:n_new) {
      pred_output_list_lf[[k]] <- list()
      pred_output_list_discrep[[k]] <- list()
    }

    # Raw per-rank-component coefficient draws (gamma_star/zeta_star), BEFORE
    # the spatial basis expansion (gamma_part = B_eta_tile_x_G_unfold %*%
    # gamma_star) that turns them into a field
    pred_output_list_gamma <- vector("list", n_new)
    pred_output_list_zeta <- vector("list", n_new)
    for (k in 1:n_new) {
      pred_output_list_gamma[[k]] <- list()
      pred_output_list_zeta[[k]] <- list()
    }

    count <- 0
    attempts <- 0
    max_attempts <- 3 * n_post_draws # Avoid infinite loops

    if (verbose) {
      pb <- progress_bar$new(
        format = "  MF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
        total = n_post_draws,
        width = 70
      )
    }

    while (count < n_post_draws && attempts < max_attempts) {
      attempts <- attempts + 1

      # Sample chain and draw index
      chain_i <- sample(1:n_chains, 1)
      draw_i <- sample(1:chain_lengths[chain_i], 1)
      chain <- mf_mcmc_outputs[[chain_i]]$post_chains

      # Extract parameters
      lambda_eta <- chain$lambda_eta_chain[draw_i]
      lambda_w_vec <- chain$lambda_w_vec_chain[draw_i, ]
      rho_w_mat <- drop_first_dim(chain$rho_w_mat_chain[draw_i, , , drop = FALSE])

      lambda_delta <- chain$lambda_delta_chain[draw_i]
      lambda_v_vec <- chain$lambda_v_vec_chain[draw_i, ]
      rho_v_mat <- drop_first_dim(chain$rho_v_mat_chain[draw_i, , , drop = FALSE])

      # Construct covariance blocks
      V_11 <- Sigma_general_from_raw(raw_lf_lf, lambda_w_vec, rho_w_mat, r_d_lf)
      V_12 <- Sigma_general_from_raw(raw_lf_hf, lambda_w_vec, rho_w_mat, r_d_lf)
      V_21 <- t(V_12)
      V_22 <- Sigma_general_from_raw(raw_hf_hf, lambda_w_vec, rho_w_mat, r_d_lf)
      V <- rbind(cbind(V_11, V_12), cbind(V_21, V_22))

      A <- bdiag(
        V,
        Sigma_general_from_raw(raw_hf_hf, lambda_v_vec, rho_v_mat, r_d_discrep)
      )
      B <- bdiag(1 / lambda_eta * solve_C_t_C, 1 / lambda_delta * solve_K_t_K)

      L <- tryCatch(chol(A + B), error = function(e) NULL)
      if (is.null(L)) {
        next
      }

      y <- backsolve(L, forwardsolve(t(L), hat_vals)) # test-input-independent

      # All n_new test points' cross-covariance blocks built in ONE shot per
      # piece, then reassembled per test point (bdiag() is cheap here -- it's
      # just zero-padding already-computed small blocks, not touching the
      # expensive A+B factorization above)
      C_lf_all <- Sigma_general_from_raw(raw_lf_theta, lambda_w_vec, rho_w_mat, r_d_lf)
      C_hf_all <- Sigma_general_from_raw(raw_hf_theta, lambda_w_vec, rho_w_mat, r_d_lf)
      C_discrep_all <- Sigma_general_from_raw(raw_hf_theta, lambda_v_vec, rho_v_mat, r_d_discrep)
      D_lf_all <- Sigma_general_from_raw(raw_theta_theta, lambda_w_vec, rho_w_mat, r_d_lf)
      D_discrep_all <- Sigma_general_from_raw(raw_theta_theta, lambda_v_vec, rho_v_mat, r_d_discrep)

      # Sigma_general_from_raw()'s block structure is grouped by RANK COMPONENT
      # (r_d_lf or r_d_discrep blocks of n_new columns each), not by test point
      C_all <- do.call(cbind, lapply(1:n_new, function(k) {
        idx_lf <- k + (0:(r_d_lf - 1)) * n_new
        idx_dc <- k + (0:(r_d_discrep - 1)) * n_new
        bdiag(
          rbind(C_lf_all[, idx_lf, drop = FALSE], C_hf_all[, idx_lf, drop = FALSE]),
          C_discrep_all[, idx_dc, drop = FALSE]
        )
      }))
      Y_all <- backsolve(L, forwardsolve(t(L), C_all))

      D_list <- lapply(1:n_new, function(k) {
        idx_lf <- k + (0:(r_d_lf - 1)) * n_new
        idx_dc <- k + (0:(r_d_discrep - 1)) * n_new
        bdiag(D_lf_all[idx_lf, idx_lf, drop = FALSE], D_discrep_all[idx_dc, idx_dc, drop = FALSE])
      })

      # A failed draw at any test point skips this posterior draw
      draws_this_iter <- vector("list", n_new)
      lf_draws_this_iter <- vector("list", n_new)
      discrep_draws_this_iter <- vector("list", n_new)
      gamma_draws_this_iter <- vector("list", n_new)
      zeta_draws_this_iter <- vector("list", n_new)
      draw_failed <- FALSE

      for (k in 1:n_new) {
        idx_k <- ((k - 1) * col_width + 1):(k * col_width)
        C_k <- C_all[, idx_k, drop = FALSE]
        Y_k <- Y_all[, idx_k, drop = FALSE]
        D_k <- D_list[[k]]

        mu <- t(C_k) %*% y
        sigma <- as.matrix(D_k - t(C_k) %*% Y_k)

        # Check PD
        eigvals <- eigen(forceSymmetric(sigma), symmetric = TRUE, only.values = TRUE)$values
        if (any(eigvals <= 0)) {
          draw_failed <- TRUE
          break
        }

        out <- rmvnorm(1, mu, sigma)
        gamma_star <- out[1:r_d_lf]
        zeta_star <- out[-(1:r_d_lf)]

        gamma_part <- interpolate_bases_outputs$B_eta_tile_x_G_unfold %*% gamma_star
        if (r_d_discrep == 1) {
          zeta_part <- D_G * zeta_star
        } else {
          zeta_part <- D_G %*% zeta_star
        }
        noise_part <- rnorm(prod(output_dim), 0, sd = sqrt(1 / lambda_delta)) +
          rnorm(prod(output_dim), 0, sd = sqrt(1 / lambda_eta))

        pred <- gamma_part + zeta_part
        if (include_noise) {
          pred <- pred + noise_part
        }

        # LF and discrepancy components of the draw
        lf_draws_this_iter[[k]] <- gamma_part
        discrep_draws_this_iter[[k]] <- if (include_noise) zeta_part + noise_part else zeta_part
        gamma_draws_this_iter[[k]] <- gamma_star
        zeta_draws_this_iter[[k]] <- zeta_star

        draws_this_iter[[k]] <- pred
      }

      if (draw_failed) next

      count <- count + 1
      for (k in 1:n_new) {
        pred_output_list[[k]][[count]] <- draws_this_iter[[k]]
        pred_output_list_lf[[k]][[count]] <- lf_draws_this_iter[[k]]
        pred_output_list_discrep[[k]][[count]] <- discrep_draws_this_iter[[k]]
        pred_output_list_gamma[[k]][[count]] <- gamma_draws_this_iter[[k]]
        pred_output_list_zeta[[k]][[count]] <- zeta_draws_this_iter[[k]]
      }
      if (verbose) pb$tick()
    }

    if (count == 0) {
      stop("❌ No valid posterior draws were collected.")
    }

    list(
      combined = lapply(pred_output_list, function(l) do.call(cbind, l)),
      lf = lapply(pred_output_list_lf, function(l) do.call(cbind, l)),
      discrepancy = lapply(pred_output_list_discrep, function(l) do.call(cbind, l)),
      # [n_post_draws x r_d_lf] / [n_post_draws x r_d_discrep] per test point --
      # rbind (not cbind) since each draw contributes one ROW of per-component
      # coefficients, unlike the spatial fields above
      gamma = lapply(pred_output_list_gamma, function(l) do.call(rbind, l)),
      zeta = lapply(pred_output_list_zeta, function(l) do.call(rbind, l))
    )
  }

  crps_from_samples <- function(y, draws) {
    # y: vector of observations (length N)
    # draws: matrix [N x S] (S posterior samples per location)

    S <- ncol(draws)

    term1 <- rowMeans(abs(draws - y))

    term2 <- numeric(nrow(draws))
    for (i in seq_len(nrow(draws))) {
      di <- draws[i, ]
      term2[i] <- mean(abs(outer(di, di, "-")))
    }

    term1 - 0.5 * term2
  }

  # Run post pred
  post_pred_result <- post_pred_mf(untested_input, n_post_draws)
  post_pred_draws_mats <- post_pred_result$combined
  post_pred_lf_mats <- post_pred_result$lf
  post_pred_discrep_mats <- post_pred_result$discrepancy
  post_pred_gamma_mats <- post_pred_result$gamma
  post_pred_zeta_mats <- post_pred_result$zeta
  n_new <- length(post_pred_draws_mats)

  # Per-test-input mean/sd, stacked along a trailing test-input dimension
  draws_list <- vector("list", n_new)
  mean_list <- vector("list", n_new)
  sd_list <- vector("list", n_new)
  mean_lf_list <- vector("list", n_new)
  sd_lf_list <- vector("list", n_new)
  mean_discrep_list <- vector("list", n_new)
  sd_discrep_list <- vector("list", n_new)
  # Per test point: a length-r_d_lf (gamma) / length-r_d_discrep (zeta) vector
  # of posterior mean/sd, one entry per rank component -- stacked into [n_new x
  # r_d_lf]/[n_new x r_d_discrep] matrices below, unlike the spatiotemporal-
  # array outputs above
  gamma_mean_list <- vector("list", n_new)
  gamma_sd_list <- vector("list", n_new)
  zeta_mean_list <- vector("list", n_new)
  zeta_sd_list <- vector("list", n_new)

  for (idx in 1:n_new) {
    post_pred_draws_mat <- post_pred_draws_mats[[idx]]
    if (inherits(post_pred_draws_mat, "dgeMatrix")) {
      post_pred_draws_mat <- as.matrix(post_pred_draws_mat)
    }
    post_pred_lf_mat <- as.matrix(post_pred_lf_mats[[idx]])
    post_pred_discrep_mat <- as.matrix(post_pred_discrep_mats[[idx]])
    post_pred_gamma_mat <- as.matrix(post_pred_gamma_mats[[idx]])
    post_pred_zeta_mat <- as.matrix(post_pred_zeta_mats[[idx]])

    post_pred_draws_arr <- array(post_pred_draws_mat, c(output_dim, n_post_draws))

    spatiotemp_modes <- 1:(length(dim(post_pred_draws_arr)) - 1)
    spatiotemp_dim_idx <- dim(post_pred_draws_arr)[spatiotemp_modes]

    draws_list[[idx]] <- post_pred_draws_arr
    mean_list[[idx]] <- array(rowMeans(post_pred_draws_mat), dim = spatiotemp_dim_idx)
    sd_list[[idx]] <- array(rowSds(post_pred_draws_mat), dim = spatiotemp_dim_idx)
    mean_lf_list[[idx]] <- array(rowMeans(post_pred_lf_mat), dim = spatiotemp_dim_idx)
    sd_lf_list[[idx]] <- array(rowSds(post_pred_lf_mat), dim = spatiotemp_dim_idx)
    mean_discrep_list[[idx]] <- array(rowMeans(post_pred_discrep_mat), dim = spatiotemp_dim_idx)
    sd_discrep_list[[idx]] <- array(rowSds(post_pred_discrep_mat), dim = spatiotemp_dim_idx)
    gamma_mean_list[[idx]] <- colMeans(post_pred_gamma_mat)
    gamma_sd_list[[idx]] <- colSds(post_pred_gamma_mat)
    zeta_mean_list[[idx]] <- colMeans(post_pred_zeta_mat)
    zeta_sd_list[[idx]] <- colSds(post_pred_zeta_mat)
  }

  spatiotemp_dim <- dim(mean_list[[1]])

  post_pred_draws_arr <- array(unlist(draws_list), dim = c(dim(draws_list[[1]]), n_new))
  post_pred_mean <- array(unlist(mean_list), dim = c(spatiotemp_dim, n_new))
  post_pred_sd <- array(unlist(sd_list), dim = c(spatiotemp_dim, n_new))
  post_pred_lf_mean <- array(unlist(mean_lf_list), dim = c(spatiotemp_dim, n_new))
  post_pred_lf_sd <- array(unlist(sd_lf_list), dim = c(spatiotemp_dim, n_new))
  post_pred_discrep_mean <- array(unlist(mean_discrep_list), dim = c(spatiotemp_dim, n_new))
  post_pred_discrep_sd <- array(unlist(sd_discrep_list), dim = c(spatiotemp_dim, n_new))

  # [n_new x r_d_lf] / [n_new x r_d_discrep] -- one column per rank component,
  # one row per test input
  gamma_mean <- do.call(rbind, gamma_mean_list)
  gamma_sd <- do.call(rbind, gamma_sd_list)
  zeta_mean <- do.call(rbind, zeta_mean_list)
  zeta_sd <- do.call(rbind, zeta_sd_list)

  return(list(
    draws = post_pred_draws_arr,
    mean = post_pred_mean,
    sd = post_pred_sd,
    mean_lf = post_pred_lf_mean,
    sd_lf = post_pred_lf_sd,
    mean_discrepancy = post_pred_discrep_mean,
    sd_discrepancy = post_pred_discrep_sd,
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    zeta_mean = zeta_mean,
    zeta_sd = zeta_sd,
    D_G = D_G
  ))
}

# Format a whole number of seconds as e.g. "45s", "3m 12s", "1h 05m".
format_eta_secs <- function(secs) {
  secs <- max(0, round(secs))
  if (secs < 60) {
    return(sprintf("%ds", secs))
  }
  mins <- secs %/% 60
  rem_s <- secs %% 60
  if (mins < 60) {
    return(sprintf("%dm %02ds", mins, rem_s))
  }
  hrs <- mins %/% 60
  rem_m <- mins %% 60
  sprintf("%dh %02dm", hrs, rem_m)
}

# Runs independent MCMC chains in parallel
run_chains_parallel <- function(chain_fn, n_chains, n_iter_total,
                                n_cores = min(n_chains, max(1, parallel::detectCores() - 1)),
                                poll_interval = 1) {
  progress_files <- vapply(seq_len(n_chains), function(i) tempfile(), character(1))
  on.exit(unlink(progress_files), add = TRUE)

  results <- vector("list", n_chains)
  queued <- seq_len(n_chains)
  pending <- list()
  start_time <- vector("list", n_chains)

  launch <- function(i) {
    job <- parallel::mcparallel(chain_fn(i, progress_files[i]), name = as.character(i), silent = TRUE)
    pending[[as.character(i)]] <<- job
    start_time[[i]] <<- Sys.time()
  }

  # Redraws the block in place: move the cursor up over the previous block, then
  # clear-and-rewrite each line (robust even if a line gets shorter)
  in_rstudio <- Sys.getenv("RSTUDIO") == "1"
  n_lines_shown <- 0
  prev_len <- 0 # RStudio path only: since \033[K isn't supported either, a
  # shorter line must be padded to at least the previous
  # line's length or it'll leave stray trailing characters.
  redraw <- function(lines) {
    if (in_rstudio) {
      chain_status <- paste(trimws(lines[-1]), collapse = " | ")
      content <- paste0(format(Sys.time(), "%H:%M:%S"), " ", chain_status)
      pad_len <- max(nchar(content), prev_len)
      cat("\r", formatC(content, width = -pad_len), sep = "")
      flush.console()
      prev_len <<- pad_len
    } else {
      if (n_lines_shown > 0) cat(sprintf("\033[%dA", n_lines_shown))
      for (ln in lines) cat("\r\033[K", ln, "\n", sep = "")
      n_lines_shown <<- length(lines)
    }
  }

  n_start <- min(n_cores, length(queued))
  for (k in seq_len(n_start)) launch(queued[k])
  queued <- queued[-seq_len(n_start)]

  repeat {
    Sys.sleep(poll_interval)

    status <- vapply(seq_len(n_chains), function(i) {
      if (file.exists(progress_files[i])) {
        v <- suppressWarnings(as.integer(readLines(progress_files[i], warn = FALSE)))
        if (length(v) == 0 || is.na(v)) 0L else v
      } else {
        0L
      }
    }, integer(1))

    lines <- sprintf("[%s] MCMC progress:", format(Sys.time(), "%H:%M:%S"))
    for (i in seq_len(n_chains)) {
      pct <- round(100 * status[i] / n_iter_total)
      eta_str <- "--"
      if (!is.null(start_time[[i]]) && status[i] > 0) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time[[i]], units = "secs"))
        remaining <- n_iter_total - status[i]
        if (remaining <= 0) {
          eta_str <- "0s"
        } else if (elapsed > 0) {
          eta_str <- format_eta_secs(remaining / (status[i] / elapsed))
        }
      }
      lines <- c(lines, sprintf("  Chain %d: %d/%d (%d%%) ETA: %s", i, status[i], n_iter_total, pct, eta_str))
    }
    redraw(lines)

    if (length(pending) > 0) {
      finished <- parallel::mccollect(pending, wait = FALSE)
      if (!is.null(finished)) {
        for (nm in names(finished)) {
          results[[as.integer(nm)]] <- finished[[nm]]
        }
        pending[names(finished)] <- NULL

        n_free <- n_cores - length(pending)
        while (n_free > 0 && length(queued) > 0) {
          launch(queued[1])
          queued <- queued[-1]
          n_free <- n_free - 1
        }
      }
    }

    if (length(pending) == 0 && length(queued) == 0) break
  }

  results
}
