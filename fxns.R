setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(progress)
library(rTensor)
library(Matrix)
library(mvtnorm)
library(FNN)
library(matrixStats)
library(pbmcapply)
library(coda)
library(truncnorm)
library(parallel)
library(pbapply)

theta2unif<- function(A){
  library(KScorrect) # loguniform distribution
  A <- as.matrix(A)
  B <- matrix(NA,nr = nrow(A),nc = ncol(A))    
  
  
  # need to transform based on the prior
  # ref https://www.johndcook.com/quantiles_parameters.pdf to find parameters based on quantiles
  # set min = 1%, max = 99%
  p2 = 0.99; p1 = 0.01
  # log normal for dragio
  min <- log(0.2e-3); max <- log(160e-3)
  sigma <- (max-min)/(qnorm(p2)-qnorm(p1))
  mu <- (min*qnorm(p2)-max*qnorm(p1))/(qnorm(p2)-qnorm(p1))
  # check qlnorm(0.01,mu,sigma) == 0.2e-3 &  qlnorm(0.99,mu,sigma) == 160e-3
  # check plnorm(0.2e-3,mu,sigma) == 0.05 & plnorm(160e-3,mu,sigma) == 0.95
  
  # log u for lambda_pond
  a <- log(1/1.15e-4); b <- log(1/1.15e-8); p2 = 0.99; p1 = 0.01
  r <- p2/p1; u1 <- (b-a*r)/(1-r); u2 <- u1 + (b-a)/(p2-p1)
  
  # check
  # plunif(1/1.15e-4,exp(u1),exp(u2)) == 0.01
  # > plunif(1/1.15e-8,exp(u1),exp(u2)) 
  # [1] 0.99
  # qlunif(p1,exp(u1),exp(u2)) == 1/1.15e-4
  # qlunif(p2,exp(u1),exp(u2)) == 1/1.15e-8
  
  # n distribution for rsnw_mlt
  min <- 250; max <- 3000
  sig <- (max-min)/(qnorm(p2)-qnorm(p1))
  m <- (min*qnorm(p2)-max*qnorm(p1))/(qnorm(p2)-qnorm(p1))
  # check qnorm(p1,m,sig) == 250 &  qnorm(p2,m,sig) == 3000
  # check pnorm(250,m,sig) == p1 & pnorm(3000,m,sig) == p2
  
  B[,1] <- punif(A[,1], min = 0.03, max = 0.65)
  B[,2] <- plnorm(A[,2], mean=mu, sd=sigma) # qunif(A[,2], min = 0.2e-3, max = 160e-3)
  B[,3] <- plunif(A[,3], min = exp(u1), max = exp(u2)) # qunif(A[,3], min = 1/1.15e-4, max = 1/1.15e-8)
  B[,4] <- pnorm(A[,4], mean=m, sd=sig) # qunif(A[,4], min = 250, max = 3000)
  B[,5] <- punif(A[,5], min = -2, max = 2)
  B[,6] <- punif(A[,6], min = 0.55, max = 0.95)
  
  colnames(B) <- c("ksno", "dragio", "lambda_pond", "rsnw_mlt",
                   "R_snw", "phi_i_mushy")
  return(B)
}


unif2theta <- function(A){
  library(KScorrect) # loguniform distribution
  A <- as.matrix(A, nc = 6)
  B <- matrix(NA,nr = nrow(A),nc = ncol(A))    
  
  
  # need to transform based on dragio prior
  # ref https://www.johndcook.com/quantiles_parameters.pdf to find parameters based on quantiles
  # set min = 1%, max = 99%
  p2 = 0.99; p1 = 0.01
  # log normal for dragio
  min <- log(0.2e-3); max <- log(160e-3)
  sigma <- (max-min)/(qnorm(p2)-qnorm(p1))
  mu <- (min*qnorm(p2)-max*qnorm(p1))/(qnorm(p2)-qnorm(p1))
  # check qlnorm(0.01,mu,sigma) == 0.2e-3 &  qlnorm(0.99,mu,sigma) == 160e-3
  # check plnorm(0.2e-3,mu,sigma) == 0.05 & plnorm(160e-3,mu,sigma) == 0.95
  
  # log u for lambda_pond
  a <- log(1/1.15e-4); b <- log(1/1.15e-8); p2 = 0.99; p1 = 0.01
  r <- p2/p1; u1 <- (b-a*r)/(1-r); u2 <- u1 + (b-a)/(p2-p1)
  
  # check
  # plunif(1/1.15e-4,exp(u1),exp(u2)) == 0.01
  # > plunif(1/1.15e-8,exp(u1),exp(u2)) 
  # [1] 0.99
  # qlunif(p1,exp(u1),exp(u2)) == 1/1.15e-4
  # qlunif(p2,exp(u1),exp(u2)) == 1/1.15e-8
  
  # n distribution for rsnw_mlt
  min <- 250; max <- 3000
  sig <- (max-min)/(qnorm(p2)-qnorm(p1))
  m <- (min*qnorm(p2)-max*qnorm(p1))/(qnorm(p2)-qnorm(p1))
  # check qnorm(p1,m,sig) == 250 &  qnorm(p2,m,sig) == 3000
  # check pnorm(250,m,sig) == p1 & pnorm(3000,m,sig) == p2
  
  B[,1] <- qunif(A[,1], min = 0.03, max = 0.65)
  B[,2] <- qlnorm(A[,2], mean=mu, sd=sigma) # qunif(A[,2], min = 0.2e-3, max = 160e-3)
  B[,3] <- qlunif(A[,3], min = exp(u1), max = exp(u2)) # qunif(A[,3], min = 1/1.15e-4, max = 1/1.15e-8)
  B[,4] <- qnorm(A[,4], mean=m, sd=sig) # qunif(A[,4], min = 250, max = 3000)
  B[,5] <- qunif(A[,5], min = -2, max = 2)
  B[,6] <- qunif(A[,6], min = 0.55, max = 0.95)
  
  colnames(B) <- c("ksno", "dragio", "lambda_pond", "rsnw_mlt",
                   "R_snw", "phi_i_mushy")
  return(B)
}


logit <- function(p, center = 0, rate = 1) {
  # if (p <= 0 | p >= 1) stop("Input values must be between 0 and 1")
  (1/rate) * log(p / (1 - p)) + center
}

logistic <- function(x, center = 0, rate = 1) {
  1 / (1 + exp(-1 * rate * (x - center)))
}

softclip <- function(y, epsilon = 0.005) {
  numerator <- exp(y / epsilon) - 1
  denominator <- exp(1 / epsilon) - 1
  numerator / denominator
}

logit_softclip <- function(y, epsilon = 0.005, center = 0.5) {
  sc <- softclip(y, epsilon)
  out <- logit(sc)
  
  out <- out - log(softclip(center, epsilon) / (1 - softclip(center, epsilon)))
  
  return(out)
}

logistic_softclip <- function(logit_y, epsilon = 0.005, center = 0.5) {
  out <- logit_y + log(softclip(center, epsilon) / (1 - softclip(center, epsilon)))
  out <- logistic(out)
  
  # Invert the softclip function
  denom <- exp(1 / epsilon) - 1
  out <- epsilon * log(out * denom + 1)
  
  return(out)
}

scale_transformed <- function(transformed_data, max_data) {
  if (max(transformed_data) != -1 * min(transformed_data)) {
    return("Error: transofrmed data not symmetric.")
  }
  return(transformed_data / max_data / 2 + 0.5)
}

unscaled_transformed <- function(dat, max_data) {
  return(max_data  * 2 * (dat - 0.5))
}


select_rank_variance <- function(tnsr, target_var = 0.95) {
  if (length(target_var) == 1) {
    target_var <- rep(target_var, length(dim(tnsr)))
  }
  
  modes <- tnsr@modes
  ranks <- integer(length(modes))
  
  pb <- progress_bar$new(
    total = length(modes),
    format = "  Selecting Ranks [:bar] :percent eta: :eta",
    clear = TRUE,
    width = 70
  )
  
  for (n in seq_along(modes)) {
    unfolded <- k_unfold(tnsr, m = n)@data
    svd_res <- svd(unfolded)
    singular_values_sq <- svd_res$d^2
    total_var <- sum(singular_values_sq)
    cumulative_var <- cumsum(singular_values_sq) / total_var
    
    
    r <- which(cumulative_var >= target_var[n])[1]
    ranks[n] <- r
    pb$tick()
  }
  
  return(ranks)
}

# ------------------------------------------------------------------------------
# @brief Applies Tucker decomposition to a 3D output array.
#
# This function performs a Tucker decomposition on a 3D tensor (e.g., space × time × simulation),
# reconstructs the tensor using the reduced ranks, calculates explained variance,
# and generates a set of outer-product basis matrices from the factor matrices.
#
# @param output    A 3D numeric array (dimensions: space × time × simulation).
# @param ranks     A numeric vector of length 3 specifying Tucker ranks for each mode.
# @param threshold A numeric value (0–100) specifying the minimum required explained variance.
#
# @return A named list containing:
#         - ranks: Vector of Tucker ranks used.
#         - core_tensor: The core tensor (reduced representation).
#         - factor_matrices: List containing U1, U2, U3 (mode-specific factor matrices).
#         - output_array_reduced: Reconstructed array from Tucker components.
#         - explained_var_tucker: Variance explained (%) by Tucker decomposition.
#         - explained_var_pca: Variance explained (%) by PCA on mode-3 unfolded matrix.
#         - n_bases: Number of outer-product basis matrices generated.
#         - basis_matrices: List of basis matrices (outer products of U1 and U2 columns).
#         - B_mat: Matrix with vectorized basis matrices as columns.
# ------------------------------------------------------------------------------
apply_tucker_decomp <- function(output, ranks, threshold, method, calculate_bases) {
  # Convert array to tensor for tucker()
  output_tens <- as.tensor(output)
  
  n_modes <- length(dim(output))
  
  # Check that ranks length matches tensor order
  if (length(ranks) != n_modes) {
    stop("Length of `ranks` must match the number of tensor modes.")
  }
  
  # Perform tucker decomposition
  if (method == "hosvd") {
    tucker_decomp <- hosvd(output_tens, ranks) |> suppressMessages()
  } else if (method == "hooi") {
    tucker_decomp <- tucker(output_tens, ranks) |> suppressMessages() 
  }
  
  G   <- tucker_decomp$Z             # core tensor
  U_list <- tucker_decomp$U
  
  # Reconstruct the tensor from Tucker components
  output_tens_reduced     <- ttl(G, U_list, 1:1:n_modes)
  output_array_reduced    <- output_tens_reduced@data
  
  # Calculate the explained variance
  norm_original <- fnorm(output_tens)
  norm_approx   <- fnorm(output_tens_reduced)
  norm_diff     <- fnorm(output_tens - output_tens_reduced)
  
  explained_var_tucker <- 100 * (1 - norm_diff^2 / norm_original^2)
  
  # Warn if explained variance is below threshold
  if (explained_var_tucker < threshold) {
    cat(sprintf(
      "\t\t⚠️ Explained variance by Tucker decomposition is %.2f%%, below the threshold of %.2f%%.\n",
      explained_var_tucker, threshold
    ))
  } else {
    cat("\t\t✅ Tucker met explained variance threshold.\n")
  }
  
  # PCA explained variance
  explained_var_pca <- NA
  if (n_modes >= 2) {
    last_mode <- n_modes
    output_mat <- matrix(output, ncol = dim(output)[last_mode])
    svd_result <- svd(output_mat)
    sdev <- svd_result$d
    explained_vars_pca <- (sdev^2) / sum(sdev^2)
    explained_var_pca <- sum(explained_vars_pca[1:ranks[last_mode]]) * 100
  }
  
  # Basis construction over first (n_modes - 1) modes
  basis_ranks <- ranks[1:(n_modes - 1)]
  bases <- list()
  total_iters <- prod(basis_ranks)
  
  if (!calculate_bases) {
    return(list(
      ranks = ranks,
      core_tensor = G@data,
      factor_matrices = U_list,
      output_array_reduced = output_array_reduced,
      explained_var_tucker = explained_var_tucker,
      explained_var_pca = explained_var_pca,
      bases = NULL
    ))
  }
  
  pb <- progress_bar$new(
    format = "  Creating basis matrices [:bar] :percent ETA: :eta",
    total = total_iters,
    clear = TRUE,
    width = 70
  )
  
  index_grid <- as.matrix(expand.grid(lapply(basis_ranks, seq_len)))
  
  for (i in 1:nrow(index_grid)) {
    idx <- index_grid[i, ]
    # Start with first mode's vector
    tensor_basis <- U_list[[1]][, idx[1]]
    
    for (m in 2:(n_modes - 1)) {
      vec_m <- U_list[[m]][, idx[m]]
      tensor_basis <- outer(tensor_basis, vec_m)
    }
    
    # Flatten to matrix or keep shape depending on your needs
    key <- paste0("B_", paste(idx, collapse = "_"))
    bases[[key]] <- tensor_basis
    pb$tick()
  }
  
  gc()
  
  # Return results as a named list
  return(list(
    ranks = ranks,
    core_tensor = G@data,
    factor_matrices = U_list,
    output_array_reduced = output_array_reduced,
    explained_var_tucker = explained_var_tucker,
    explained_var_pca = explained_var_pca,
    n_bases = length(bases),
    bases = bases
  ))
}


calculate_basis <- function(U_list, idx) {
  Reduce(function(x, y) kronecker(y, x),  # Reverse order here
         Map(function(U, i) U[, i], U_list[-length(U_list)], idx))
}



# ------------------------------------------------------------------------------
# @brief Computes the Gaussian (RBF) covariance matrix between two design matrices.
#
# Given two input design matrices, this function computes a scaled radial basis
# function (RBF) kernel (also known as the squared exponential kernel),
# where each input feature is scaled by a per-dimension range parameter.
#
# @param design_mat_1 A numeric matrix (n1 × d), where each row is a d-dimensional point.
# @param design_mat_2 A numeric matrix (n2 × d), to be compared with design_mat_1.
# @param lambda_w     A positive scalar that scales the overall kernel (inverse variance).
# @param rho_vec      A numeric vector of length d containing positive range parameters
#                     for each input dimension (used to scale the columns).
#
# @return A numeric matrix of dimensions (n1 × n2), representing the scaled
#         Gaussian covariance between the rows of the input matrices.
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# @brief Constructs a block-diagonal covariance matrix using group-specific kernels.
#
# This function builds a large covariance matrix by applying the Gaussian kernel
# (via R_mat) to each group in the third mode. Each group uses its own set of
# kernel hyperparameters (lambda and rho).
#
# @param design_mat_1 A numeric matrix (n1 × d), representing the first design set.
# @param design_mat_2 A numeric matrix (n2 × d), representing the second design set.
# @param lambda_vec   A numeric vector of length r3, where each element is a
#                     lambda (inverse variance) for the corresponding group.
# @param rho_mat      A numeric matrix of shape (r3 × d), where each row contains
#                     the per-dimension range parameters for one group.
# @param r3           (Optional) Integer specifying the number of groups (defaults
#                     to nrow(rho_mat) if NULL).
#
# @return A numeric matrix of dimensions (n1 × r3) × (n2 × r3), where each r3 × r3
#         block is a group-specific covariance matrix between the input design points.
# ------------------------------------------------------------------------------
Sigma_general <- function(design_mat_1, design_mat_2, lambda_vec, rho_mat, r3 = NULL) {
  # Set r3 to number of rows in rho_mat if not provided
  if (is.null(r3)) {
    r3 <- nrow(rho_mat)
  }
  
  d1 <- nrow(design_mat_1)  # Number of input points in design_mat_1
  d2 <- nrow(design_mat_2)  # Number of input points in design_mat_2
  
  # Pre-allocate output matrix for all blocks
  out <- matrix(0, d1 * r3, d2 * r3)
  
  # For each group (i), compute the covariance block using group-specific parameters
  for (i in 1:r3) {
    lambda_i <- lambda_vec[i]     # Inverse variance for group i
    rho_i    <- rho_mat[i, ]      # Range parameters for group i
    
    # Determine row and column block indices for this group
    row_idx_1 <- ((i - 1) * d1 + 1):(i * d1)
    row_idx_2 <- ((i - 1) * d2 + 1):(i * d2)
    
    # Fill block (i, i) with group-specific kernel values
    out[row_idx_1, row_idx_2] <- R_mat(design_mat_1, design_mat_2, lambda_i, rho_i)
  }
  
  return(out)
}

# ------------------------------------------------------------------------------
# @brief Constructs a sparse block-diagonal covariance matrix using group-specific kernels.
#
# This function creates a sparse block-diagonal matrix where each block corresponds to
# a group-specific Gaussian (RBF) kernel computed from two design matrices. It is a
# memory-efficient alternative to the dense version in `Sigma_general()`.
#
# @param design_mat_1 A numeric matrix (n1 × d), representing the first set of design points.
# @param design_mat_2 A numeric matrix (n2 × d), representing the second set of design points.
# @param lambda_vec   A numeric vector of length r3, with each value controlling the scale
#                     (inverse variance) of the corresponding group's kernel.
# @param rho_mat      A numeric matrix (r3 × d), where each row contains per-dimension
#                     range parameters for each group.
# @param r3           (Optional) Integer. Number of groups. Defaults to nrow(rho_mat).
#
# @return A sparse block-diagonal matrix of class "dgCMatrix" (from Matrix package),
#         with dimensions (n1 × r3) × (n2 × r3), where each block is a kernel matrix
#         computed with group-specific parameters.
# ------------------------------------------------------------------------------
Sigma_general_sparse <- function(design_mat_1, design_mat_2, lambda_vec, rho_mat, r3 = NULL) {
  # Default to number of groups equal to number of rows in rho_mat
  if (is.null(r3)) r3 <- nrow(rho_mat)
  
  # Build a list of group-specific kernel matrices (dense blocks)
  block_list <- lapply(seq_len(r3), function(i) {
    R_mat(design_mat_1, design_mat_2, lambda_vec[i], rho_mat[i, ])
  })
  
  # Combine the blocks into a sparse block-diagonal matrix
  Matrix::bdiag(block_list)
}

# ------------------------------------------------------------------------------
# @brief Constructs a block-diagonal covariance matrix from precomputed squared distances.
#
# This function efficiently builds a covariance matrix by applying Gaussian kernels
# to a list of precomputed squared Euclidean distance matrices, one per group.
# Each group's kernel is scaled by the inverse of a corresponding lambda parameter.
#
# @param dists_sq_list A list of numeric matrices, each containing squared pairwise
#                      distances between design points for one group.
# @param lambda_vec    A numeric vector of length r3, where each element scales
#                      the corresponding group's kernel (inverse variance).
#
# @return A numeric matrix of size (n1 × r3) × (n2 × r3) representing a block-diagonal
#         covariance matrix constructed from the group-wise kernels.
# ------------------------------------------------------------------------------
Sigma_general_fast <- function(dists_sq_list, lambda_vec) {
  r3 <- length(lambda_vec)                    # Number of groups
  d1 <- nrow(dists_sq_list[[1]])              # Number of points in design_mat_1 per group
  d2 <- ncol(dists_sq_list[[1]])              # Number of points in design_mat_2 per group
  
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

# ------------------------------------------------------------------------------
# @brief Precomputes and caches squared Euclidean distance matrices scaled by range parameters.
#
# This function computes a list of squared distance matrices between two design matrices,
# where each distance matrix corresponds to a different set of per-dimension scaling
# parameters (rho_i) from the rows of rho_mat.
#
# @param design_1 A numeric matrix (n1 × d), representing the first set of input points.
# @param design_2 A numeric matrix (n2 × d), representing the second set of input points.
# @param rho_mat  A numeric matrix (r3 × d), where each row contains scaling parameters
#                 (range parameters) for a particular component/group.
#
# @return A list of length r3, where each element is an (n1 × n2) matrix containing the
#         squared Euclidean distances between scaled points in design_1 and design_2
#         for the corresponding component.
# ------------------------------------------------------------------------------
initialize_distance_cache_general <- function(design_1, design_2, rho_mat) {
  r3 <- nrow(rho_mat)                # Number of components/groups
  dists_sq_list <- vector("list", r3)  # Preallocate list for squared distance matrices
  
  for (i in 1:r3) {
    rho_i <- rho_mat[i, ]            # Range parameters for component i
    
    # Scale columns of design matrices by 1 / rho_i (element-wise division)
    scaled_1 <- sweep(design_1, 2, rho_i, "/")  # (n1 × d)
    scaled_2 <- sweep(design_2, 2, rho_i, "/")  # (n2 × d)
    
    # Compute pairwise squared Euclidean distances between scaled points
    dists_sq_list[[i]] <- fields::rdist(scaled_1, scaled_2)^2
  }
  
  return(dists_sq_list)
}

# ------------------------------------------------------------------------------
# @brief Updates the squared distance matrix cache for a single component.
#
# Given updated range parameters (rho) for a specific component, this function
# recalculates the squared Euclidean distances between scaled design points,
# updating only the specified element of the cached list.
#
# @param design_1     A numeric matrix (n1 × d), first set of design points.
# @param design_2     A numeric matrix (n2 × d), second set of design points.
# @param rho_mat      A numeric matrix (r3 × d), each row contains range parameters
#                     for each component.
# @param dists_sq_list A list of length r3, the cached squared distance matrices.
# @param i            Integer index specifying which element of the cache to update.
#
# @return The updated list of squared distance matrices with the i-th element recalculated.
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# @brief Computes the density of the half-normal distribution at given points.
#
# The half-normal distribution is the absolute value of a normal distribution with
# mean zero and standard deviation sigma. The density is zero for negative inputs.
#
# @param x      Numeric vector or value at which to evaluate the density.
# @param sigma  Positive numeric value representing the standard deviation of the
#               underlying normal distribution (default = 1).
# @param log    Logical; if TRUE, returns the log-density instead of the density.
#
# @return A numeric vector of the same length as x containing the (log-)density
#         of the half-normal distribution evaluated at x.
# ------------------------------------------------------------------------------
dhalfnorm <- function(x, sigma = 1, log = FALSE) {
  # Density is zero for negative values (half-normal is zero for x < 0)
  if (any(x < 0)) return(0)
  
  # Compute log density using normal density and adding log(2)
  log_density <- log(2) + dnorm(x, mean = 0, sd = sigma, log = TRUE)
  
  # Return log density or density according to user request
  if (log) return(log_density)
  return(exp(log_density))
}

# ------------------------------------------------------------------------------
# @brief Computes the density of the Inverse Gamma distribution at given points.
#
# The Inverse Gamma distribution with shape parameter alpha and scale parameter beta
# is defined only for positive x. The density is zero for non-positive values.
#
# @param x      Numeric vector or value at which to evaluate the density.
# @param alpha  Positive numeric shape parameter of the Inverse Gamma distribution.
# @param beta   Positive numeric scale parameter of the Inverse Gamma distribution.
#
# @return A numeric vector of densities evaluated at x.
# ------------------------------------------------------------------------------
dinvgamma <- function(x, alpha, beta) {
  # Density formula:
  # f(x) = (beta^alpha / Gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x),  x > 0
  ifelse(x > 0,
         (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x),
         0)
}

# ------------------------------------------------------------------------------
# @brief Computes the log-density of a multivariate normal distribution using Cholesky decomposition.
#
# This function efficiently computes the log-density of a multivariate normal
# distribution by avoiding direct matrix inversion. It uses the Cholesky
# decomposition of the covariance matrix for numerical stability and speed.
#
# @param x      A numeric vector representing the observation (length d).
# @param mu     A numeric vector of the same length as x, the mean vector.
# @param Sigma  A positive-definite covariance matrix (d × d).
#
# @return A scalar numeric value: the log-density of the multivariate normal at x.
# ------------------------------------------------------------------------------
log_mvnorm_cholesky <- function(x, mu, Sigma) {
  # Cholesky decomposition: Sigma = L * L^T
  L <- base::chol(Sigma)
  
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

# ------------------------------------------------------------------------------
# @brief Computes the log-density of a multivariate normal distribution using a sparse Cholesky decomposition.
#
# This function evaluates the log-density of a multivariate normal distribution
# assuming a sparse symmetric positive-definite covariance matrix. It leverages
# the sparse Cholesky factorization from the Matrix package for computational
# efficiency with large, sparse matrices.
#
# @param x      Numeric vector of length n, the observation point.
# @param mu     Numeric vector of length n, the mean of the multivariate normal distribution.
# @param Sigma  A sparse symmetric positive-definite covariance matrix of class "dsCMatrix".
#
# @return A scalar numeric value representing the log-density at point x.
# ------------------------------------------------------------------------------
log_mvnorm_cholesky_sparse <- function(x, mu, Sigma) {
  # Ensure Sigma is a sparse symmetric matrix in compressed column format
  stopifnot(class(Sigma)[1] == "dsCMatrix")
  
  # Perform sparse Cholesky decomposition: Sigma = L %*% t(L)
  # Returns a CHMfactor object used for solving systems
  L <- Matrix::Cholesky(Sigma, LDL = FALSE)
  
  # Solve the system L * y = (x - mu), where L is the lower Cholesky factor
  # system = "L" specifies we solve the lower-triangular system
  y <- Matrix::solve(L, x - mu, system = "L")
  
  # Compute log-determinant from Sigma (Matrix handles this safely for sparse matrices)
  log_det_Sigma <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  
  # Compute quadratic form (x - mu)^T * Sigma^{-1} * (x - mu) = y^T * y
  quadratic_form <- sum(y^2)
  
  n <- length(x)  # Dimensionality of the distribution
  
  # Final log-density expression
  log_density <- -0.5 * (n * log(2 * pi) + log_det_Sigma + quadratic_form)
  
  return(log_density)
}

# ------------------------------------------------------------------------------
# @brief Computes the log-density of a multivariate normal by summing over blockwise components.
#
# This function assumes a block-diagonal covariance matrix and evaluates the
# total log-density by computing the log-density of each block separately using
# Cholesky decomposition. This is efficient and numerically stable when the
# covariance structure is block-diagonal or approximated as such.
#
# @param x          Numeric vector of length n, the observation point.
# @param mu         Numeric vector of length n, the mean vector.
# @param Sigma      Full covariance matrix (n × n), assumed to be block-diagonal
#                   or treated as such. Can be dense or sparse.
# @param block_size Integer size of each block (must evenly divide n).
#
# @return A scalar numeric value representing the total log-density of x under
#         the block-structured multivariate normal.
# ------------------------------------------------------------------------------
log_mvnorm_blockwise <- function(x, mu, Sigma, block_size) {
  n <- length(x)
  
  # Input checks
  stopifnot(length(mu) == n, nrow(Sigma) == n, ncol(Sigma) == n)
  stopifnot(n %% block_size == 0)  # Ensure x can be evenly split into blocks
  
  num_blocks <- n / block_size  # Total number of blocks
  
  # Create index list: each element contains indices of a block
  idx_list <- split(seq_len(n), rep(seq_len(num_blocks), each = block_size))
  
  total_logdens <- 0  # Accumulator for total log-density
  
  for (idx in idx_list) {
    # Extract block components of x, mu, and Sigma
    x_block <- x[idx]
    mu_block <- mu[idx]
    Sigma_block <- Sigma[idx, idx]  # Extract block covariance matrix
    
    # Compute and accumulate log-density for the block
    total_logdens <- total_logdens + log_mvnorm_cholesky(x_block, mu_block, Sigma_block)
  }
  
  return(total_logdens)
}

# ------------------------------------------------------------------------------
# @brief Computes the matrix product C = B %*% A using a custom block-wise formulation.
#
# This function builds a large matrix by computing the product between a list of vectors
# and a coefficient matrix C, distributing results across blocks. It assumes that
# each vector in `vectors` is repeated across `n` blocks and scaled by entries in C.
# Progress is tracked with a progress bar.
#
# @param vectors A list of k numeric vectors, each of length m. These are the columns of B.
# @param C       A numeric matrix of dimensions (k × q). Coefficients that scale vectors.
# @param n       Integer. The number of blocks (e.g., spatial or temporal replication).
#
# @return A numeric matrix of size (n * m) × (n * q), where each (i, l)-th column is formed
#         by scaling the i-th vector in `vectors` by C[i, l], and distributing it across blocks.
# ------------------------------------------------------------------------------
# matrix_calc <- function(vectors, C, n) {
#   k <- length(vectors)         # Number of vectors (length of list)
#   m <- length(vectors[[1]])    # Length of each vector
#   q <- ncol(C)                 # Number of output columns per block
#   
#   total_iters <- k * q
#   
#   # Initialize a progress bar to monitor the loop
#   pb <- progress_bar$new(
#     total = total_iters,
#     format = "  Matrix Calculation [:bar] :percent eta: :eta Rate: :tick_rate",
#     clear = TRUE,
#     width = 70
#   )
#   
#   # Initialize the result matrix with appropriate size
#   result <- matrix(0, nrow = n * m, ncol = n * q)
#   
#   # Outer loop over vectors (rows of C)
#   for (i in seq_len(k)) {
#     v <- as.numeric(vectors[[i]])  # Convert i-th vector to numeric
#     for (l in seq_len(q)) {
#       scalar <- C[i, l]            # Get coefficient for this vector-column combination
#       
#       col_start <- (l - 1) * n + 1  # Start index in result columns for block l
#       
#       # Loop over n blocks
#       for (j in seq_len(n)) {
#         row_block_idx <- ((j - 1) * m + 1):(j * m)   # Row indices for block j
#         col_idx <- col_start + j - 1                # Column index for (j,l)-th output
#         
#         # Accumulate scaled vector into correct position
#         result[row_block_idx, col_idx] <- result[row_block_idx, col_idx] + scalar * v
#       }
#       
#       pb$tick()  # Update progress bar after each (i, l) iteration
#     }
#   }
#   
#   return(result)
# }
matrix_calc <- function(vectors, A, x = NULL, d) {
  r <- length(vectors)         # Number of vectors
  n <- length(vectors[[1]])    # Length of each vector (s x t)
  r3 <- ncol(A)                # Number of output columns per block
  
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
        dot_product <- sum(results[[i]] * x[idx])  # scalar
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

# parallel_weighted_matrix_sum <- function(X_list, weights, n_cores) {
#   stopifnot(length(X_list) == length(weights))
#   
#   chunks <- split(seq_along(X_list), cut(seq_along(X_list), n_cores, labels = FALSE))
#   
#   partial_sums <- pbmclapply(chunks, function(idx) {
#     Reduce(`+`, Map(function(i) X_list[[i]] * weights[i], idx))
#   }, mc.cores = n_cores)
#   
#   Reduce(`+`, partial_sums)
# }

matrix_calc_mf <- function(matrices_interp, tuck_discrep, 
                           A_interp, A_discrep, 
                           x = NULL, 
                           d) {
  r_lf <- length(matrices_interp)         # Number of vectors
  r_d_lf <- ncol(A_interp)                # Number of output columns per block
  
  r_discrep <- prod(tuck_discrep$ranks[-length(tuck_discrep$ranks)])
  r_d_discrep <- ncol(A_discrep)                # Number of output columns per block
  
  n <- prod(dim(matrices_interp[[1]]))    # Length of each vector (s_hf x t_hf)
  
  
  results <- vector("list", r_d_lf + r_d_discrep)
  for (j in 1:r_d_lf) {
    weights <- as.vector(A_interp[, j])          # vector of length k
    out <- array(0, dim = dim(matrices_interp[[1]]))
    for (k in 1:length(matrices_interp)) {
      out <- out + weights[k] * matrices_interp[[k]]
      if (k %% 100 == 0) gc()
    }
    cat("\r")
    results[[j]] <- as.vector(out)
  }
  
  U_list <- tuck_discrep$factor_matrices
  
  # Define your original variables
  # r_y_discrep, r_m_discrep, r_s_discrep, A_discrep, U_list, r_d_lf
  
  # 1. Generate all idx_combinations once (this matrix itself should be small enough)
  idx_combinations <- expand.grid(s = 1:r_s_discrep,
                                  m = 1:r_m_discrep,
                                  y = 1:r_y_discrep)
  total_combinations <- nrow(idx_combinations)
  basis_len <- prod(sapply(U_list[-length(U_list)], nrow))
  
  # 2. Determine Block Size (adjust as needed based on RAM and performance)
  # Example: If each basis matrix is 30MB, and you want to load 10GB (10,000 MB) at a time:
  # block_size_MB <- 10000
  # matrices_per_block <- floor(block_size_MB / 30) # ~333 matrices per block
  matrices_per_block <- 300 # Start with a reasonable guess, then fine-tune
  
  num_blocks <- ceiling(total_combinations / matrices_per_block)
  
  # Initialize results. Instead of list of vectors, maybe a pre-allocated matrix for efficiency
  # If results[[r_d_lf + j]] should be a vector of length 'basis_len'
  results_discrep_matrix <- matrix(0, nrow = basis_len, ncol = r_d_discrep)
  
  
  # 3. Loop through blocks
  # (Optional) Progress bar for blocks
  pb_block <- progress_bar$new(
    format = "Processing blocks [:bar] :current/:total (:percent) ETA: :eta",
    total = num_blocks, width = 70
  )
  
  for (b in 1:num_blocks) {
    start_idx <- (b - 1) * matrices_per_block + 1
    end_idx <- min(b * matrices_per_block, total_combinations)
    current_block_indices <- start_idx:end_idx
    
    # Get the A_discrep weights for the current block (subset the counter)
    # This corresponds to your 'weights[counter]' in the original code
    # A_discrep has 'total_combinations' rows (your original 'counter' values)
    # and 'r_d_discrep' columns (your original 'j' values)
    weights_for_block <- A_discrep[current_block_indices, , drop = FALSE]
    
    # Generate all basis matrices for the current block
    # This is the 'vapply' step from earlier, but on a subset of combinations
    current_basis_matrices_flat <- vapply(
      current_block_indices,
      FUN = function(i) {
        idx <- as.numeric(idx_combinations[i, ])
        as.vector(calculate_basis(U_list, idx))
      },
      FUN.VALUE = numeric(basis_len)
    )
    # current_basis_matrices_flat is (basis_len x length(current_block_indices))
    
    # Perform the weighted sum for the current block using matrix multiplication
    # (basis_len x N_block) %*% (N_block x r_d_discrep) -> (basis_len x r_d_discrep)
    block_sum_results <- current_basis_matrices_flat %*% weights_for_block
    
    # Accumulate results for each 'j' (column)
    results_discrep_matrix <- results_discrep_matrix + block_sum_results
    
    rm(current_basis_matrices_flat, weights_for_block, block_sum_results) # Free memory
    gc(verbose = FALSE) # Force garbage collection after each block to manage memory
    
    pb_block$tick()
  }
  
  # Assign accumulated results back to your 'results' list
  for (j in 1:r_d_discrep) {
    results[[r_d_lf + j]] <- results_discrep_matrix[, j]
  }
  rm(results_discrep_matrix) # Clean up if it was a large temporary object
  gc(verbose = FALSE)
  
  
  # pb <- progress_bar$new(
  #   format = "  Matrix Calc [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
  #   total = r_d_discrep * r_y_discrep * r_m_discrep * r_s_discrep,
  #   width = 70
  # )
  # 
  # for (j in 1:r_d_discrep) {
  #   weights <- as.vector(A_discrep[, j])
  #   out <- array(0, dim = lapply(U_list, nrow)[-length(U_list)]) |> as.vector()
  #   print(dim(out))
  #   counter <- 1
  # 
  #   for (y in 1:r_y_discrep) {
  #     for (m in 1:r_m_discrep) {
  #       for (s in 1:r_s_discrep) {
  #         basis_matrix <- calculate_basis(U_list,
  #                                         idx = c(s, m, y))
  #         out <- out + weights[counter] * basis_matrix
  #         counter <- counter + 1
  # 
  #         pb$tick()
  #       }
  #       gc()
  #     }
  #   }
  # 
  #   results[[r_d_lf + j]] <- as.vector(out)
  # 
  #   gc()
  # }
  
  
  
  
  
  # library(parallel)
  # library(progress)
  # 
  # # ---- Setup ----
  # U_list <- tuck_discrep$factor_matrices
  # n_cores <- 4  # Adjust this for your system
  # r_j <- r_d_discrep
  # 
  # # Split j's into chunks (one per core)
  # chunks <- split(1:r_j, cut(1:r_j, n_cores, labels = FALSE))
  # 
  # # ---- Parallel worker function ----
  # process_chunk <- function(j_list, chunk_id, A_discrep, U_list) {
  #   library(progress)
  #   
  #   pb <- progress_bar$new(
  #     format = paste0("Chunk ", chunk_id, " [:bar] :current/:total (:percent) ETA: :eta RATE: :rate"),
  #     total = length(j_list) * r_y_discrep * r_m_discrep * r_s_discrep,
  #     clear = FALSE, width = 70
  #   )
  #   
  #   results_local <- list()
  #   
  #   for (j in j_list) {
  #     weights <- as.vector(A_discrep[, j])
  #     out_dim <- sapply(U_list[-length(U_list)], nrow)
  #     out <- numeric(prod(out_dim))
  #     
  #     counter <- 1
  #     for (y in 1:r_y_discrep) {
  #       for (m in 1:r_m_discrep) {
  #         for (s in 1:r_s_discrep) {
  #           idx <- c(s, m, y)
  #           basis_matrix <- Reduce(kronecker, Map(function(U, i) U[, i], U_list[-length(U_list)], idx))
  #           out <- out + weights[counter] * basis_matrix
  #           counter <- counter + 1
  #           if (counter %% 100 == 0) pb$tick(100)          }
  #         gc(FALSE)
  #       }
  #     }
  #     
  #     results_local[[as.character(j)]] <- out
  #     gc(FALSE)
  #   }
  #   
  #   return(results_local)
  # }
  # 
  # # ---- Run parallel jobs ----
  # jobs <- vector("list", n_cores)
  # 
  # for (i in seq_len(n_cores)) {
  #   jobs[[i]] <- mcparallel(process_chunk(
  #     chunks[[i]], i, A_discrep, U_list
  #   ))
  # }
  # results_chunks <- mccollect(jobs)
  # 
  # # ---- Merge results ----
  # for (chunk in results_chunks) {
  #   for (j_name in names(chunk)) {
  #     j <- as.integer(j_name)
  #     results[[r_d_lf + j]] <- chunk[[j_name]]
  #   }
  # }
  
  # Progress bar for coef_matrix calculation (r3 x r3 iterations)
  pb3 <- progress_bar$new(
    format = "  Computing coef_matrix [:bar] :percent eta: :eta",
    total = (r_d_lf + r_d_discrep) * (r_d_lf + r_d_discrep), clear = TRUE, width = 70
  )
  
  coef_matrix <- matrix(0, nrow = (r_d_lf + r_d_discrep), ncol = (r_d_lf + r_d_discrep))
  # Compute dot product for each pair
  for (i in 1:(r_d_lf + r_d_discrep)) {
    for (j in 1:(r_d_lf + r_d_discrep)) {
      coef_matrix[i, j] <- sum(results[[i]] * results[[j]])
      pb3$tick()
    }
  }
  
  solve_Kt_K <- kronecker(solve(coef_matrix), Diagonal(d))
  
  if (!is.null(x)) {
    total_iters <- (r_d_lf + r_d_discrep) * d
    pb4 <- progress_bar$new(
      format = "  Computing Kt_x [:bar] :percent eta: :eta",
      total = total_iters, clear = TRUE, width = 70
    )
    results_2 <- vector("numeric", total_iters)
    counter <- 1
    for (i in 1:(r_d_lf + r_d_discrep)) {
      for (j in 1:d) {
        idx <- ((j - 1) * n + 1):(j * n)
        dot_product <- sum(results[[i]] * x[idx])  # scalar
        results_2[counter] <- dot_product
        counter <- counter + 1
        pb4$tick()
      }
    }
    
    Kt_x <- unlist(results_2)
    
    return(list(
      solve_Kt_K = solve_Kt_K,
      Kt_x       = Kt_x
    ))
  }
}

# ------------------------------------------------------------------------------
# @brief Interpolates values at new spatial locations using k-nearest neighbors regression.
#
# This function uses k-NN regression to interpolate output values from known coordinates
# (`coords_from`) to new locations (`coords_to`) based on spatial proximity. It leverages
# the `FNN::knn.reg` function for fast nearest neighbor lookup and prediction.
#
# @param coords_from A numeric matrix or data frame (n × d), the known spatial locations.
# @param outputs     A numeric vector or matrix of length n (or n × p), the output values
#                    associated with `coords_from`.
# @param coords_to   A numeric matrix or data frame (m × d), the target spatial locations
#                    where interpolation is desired.
# @param k           Integer. The number of nearest neighbors to use for interpolation.
#
# @return A numeric vector or matrix of interpolated values at the new coordinates (`coords_to`).
#         Dimensions match those of `outputs` if matrix-valued.
# ------------------------------------------------------------------------------
interpolate_locations <- function(coords_from, outputs, coords_to, k) {
  knn_result <- knn.reg(
    train = coords_from,  # Known spatial coordinates
    test  = coords_to,    # Target locations to interpolate
    y     = outputs,      # Known outputs
    k     = k             # Number of neighbors
  )
  return(knn_result$pred)  # Return predicted interpolated values
}

# ------------------------------------------------------------------------------
# @brief Performs Gram-Schmidt orthogonalization on the columns of a matrix.
#
# This function orthogonalizes the columns of a matrix A using the classical
# Gram-Schmidt process. The result is a matrix Q whose columns are mutually
# orthogonal vectors spanning the same column space as A.
#
# Note: This implementation does NOT normalize the resulting vectors.
#
# @param A  A numeric matrix of size (m × n), where each column represents a vector.
#
# @return A numeric matrix Q of size (m × n) containing orthogonalized (but not normalized)
#         columns corresponding to the input matrix A.
# ------------------------------------------------------------------------------
gram_schmidt <- function(A) {
  n <- ncol(A)  # Number of vectors to orthogonalize (columns)
  m <- nrow(A)  # Dimensionality of each vector
  
  Q <- matrix(0, nrow = m, ncol = n)  # Initialize orthogonal matrix
  
  for (i in 1:n) {
    v <- A[, i]  # Start with the i-th column of A
    
    if (i != 1) {
      # Subtract projections of v onto all previous orthogonal vectors
      for (j in 1:(i - 1)) {
        proj <- sum(v * Q[, j]) / sum(Q[, j]^2) * Q[, j]  # Projection of v onto Q[,j]
        v <- v - proj  # Remove component along Q[,j]
      }
    }
    
    Q[, i] <- v  # Store the orthogonalized vector
  }
  
  return(Q)
}

# ------------------------------------------------------------------------------
# @brief Interpolates a list of basis matrices (from Tucker decomposition) to new coordinates.
#
# This function takes the output of a Tucker decomposition (typically from a low-fidelity model),
# and interpolates each basis matrix column-wise from low-fidelity locations to high-fidelity
# locations using k-nearest neighbors regression.
#
# It assumes that `coord_lf` (low-fidelity coordinates) and `coord_hf` (high-fidelity coordinates)
# exist in the enclosing environment (or global scope).
#
# @param tucker_output_lf A list containing Tucker decomposition output, including basis matrices.
#                         This must include `basis_matrices`, a list of matrix-valued bases.
# @param k                Integer. The number of nearest neighbors to use for interpolation.
#
# @return A list containing:
#         - basis_lf_interp: A list of interpolated basis matrices.
#         - B_eta_tilde: A single matrix where all interpolated bases are vectorized and
#                        stacked column-wise.
# ------------------------------------------------------------------------------
interpolate_bases <- function(tucker_output_lf, k) {
  bases_lf <- tucker_output_lf$bases  # List of basis matrices
  
  G <- tucker_output_lf$core_tensor
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  
  if (is.null(bases_lf)) {
    message("LF bases must be calculated.")
    return(NULL)
  }
  
  # Number of rows in output = space_dim
  # Number of columns in output = ncol(G_unfold)
  n_out_cols <- ncol(G_unfold)
  n_out_rows <- prod(nrow(coord_hf), dim(bases_lf[[1]])[-1])
  
  # Initialize output matrix with zeros
  B_eta_tile_x_G_unfold <- matrix(0, nrow = n_out_rows, ncol = n_out_cols)
  
  # Keep track of column indices for each basis in B_eta_tilde
  start_col <- 1
  
  pb <- progress_bar$new(
    total = length(bases_lf),
    format = " Calculating B_eta_tilde %*% G_unfold [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
    clear = TRUE,
    width = 100
  )
  
  for (i in seq_along(bases_lf)) {
    # Interpolate this basis (array with dims [space_dim, ...])
    interp_basis <- fast_knn_interpolate_first_dim(
      data_array = bases_lf[[i]],
      coords_old = coord_lf,
      coords_new = coord_hf,
      k = k
    )
    
    vec_interp_basis <- matrix(as.vector(interp_basis), ncol = 1)
    
    
    G_sub <- G_unfold[i, , drop = FALSE]
    
    # Multiply: (space_dim x vec_len) %*% (vec_len x columns_subset_of_G) → space_dim x columns_subset_of_G
    # But G_sub columns correspond to basis vectorization — so this is a matrix multiplication for partial cols
    
    # Here G_unfold columns are vectorized bases * modes, so multiplication is valid
    partial_res <- vec_interp_basis %*% G_sub
    
    # Accumulate into output matrix
    B_eta_tile_x_G_unfold <- B_eta_tile_x_G_unfold + partial_res
    
    rm(interp_basis, vec_interp_basis, G_sub, partial_res)
    
    if (i %% 10 == 0) {
      gc()
    }
    
    pb$tick()
  }
  
  # Interpolate factor matrices first mode
  factor_matrices_interp <- tucker_output_lf$factor_matrices
  factor_matrices_interp[[1]] <- apply(
    factor_matrices_interp[[1]],
    MARGIN = 2,
    interpolate_locations,
    coords_from = coord_lf,
    coords_to = coord_hf,
    k = k
  )
  
  return(list(
    factor_matrices_interp = factor_matrices_interp,
    B_eta_tile_x_G_unfold = B_eta_tile_x_G_unfold
  ))
}

prepare_sf_mcmc <- function(output, tucker_outputs, ranks) {
  # Number of simulations (mode-3 dimension)
  d_train <- dim(output)[length(dim(output))]
  
  # Extract core tensor and vectorized basis matrix
  G <- tucker_outputs$core_tensor  # Core tensor (3D array)
  bases <- tucker_outputs$bases  
  
  # Unfold the core tensor along last mode (simulations) and transpose
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  
  # Flatten the 3D output array into a long vector
  output_vec <- as.vector(output)
  
  r_d <- ranks[length(ranks)]
  
  n_out <- prod(dim(output)[-length(dim(output))])
  
  # if bases already calculated
  if (!is.null(tucker_outputs$bases)) {
    matrix_calc_results <- matrix_calc(lapply(bases, as.vector), 
                                       G_unfold, 
                                       output_vec, 
                                       d_train)
    
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
        
        out <- out + G_unfold[i , k] * b
        
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
        z_d <- as.vector(output[, , , d])
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
  b_eta_prime <- as.numeric(b_eta_prime)  # Ensure scalar output
  
  
  
  return(list(
    solve_C_t_C  = solve_Ct_C,
    gamma_hat    = gamma_hat,
    a_eta_prime  = a_eta_prime,
    b_eta_prime  = b_eta_prime
  ))
}

# ------------------------------------------------------------------------------
# @brief Prepares matrices and posterior parameters for multifidelity MCMC sampling.
#
# This function builds a joint model combining low-fidelity output and discrepancy terms 
# (between high- and low-fidelity models). It constructs the design matrix `K` and computes
# posterior estimates of latent coefficients, as well as updated inverse-gamma parameters.
#
# Global dependencies:
#   - Uses: a_delta, b_delta, s_hf, t_hf, r_d_lf, r_d_discrep
#   - These should ideally be passed in for clarity and testability.
#
# @param output_hf             A 3D numeric array (space × time × simulations) for high-fidelity output.
# @param tucker_outputs_lf     Tucker decomposition result from low-fidelity model.
# @param tucker_outputs_discrep Tucker decomposition result from discrepancy model (HF - LF).
# @param interp_bases_outputs  List containing interpolated low-fidelity basis matrices, including:
#                                - B_eta_tilde: Vectorized interpolated LF bases (used for C_tilde).
#
# @return A named list with:
#         - solve_K_t_K      : Inverse of KᵗK (used in posterior updates).
#         - gamma_hat_hf     : Posterior mean estimate for LF coefficients.
#         - zeta_hat         : Posterior mean estimate for discrepancy coefficients.
#         - a_delta_prime    : Updated shape parameter for inverse-gamma prior.
#         - b_delta_prime    : Updated scale parameter for inverse-gamma prior.
# ------------------------------------------------------------------------------
prepare_mf_mcmc <- function(output_hf, tucker_outputs_lf, tucker_outputs_discrep, interp_bases_outputs, ranks_lf, ranks_discrep) {
  
  # Number of high-fidelity simulations
  d_train_hf <- dim(output_hf)[length(dim(output_hf))]
  
  # Extract core tensors
  G         <- tucker_outputs_lf$core_tensor
  G_discrep <- tucker_outputs_discrep$core_tensor
  
  # Unfold core tensors along last mode (design) and transpose
  G_unfold         <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  G_unfold_discrep <- t(k_unfold(as.tensor(G_discrep), m = length(dim(G_discrep)))@data)
  
  
  # Extract interpolated LF bases and discrepancy bases
  # bases_interp <-  interp_bases_outputs$basis_lf_interp # Interpolated low-fidelity basis
  bases_discrep <- tucker_outputs_discrep$basis_matrices
  
  # if(is.null(bases_interp)) {
  #   message("LF (and interpolated LF) bases must be calculated.")
  #   return(NULL)  # or return(some default value)
  # }
  
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
      total = r_d_lf * nrow(index_grid_lf) + r_d_discrep * nrow(index_grid_discrep) + (r_d_lf + r_d_discrep)* d_train_hf,
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
        
        out <- out + G_unfold[i , k] * b
        
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
        
        out <- out + G_unfold_discrep[i , k] * b
        
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
        z_d <- as.vector(output_hf[, , , d])
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
  zeta_hat     <- u_hat[-(1:(r_d_lf * d_train_hf))]
  
  # Update inverse-gamma hyperparameters for variance (discrepancy noise)
  a_delta_prime <- a_delta + (d_train_hf * (s_hf * t_hf - r_d_lf - r_d_discrep)) / 2
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
    b_delta_prime = b_delta_prime
  ))
}


is_pos_def <- function(mat) {
  res <- tryCatch({
    chol(mat)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(res)
}

# ------------------------------------------------------------------------------
# @brief Runs adaptive Metropolis-Hastings MCMC for single-fidelity model parameters.
#
# This function performs Metropolis-Hastings updates on parameters lambda_eta, 
# lambda_w_vec, and rho_mat using adaptive tuning of proposal variances to target
# a desired acceptance rate. It samples from the posterior of these parameters
# given observed output and design.
#
# Global dependencies:
#   - Requires functions: initialize_distance_cache_general, update_distance_cache_element_general,
#     Sigma_general_fast, log_mvnorm_cholesky, dhalfnorm, dinvgamma
#   - Requires hyperparameters: a_eta_prime, b_eta_prime, sigma_lambda_w, a_rho_w, b_rho_w
#
# @param output               3D numeric array (space × time × simulations) of observed outputs.
# @param design               Numeric matrix of design points (simulations × features).
# @param mcmc_setup_outputs   List containing:
#                               - solve_C_t_C: matrix for precision computations,
#                               - gamma_hat: posterior mean vector,
#                               - a_eta_prime: shape parameter for inverse-gamma prior,
#                               - b_eta_prime: rate parameter for inverse-gamma prior.
# @param n_iter               Integer, total number of MCMC iterations.
# @param burn_in              Integer, number of initial iterations to discard.
# @param r3                   Integer, latent dimension size.
# @param init_lambda_eta      Optional numeric scalar, initial value for lambda_eta.
# @param init_lambda_w_vec    Optional numeric vector, initial values for lambda_w_vec.
# @param init_rho_mat         Optional numeric matrix, initial values for rho_mat.
# @param omega_init           List with initial proposal standard deviations for MH updates:
#                               - lambda_eta: numeric scalar,
#                               - lambda_w: numeric scalar or vector,
#                               - rho: numeric scalar or matrix.
#
# @return A list with:
#           - n_iter: total MCMC iterations,
#           - burn_in: burn-in iterations,
#           - init_params: initial parameter values,
#           - post_chains: posterior chains for lambda_eta, lambda_w_vec, rho_mat, and log likelihood.
# ------------------------------------------------------------------------------
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
    verbose = F
) {
  d_train <- dim(output)[length(dim(output))]
  
  # Extract fixed components
  solve_C_t_C <- mcmc_setup_outputs$solve_C_t_C
  gamma_hat   <- mcmc_setup_outputs$gamma_hat
  a_eta_prime <- mcmc_setup_outputs$a_eta_prime
  b_eta_prime <- mcmc_setup_outputs$b_eta_prime
  
  # Initialize parameters
  lambda_eta   <- init_lambda_eta   %||% (a_eta_prime / b_eta_prime)
  lambda_w_vec <- init_lambda_w_vec %||% rep(1, r3)
  rho_mat      <- init_rho_mat      %||% matrix(2, nrow = r3, ncol = ncol(design))
  
  # Initialize squared-distance cache
  dists_sq_list <- initialize_distance_cache_general(design, design, rho_mat)
  
  # Log-posterior function
  log_posterior <- function(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat) {
    Sigma <- lambda_eta^(-1) * solve_C_t_C + Sigma_general_fast(dists_sq_list, lambda_w_vec)
    
    # Check if Sigma is positive definite
    if (!is_pos_def(Sigma)) {
      print("⚠️ Sigma not positive definite — applying nearPD fix.")
      
      Sigma <- as.matrix(forceSymmetric(
        nearPD(Sigma, ensureSymmetry = TRUE)$mat
      ))
      
      # Add small jitter to diagonal to ensure strict PD
      # Sigma <- Sigma + diag(1e-8, nrow(Sigma))
    }
    
    ll <- log_mvnorm_cholesky(
      as.vector(gamma_hat),
      mu = rep(0, d_train * r3),
      Sigma = Sigma
    )
    prior_lambda_eta <- dgamma(lambda_eta, a_eta_prime, b_eta_prime, log = TRUE)
    prior_lambda_w <- sum(dhalfnorm(lambda_w_vec, sigma_lambda_w, log = TRUE))
    # prior_rho <- sum(log(dinvgamma(rho_mat, a_rho_w, b_rho_w)))
    prior_rho <- sum(dlnorm(rho_mat, a_rho_w, b_rho_w, log = TRUE))
    # prior_rho <- sum(dunif(rho_mat, 0, 100, log = T))
    
    # if (any(rho_mat > 12)) return(-Inf)
    
    return(as.numeric(ll + prior_lambda_eta + prior_lambda_w + prior_rho))
  }
  
  # --- Metropolis updates ---
  update_lambda_eta_MH <- function(lambda_eta, ...) {
    
    # Propose from truncated normal (truncated at 0 from below)
    lambda_eta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_eta, sd = omega_lambda_eta)
    
    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_eta_prop, a = 0, b = Inf, mean = lambda_eta, sd = omega_lambda_eta) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_eta, a = 0, b = Inf, mean = lambda_eta_prop, sd = omega_lambda_eta) |> log()
    
    # Compute log posterior difference
    log_post_current <- log_posterior(lambda_eta, ...)
    log_post_prop    <- log_posterior(lambda_eta_prop, ...)
    
    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = lambda_eta_prop, accept = TRUE))
    } else {
      return(list(value = lambda_eta, accept = FALSE))
    }
  }
  
  update_lambda_w_vec_i_MH <- function(i, lambda_eta, lambda_w_vec, ...) {
    prop <- lambda_w_vec
    
    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega_lambda_w[i])
    prop[i] <- prop_i_new
    
    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega_lambda_w[i]) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_w_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega_lambda_w[i]) |> log()
    
    # Calculate log posterior difference
    log_post_current <- log_posterior(lambda_eta, lambda_w_vec, ...)
    log_post_prop    <- log_posterior(lambda_eta, prop, ...)
    
    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE))
    } else {
      return(list(value = lambda_w_vec, accept = FALSE))
    }  
  }
  
  update_rho_mat_ij_MH <- function(i, j, lambda_eta, lambda_w_vec, rho_mat, dists_sq_list) {
    rho_pro <- rho_mat
    
    # Propose new rho[i,j] from truncated normal > 0
    rho_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_mat[i, j], sd = omega_rho[i, j])
    rho_pro[i, j] <- rho_prop_ij
    
    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_prop_ij, a = 0, b = Inf, mean = rho_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_mat[i, j], a = 0, b = Inf, mean = rho_prop_ij, sd = omega_rho[i, j]) |> log()
    
    # Update the distance cache with the proposed rho
    dists_sq_list_pro <- update_distance_cache_element_general(design, design, rho_pro, dists_sq_list, i)
    
    # Compute log posterior difference
    log_post_current <- log_posterior(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat)
    log_post_prop    <- log_posterior(lambda_eta, lambda_w_vec, dists_sq_list_pro, rho_pro)
    
    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = rho_pro, dists_sq_list = dists_sq_list_pro, accept = TRUE))
    } else {
      return(list(value = rho_mat, dists_sq_list = dists_sq_list, accept = FALSE))
    }  
  }
  
  # Storage
  n_save <- n_iter - burn_in
  lambda_eta_chain   <- numeric(n_save)
  lambda_w_vec_chain <- matrix(NA, nrow = n_save, ncol = r3)
  rho_mat_chain      <- array(NA, dim = c(n_save, r3, ncol(design)))
  log_lik_chain      <- numeric(n_save)
  
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
  
  for (iter in 1:n_iter) {
    
    # lambda_eta
    res <- update_lambda_eta_MH(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat)
    lambda_eta <- res$value
    if (iter <= burn_in) {
      accept_lambda_eta <- accept_lambda_eta + res$accept
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_eta / update_omega_iter
        old_omega <- omega_lambda_eta
        omega_lambda_eta <- omega_lambda_eta * exp((rate - target_accept) / sqrt(iter) * learning_rate)
        accept_lambda_eta <- 0
        if (verbose) cat(paste0("\n\t\tUpdate omega_lambda_eta:\t", round(old_omega, 2), "\t->\t", round(omega_lambda_eta, 4), "\t(Accept = ", round(rate, 2), ")\t(Current = ", round(lambda_eta, 2), ")\n"))
      }
    }
    
    # lambda_w
    for (j in 1:r3) {
      res <- update_lambda_w_vec_i_MH(j, lambda_eta, lambda_w_vec, dists_sq_list, rho_mat)
      lambda_w_vec <- res$value
      if (iter <= burn_in) {
        accept_lambda_w[j] <- accept_lambda_w[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_w[j] / update_omega_iter
          old_omega <- omega_lambda_w[j]
          omega_lambda_w[j] <- omega_lambda_w[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_w[j] <- 0
          if (verbose) cat(paste0("\t\tUpdate omega_lambda_w[", j, "]:\t", round(old_omega, 2), "\t->\t", round(omega_lambda_w[j], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(lambda_w_vec[j], 2), ")\n"))
        }
      }
    }
    
    # rho
    for (i_rho in 1:r3) {
      for (j_rho in 1:ncol(design)) {
        res <- update_rho_mat_ij_MH(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_mat, dists_sq_list)
        rho_mat <- res$value
        dists_sq_list <- res$dists_sq_list
        if (iter <= burn_in) {
          accept_rho[i_rho, j_rho] <- accept_rho[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho[i_rho, j_rho] / update_omega_iter
            old_omega <- omega_rho[i_rho, j_rho]
            omega_rho[i_rho, j_rho] <- omega_rho[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho[i_rho, j_rho] <- 0
            if (verbose) cat(paste0("\t\tUpdate omega_rho[", i_rho, ", ", j_rho, "]:\t\t", round(old_omega, 2), "\t->\t", round(omega_rho[i_rho, j_rho], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(rho_mat[i_rho, j_rho], 2), ")\n"))
          }
        }
      }
    }
    
    # Store samples
    if (iter > burn_in) {
      idx <- iter - burn_in
      lambda_eta_chain[idx]   <- lambda_eta
      lambda_w_vec_chain[idx, ] <- lambda_w_vec
      rho_mat_chain[idx, , ]  <- rho_mat
      log_lik_chain[idx]      <- log_posterior(lambda_eta, lambda_w_vec, dists_sq_list, rho_mat)
    }
    
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

# ------------------------------------------------------------------------------
# @brief Performs multifidelity MCMC sampling for joint posterior estimation.
#
# This function runs MCMC sampling combining low-fidelity and high-fidelity models,
# including discrepancy terms. It updates hyperparameters controlling covariance
# structures via Metropolis-Hastings steps and adapts proposal step sizes during
# burn-in. The function returns posterior chains for precision parameters and
# correlation matrices for both fidelity levels.
#
# Global dependencies:
#   - Uses: r_d_lf, r_d_discrep, design_train_lf, design_train_hf, r_d_lf, r_d_discrep,
#           sigma_lambda_w, sigma_lambda_v, a_rho_w, b_rho_w, a_rho_v, b_rho_v,
#           progress_bar, log_mvnorm_cholesky, Sigma_general, dhalfnorm, dinvgamma
#
# @param mcmc_lf_setup_outputs  Named list with low-fidelity MCMC setup outputs, including
#                               solve_C_t_C, gamma_hat, a_eta_prime, b_eta_prime.
# @param mcmc_mf_setup_outputs  Named list with multifidelity MCMC setup outputs, including
#                               solve_K_t_K, gamma_hat_hf, zeta_hat, a_delta_prime, b_delta_prime.
# @param n_iter                 Integer specifying total number of MCMC iterations.
# @param burn_in                Integer specifying number of burn-in iterations.
#
# @return A named list containing:
#         - n_iter            : Total iterations run.
#         - burn_in           : Burn-in iterations.
#         - init_params       : Initial values of precision and correlation parameters.
#         - post_chains       : Posterior samples for precision parameters (lambda_eta, lambda_delta),
#                               correlation vectors/matrices (lambda_w_vec, rho_w_mat, lambda_v_vec, rho_v_mat),
#                               and log-posterior likelihood values.
# ------------------------------------------------------------------------------
perform_mf_mcmc <- function(
    mcmc_lf_setup_outputs, mcmc_mf_setup_outputs, n_iter, burn_in,
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
    verbose = FALSE
) {
  
  solve_C_t_C <- mcmc_lf_setup_outputs$solve_C_t_C
  solve_K_t_K <- mcmc_mf_setup_outputs$solve_K_t_K
  
  gamma_hat_lf <- mcmc_lf_setup_outputs$gamma_hat
  gamma_hat_hf <- mcmc_mf_setup_outputs$gamma_hat_hf
  zeta_hat    <- mcmc_mf_setup_outputs$zeta_hat
  hat_vals <- c(as.vector(gamma_hat_lf), as.vector(gamma_hat_hf), as.vector(zeta_hat))
  
  a_eta_prime   <- mcmc_lf_setup_outputs$a_eta_prime
  b_eta_prime   <- mcmc_lf_setup_outputs$b_eta_prime
  a_delta_prime <- mcmc_mf_setup_outputs$a_delta_prime
  b_delta_prime <- mcmc_mf_setup_outputs$b_delta_prime
  
  log_posterior_mfm <- function(lambda_eta, lambda_w_vec, rho_w_mat, 
                                lambda_delta, lambda_v_vec, rho_v_mat) {
    # --- Build Sigma.1 (dense) ---
    Sigma_C <- lambda_eta^(-1) * solve_C_t_C
    Sigma_K <- lambda_delta^(-1) * solve_K_t_K
    Sigma.1 <- matrix(0, nrow = nrow(Sigma_C) + nrow(Sigma_K), ncol = ncol(Sigma_C) + ncol(Sigma_K))
    Sigma.1[1:nrow(Sigma_C), 1:ncol(Sigma_C)] <- as.matrix(Sigma_C)
    Sigma.1[(nrow(Sigma_C) + 1):(nrow(Sigma_C) + nrow(Sigma_K)),
            (ncol(Sigma_C) + 1):(ncol(Sigma_C) + ncol(Sigma_K))] <- as.matrix(Sigma_K)
    
    # --- Build Sigma.2 (dense) ---
    Sigma11 <- Sigma_general(design_train_lf, design_train_lf, lambda_w_vec, rho_w_mat, r_d_lf)
    Sigma22 <- Sigma_general(design_train_hf, design_train_hf, lambda_w_vec, rho_w_mat, r_d_lf)
    Sigma33 <- Sigma_general(design_train_hf, design_train_hf, lambda_v_vec, rho_v_mat, r_d_discrep)
    
    Sigma12 <- Sigma_general(design_train_lf, design_train_hf, lambda_w_vec, rho_w_mat, r_d_lf)
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
    
    # Check if Sigma is positive definite
    if (!is_pos_def(Sigma)) {
      print("⚠️ Sigma not positive definite — applying nearPD fix.")
      
      Sigma <- as.matrix(forceSymmetric(
        nearPD(Sigma, ensureSymmetry = TRUE)$mat
      ))
      
      # Add small jitter to diagonal to ensure strict PD
      # Sigma <- Sigma + diag(1e-8, nrow(Sigma))
    }
    
    
    out.1 <- log_mvnorm_cholesky(hat_vals, rep(0, length(hat_vals)), Sigma)
    
    out.2 <- dgamma(lambda_eta, a_eta_prime, b_eta_prime, log = TRUE)
    
    out.3 <- sum(dhalfnorm(lambda_w_vec, sigma_lambda_w, log = TRUE))
    
    out.4 <- sum(dlnorm(rho_w_mat, a_rho_w, b_rho_w, log = TRUE))
    # out.4[any(rho_w_mat > 30)] <- -Inf
    
    out.5 <- dgamma(lambda_delta, a_delta_prime, b_delta_prime, log = TRUE)
    
    out.6 <- sum(dhalfnorm(lambda_v_vec, sigma_lambda_v, log = TRUE))
    
    out.7 <- sum(dlnorm(rho_v_mat, a_rho_v, b_rho_v, log = TRUE))
    # out.7[any(rho_v_mat > 30)] <- -Inf
    
    
    return(as.numeric(out.1 + out.2 + out.3 + out.4 + out.5 + out.6 + out.7))
  }
  
  update_lambda_eta_MH_MFM <- function(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega) {
    
    # Propose from truncated normal (truncated at 0 from below)
    lambda_eta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_eta, sd = omega)
    
    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_eta_prop, a = 0, b = Inf, mean = lambda_eta, sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_eta, a = 0, b = Inf, mean = lambda_eta_prop, sd = omega) |> log()
    
    # Compute log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta_prop, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    
    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = lambda_eta_prop, accept = TRUE))
    } else {
      return(list(value = lambda_eta, accept = FALSE))
    }
  }
  
  update_lambda_w_vec_i_MH_MFM <- function(i, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega) {
    prop <- lambda_w_vec
    
    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega)
    prop[i] <- prop_i_new
    
    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_w_vec[i], sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_w_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega) |> log()
    
    # Calculate log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta, prop, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    
    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE))
    } else {
      return(list(value = lambda_w_vec, accept = FALSE))
    }  
  }
  
  update_rho_w_mat_ij_MH_MFM <- function(i, j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_rho) {
    rho_w_pro <- rho_w_mat
    
    # Propose new rho[i,j] from truncated normal > 0
    rho_w_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_w_mat[i, j], sd = omega_rho[i, j])
    rho_w_pro[i, j] <- rho_w_prop_ij
    
    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_w_prop_ij, a = 0, b = Inf, mean = rho_w_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_w_mat[i, j], a = 0, b = Inf, mean = rho_w_prop_ij, sd = omega_rho[i, j]) |> log()
    
    # Compute log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_pro, lambda_delta, lambda_v_vec, rho_v_mat)
    
    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = rho_w_pro, accept = TRUE))
    } else {
      return(list(value = rho_w_mat, accept = FALSE))
    }  
  }
  
  update_lambda_delta_MH_MFM <- function(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega) {
    
    # Propose from truncated normal (truncated at 0 from below)
    lambda_delta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_delta, sd = omega)
    
    # Compute log proposal densities
    log_q_current_to_prop <- dtruncnorm(lambda_delta_prop, a = 0, b = Inf, mean = lambda_delta, sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_delta, a = 0, b = Inf, mean = lambda_delta_prop, sd = omega) |> log()
    
    # Compute log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta_prop, lambda_v_vec, rho_v_mat)
    
    # Compute log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = lambda_delta_prop, accept = TRUE))
    } else {
      return(list(value = lambda_delta, accept = FALSE))
    }
  }
  
  update_lambda_v_vec_i_MH_MFM <- function(i, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega) {
    prop <- lambda_v_vec
    
    # Propose new value for element i using truncated normal > 0
    prop_i_new <- rtruncnorm(1, a = 0, b = Inf, mean = lambda_v_vec[i], sd = omega)
    prop[i] <- prop_i_new
    
    # Compute truncated normal densities for proposal ratio
    log_q_current_to_prop <- dtruncnorm(prop_i_new, a = 0, b = Inf, mean = lambda_v_vec[i], sd = omega) |> log()
    log_q_prop_to_current <- dtruncnorm(lambda_v_vec[i], a = 0, b = Inf, mean = prop_i_new, sd = omega) |> log()
    
    # Calculate log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, prop, rho_v_mat)
    
    # Log acceptance probability (including proposal ratio)
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = prop, accept = TRUE))
    } else {
      return(list(value = lambda_v_vec, accept = FALSE))
    }  
  }
  
  update_rho_v_mat_ij_MH_MFM <- function(i, j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_rho) {
    rho_v_pro <- rho_v_mat
    
    # Propose new rho[i,j] from truncated normal > 0
    rho_v_prop_ij <- rtruncnorm(1, a = 0, b = Inf, mean = rho_v_mat[i, j], sd = omega_rho[i, j])
    rho_v_pro[i, j] <- rho_v_prop_ij
    
    # Compute forward and reverse proposal log densities
    log_q_current_to_prop <- dtruncnorm(rho_v_prop_ij, a = 0, b = Inf, mean = rho_v_mat[i, j], sd = omega_rho[i, j]) |> log()
    log_q_prop_to_current <- dtruncnorm(rho_v_mat[i, j], a = 0, b = Inf, mean = rho_v_prop_ij, sd = omega_rho[i, j]) |> log()
    
    # Compute log posterior difference
    log_post_current <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat)
    log_post_prop    <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_pro)
    
    # Log acceptance probability including proposal ratio
    lpr <- (log_post_prop - log_post_current) + (log_q_prop_to_current - log_q_current_to_prop)
    
    if (log(runif(1)) < lpr) {
      return(list(value = rho_v_pro, accept = TRUE))
    } else {
      return(list(value = rho_v_mat, accept = FALSE))
    }  
  }
  
  
  # Initialize parameters
  lambda_eta   <- init_lambda_eta   %||% (a_eta_prime / b_eta_prime)
  lambda_w_vec <- init_lambda_w_vec %||% rep(1, r_d_lf)
  rho_w_mat    <- init_rho_w_mat    %||% matrix(1, r_d_lf, ncol(design_train_lf))
  
  lambda_delta <- init_lambda_delta %||% (a_delta_prime / b_delta_prime)
  lambda_v_vec <- init_lambda_v_vec %||% rep(1, r_d_discrep)
  rho_v_mat    <- init_rho_v_mat    %||% matrix(1, r_d_discrep, ncol(design_train_hf))
  
  # Storage for posterior samples
  n_save <- n_iter - burn_in
  
  lambda_eta_chain    <- numeric(n_save)
  lambda_w_vec_chain  <- matrix(NA, nrow = n_save, ncol = r_d_lf)
  rho_w_mat_chain     <- array(NA, dim = c(n_save, r_d_lf, ncol(design_train_lf)))
  
  lambda_delta_chain  <- numeric(n_save)
  lambda_v_vec_chain  <- matrix(NA, nrow = n_save, ncol = r_d_discrep)
  rho_v_mat_chain     <- array(NA, dim = c(n_save, r_d_discrep, ncol(design_train_hf)))
  
  log_lik_chain       <- numeric(n_save)
  
  
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
  
  set.seed(4273)
  for (iter in 1:n_iter) {
    
    # Update lambda_eta
    res <- update_lambda_eta_MH_MFM(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_lambda_eta)
    lambda_eta <- res$value
    # Update omega during burn-in
    if (iter <= burn_in) {
      accept_lambda_eta <- accept_lambda_eta + res$accept
      # Adapt every few steps
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_eta / update_omega_iter # only look at local acceptance
        omega_lambda_eta <- omega_lambda_eta * exp((rate - target_accept) / sqrt(iter) * learning_rate) # https://mvihola.github.io/docs/AdaptiveMCMC.jl/adapt/?utm_source=chatgpt.com
        accept_lambda_eta <- 0
        if (verbose) cat(paste0("\n\t\tUpdate omega_lambda_eta:\t", round(omega_lambda_eta, 4), "\t(Accept = ", round(rate, 2), ")\t(Current = ", round(lambda_eta, 2), ")\n"))
      }
    }
    
    # Update lambda_w_vec
    for (j in 1:r_d_lf) {
      res <- update_lambda_w_vec_i_MH_MFM(j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_lambda_w[j])
      lambda_w_vec <- res$value
      if (iter <= burn_in) {
        accept_lambda_w[j] <- accept_lambda_w[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_w[j] / update_omega_iter
          omega_lambda_w[j] <- omega_lambda_w[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_w[j] <- 0
          if (verbose) cat(paste0("\t\tUpdate omega_lambda_w[", j, "]:\t", round(omega_lambda_w[j], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(lambda_w_vec[j], 2), ")\n"))
        }
      }
    }
    
    # Update rho_w_mat
    for (i_rho in 1:r_d_lf) {
      for (j_rho in 1:ncol(design_train_lf)) {
        res <- update_rho_w_mat_ij_MH_MFM(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_rho_w)
        rho_w_mat <- res$value
        if (iter <= burn_in) {
          accept_rho_w[i_rho, j_rho] <- accept_rho_w[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho_w[i_rho, j_rho] / update_omega_iter
            omega_rho_w[i_rho, j_rho] <- omega_rho_w[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho_w[i_rho, j_rho] <- 0
            if (verbose) cat(paste0("\t\tUpdate omega_rho_w[", i_rho, ", ", j_rho, "]:\t\t", round(omega_rho_w[i_rho, j_rho], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(rho_w_mat[i_rho, j_rho], 2), ")\n"))
          }
        }
      }
    }
    
    # Update lambda_eta
    res <- update_lambda_delta_MH_MFM(lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_lambda_delta)
    lambda_delta <- res$value
    # Update omega during burn-in
    if (iter <= burn_in) {
      accept_lambda_delta <- accept_lambda_delta + res$accept
      # Adapt every few steps
      if (iter %% update_omega_iter == 0) {
        rate <- accept_lambda_delta / update_omega_iter # only look at local acceptance
        omega_lambda_delta <- omega_lambda_delta * exp((rate - target_accept) / sqrt(iter) * learning_rate) # https://mvihola.github.io/docs/AdaptiveMCMC.jl/adapt/?utm_source=chatgpt.com
        accept_lambda_delta <- 0
        if (verbose) cat(paste0("\n\t\tUpdate omega_lambda_delta:\t", round(omega_lambda_delta, 4), "\t(Accept = ", round(rate, 2), ")\t(Current = ", round(lambda_delta, 2), ")\n"))
      }
    }
    
    # Update lambda_v_vec
    for (j in 1:r_d_discrep) {
      res <- update_lambda_v_vec_i_MH_MFM(j, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_lambda_v[j])
      lambda_v_vec <- res$value
      if (iter <= burn_in) {
        accept_lambda_v[j] <- accept_lambda_v[j] + res$accept
        if (iter %% update_omega_iter == 0) {
          rate <- accept_lambda_v[j] / update_omega_iter
          omega_lambda_v[j] <- omega_lambda_v[j] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
          accept_lambda_v[j] <- 0
          if (verbose) cat(paste0("\t\tUpdate omega_lambda_v[", j, "]:\t", round(omega_lambda_v[j], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(lambda_v_vec[j], 2), ")\n"))
        }
      }
    }
    
    # Update rho_v_mat
    for (i_rho in 1:r_d_discrep) {
      for (j_rho in 1:ncol(design_train_lf)) {
        res <- update_rho_v_mat_ij_MH_MFM(i_rho, j_rho, lambda_eta, lambda_w_vec, rho_w_mat, lambda_delta, lambda_v_vec, rho_v_mat, omega_rho_v)
        rho_v_mat <- res$value
        if (iter <= burn_in) {
          accept_rho_v[i_rho, j_rho] <- accept_rho_v[i_rho, j_rho] + res$accept
          if (iter %% update_omega_iter == 0) {
            rate <- accept_rho_v[i_rho, j_rho] / update_omega_iter
            omega_rho_v[i_rho, j_rho] <- omega_rho_v[i_rho, j_rho] * exp((rate - target_accept) / sqrt(iter) * learning_rate)
            accept_rho_v[i_rho, j_rho] <- 0
            if (verbose) cat(paste0("\t\tUpdate omega_rho_v[", i_rho, ", ", j_rho, "]:\t\t", round(omega_rho_v[i_rho, j_rho], 4), "\t(Accpet = ", round(rate, 2), ")\t(Current = ", round(rho_v_mat[i_rho, j_rho], 2), ")\n"))
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
      
      log_lik_chain[iter - burn_in] <- log_posterior_mfm(lambda_eta, lambda_w_vec, rho_w_mat,
                                                         lambda_delta, lambda_v_vec, rho_v_mat)
    }
    
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

# ------------------------------------------------------------------------------
# @brief Checks convergence of MCMC chains using the Geweke diagnostic.
#
# This function computes the Geweke z-scores for each parameter chain contained
# within the MCMC output object. It supports scalar (vector), matrix, and 3D array
# chains and flags parameters that do not appear to have converged based on a
# specified z-score threshold. Requires the 'coda' package.
#
# @param mcmc_obj    List containing MCMC posterior chains under `post_chains`.
# @param z_threshold Numeric threshold for absolute Geweke z-score to flag
#                    potential non-convergence (default: 2).
#
# @return Invisibly returns a named numeric vector of Geweke z-scores for all parameters.
#         Prints a summary message highlighting any parameters failing convergence.
# ------------------------------------------------------------------------------
check_geweke_convergence <- function(mcmc_obj, frac1, frac2, z_threshold, verbose = T) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Please install the 'coda' package to use this function.")
  }
  library(coda)
  
  all_z <- c()
  chains_list <- list()
  
  post_chains <- mcmc_obj$post_chains
  
  for (name in names(post_chains)) {
    chain <- post_chains[[name]]
    
    # Vector chain
    if (is.vector(chain) && is.numeric(chain)) {
      z <- tryCatch(
        geweke.diag(mcmc(chain), frac1, frac2)$z,
        error = function(e) NA_real_
      )
      if (!is.na(z) && !is.nan(z)) {
        names(z) <- name
        all_z <- c(all_z, z)
        chains_list[[name]] <- mcmc(chain)
      } else {
        # message("Skipping parameter with NA/NaN Geweke z: ", name)
      }
      
      # Matrix chain (2D)
    } else if (is.matrix(chain)) {
      for (j in 1:ncol(chain)) {
        param_name <- paste0(name, "_", j)
        z <- tryCatch(
          geweke.diag(mcmc(chain[, j]), frac1, frac2)$z,
          error = function(e) NA_real_
        )
        if (!is.na(z) && !is.nan(z)) {
          names(z) <- param_name
          all_z <- c(all_z, z)
          chains_list[[param_name]] <- mcmc(chain[, j])
        } else {
          message("Skipping parameter with NA/NaN Geweke z: ", param_name)
        }
      }
      
      # 3D array chain
    } else if (length(dim(chain)) == 3) {
      d1 <- dim(chain)[2]
      d2 <- dim(chain)[3]
      for (i in 1:d1) {
        for (j in 1:d2) {
          param_name <- paste0(name, "_", i, "_", j)
          chain_vec <- chain[, i, j]
          z <- tryCatch(
            geweke.diag(mcmc(chain_vec), frac1, frac2)$z,
            error = function(e) NA_real_
          )
          if (!is.na(z) && !is.nan(z)) {
            names(z) <- param_name
            all_z <- c(all_z, z)
            chains_list[[param_name]] <- mcmc(chain_vec)
          } else {
            message("Skipping parameter with NA/NaN Geweke z: ", param_name)
          }
        }
      }
    } else {
      message("Skipping unknown chain format for: ", name)
    }
  }
  
  # Find non-converged parameters
  non_converged <- abs(all_z) > z_threshold
  
  if (any(non_converged, na.rm = TRUE)) {
    if (verbose) cat("\t\t⚠️ Geweke diagnostic suggests non-convergence for parameters:\n\t\t",
                     paste(names(all_z)[non_converged], collapse = ", "), "\n")
  } else {
    if (verbose) cat("\t\t✅ Geweke diagnostic: All parameters appear to have converged.\n")
  }
  
  return(all_z)
}

trim_chain <- function(chain, discard_frac) {
  new_chain <- chain
  for (name in names(chain$post_chains)) {
    param <- chain$post_chains[[name]]
    if (is.null(dim(param))) {
      keep_start <- floor(length(param) * discard_frac) + 1
      new_chain$post_chains[[name]] <- param[keep_start:length(param)]
    } else if (length(dim(param)) == 2) {
      keep_start <- floor(nrow(param) * discard_frac) + 1
      new_chain$post_chains[[name]] <- param[keep_start:nrow(param), , drop = F]
    } else if (length(dim(param)) == 3) {
      keep_start <- floor(dim(param)[1] * discard_frac) + 1
      new_chain$post_chains[[name]] <- param[keep_start:dim(param)[1], , , drop = F]
    }
    # Add handling for 3D if needed
  }
  return(new_chain)
}
check_geweke_multiple_chains <- function(chains, frac1 = 0.1, frac2 = 0.5, alpha = 0.05, max_discard_frac = 0.5) {
  total_tests <- 0
  retained_chains <- list()
  
  for (i in seq_along(chains)) {
    cat(paste0("\n\t📦 Chain ", i, " Geweke diagnostic — ", format(Sys.time(), "%b %d %X"), "\n"))
    chain <- chains[[i]]
    original_chain <- chain  # Keep original in case trimming fails
    
    # Count number of tests
    n_tests <- 0
    for (param in chain$post_chains) {
      if (is.null(dim(param))) {
        n_tests <- n_tests + 1
      } else if (length(dim(param)) == 2) {
        n_tests <- n_tests + ncol(param)
      } else if (length(dim(param)) == 3) {
        n_tests <- n_tests + prod(dim(param)[2:3])
      }
    }
    
    # z_threshold <- qnorm(1 - alpha / 2)
    z_threshold <- qnorm(1 - alpha / (2 * n_tests))  # Bonferroni-corrected threshold
    
    geweke_results <- check_geweke_convergence(chain, frac1, frac2, z_threshold, verbose = T)
    non_converged <- abs(geweke_results) > z_threshold
    failed_params <- names(geweke_results)[non_converged]
    
    # Check if key parameters failed
    if (any(grepl("^log_lik", failed_params)) || any(grepl("^lambda_eta", failed_params))) {
      cat("\t\t\tAttempting trimming...\n")
      # Attempt trimming up to max_discard_frac
      success <- FALSE
      for (trim_frac in seq(0.1, max_discard_frac, by = 0.1)) {
        chain_trimmed <- trim_chain(chain, trim_frac)
        geweke_results <- check_geweke_convergence(chain_trimmed, frac1, frac2, z_threshold, verbose = F)
        failed_params <- names(geweke_results)[abs(geweke_results) > z_threshold]
        if (!any(grepl("^log_lik", failed_params)) && !any(grepl("^lambda_eta", failed_params))) {
          cat(paste0("\t\t\tChain salvaged by discarding first ", trim_frac * 100, "% of samples.\n"))
          chain <- chain_trimmed
          success <- TRUE
          break
        }
      }
      if (!success) {
        cat("\t\t❌ Discarding chain — unable to salvage after trimming.\n")
        next
      }
    }
    
    retained_chains[[length(retained_chains) + 1]] <- chain
    total_tests <- total_tests + n_tests
  }
  
  cat(paste0("\n\tTotal Geweke tests: ", total_tests, "\n"))
  cat(paste0("\tRetained chains: ", length(retained_chains), " of ", length(chains), "\n\n"))
  
  return(retained_chains)
}

# Find the minimum number of iterations across chains
get_min_iter <- function(param_name, chains) {
  min(sapply(chains, function(chain) {
    param <- chain$post_chains[[param_name]]
    if (is.null(dim(param))) length(param)
    else dim(param)[1]
  }))
}

# Safe wrapper that truncates to the same number of rows
safe_mcmc_list <- function(param_name, chains) {
  min_iter <- get_min_iter(param_name, chains)
  mcmc.list(lapply(chains, function(chain) {
    param <- chain$post_chains[[param_name]]
    if (is.null(dim(param))) {
      mcmc(matrix(param[1:min_iter], ncol = 1), start = 1, thin = 1)
    } else if (length(dim(param)) == 2) {
      mcmc(param[1:min_iter, , drop = FALSE], start = 1, thin = 1)
    } else {
      stop("Use custom logic for 3D parameters like rho_mat")
    }
  }))
}

check_mcmc_convergence <- function(chains) {
  # --- lambda_eta ---
  mcmc_lambda_eta_list <- safe_mcmc_list("lambda_eta_chain", chains)
  gelman_lambda_eta <- gelman.diag(mcmc_lambda_eta_list, autoburnin = FALSE)
  
  # --- lambda_w_vec ---
  mcmc_lambda_w_list <- safe_mcmc_list("lambda_w_vec_chain", chains)
  gelman_lambda_w <- gelman.diag(mcmc_lambda_w_list)
  
  # --- rho_mat ---
  # Manually truncate and reshape
  min_iter <- min(sapply(chains, function(chain) dim(chain$post_chains$rho_mat_chain)[1]))
  mcmc_rho_list <- mcmc.list(lapply(chains, function(chain) {
    rho_chain <- chain$post_chains$rho_mat_chain[1:min_iter, , , drop = FALSE]
    r3 <- dim(rho_chain)[2]
    p <- dim(rho_chain)[3]
    mat <- matrix(NA, nrow = min_iter, ncol = r3 * p)
    for (i in 1:r3) {
      for (j in 1:p) {
        col_idx <- (i - 1) * p + j
        mat[, col_idx] <- rho_chain[, i, j]
      }
    }
    mcmc(mat, start = 1, thin = 1)
  }))
  gelman_rho <- gelman.diag(mcmc_rho_list)
  
  # --- log_lik ---
  mcmc_log_lik_list <- safe_mcmc_list("log_lik_chain", chains)
  gelman_log_lik <- gelman.diag(mcmc_log_lik_list, autoburnin = FALSE)
  
  list(
    lambda_eta = gelman_lambda_eta,
    lambda_w_vec = gelman_lambda_w,
    rho_mat = gelman_rho,
    log_lik = gelman_log_lik
  )
}

check_mcmc_convergence_mf <- function(chains) {
  # --- lambda_eta ---
  mcmc_lambda_eta_list <- safe_mcmc_list("lambda_eta_chain", chains)
  gelman_lambda_eta <- gelman.diag(mcmc_lambda_eta_list, autoburnin = FALSE)
  
  # --- lambda_delta ---
  mcmc_lambda_delta_list <- safe_mcmc_list("lambda_delta_chain", chains)
  gelman_lambda_delta <- gelman.diag(mcmc_lambda_delta_list, autoburnin = FALSE)
  
  # --- lambda_w_vec ---
  mcmc_lambda_w_list <- safe_mcmc_list("lambda_w_vec_chain", chains)
  gelman_lambda_w <- gelman.diag(mcmc_lambda_w_list)
  
  # --- lambda_v_vec ---
  mcmc_lambda_v_list <- safe_mcmc_list("lambda_v_vec_chain", chains)
  gelman_lambda_v <- gelman.diag(mcmc_lambda_v_list)
  
  # --- rho_w_mat ---
  # Manually truncate and reshape
  min_iter <- min(sapply(chains, function(chain) dim(chain$post_chains$rho_w_mat_chain)[1]))
  mcmc_rho_w_list <- mcmc.list(lapply(chains, function(chain) {
    rho_w_chain <- chain$post_chains$rho_w_mat_chain[1:min_iter, , , drop = F]
    r3 <- dim(rho_w_chain)[2]
    p <- dim(rho_w_chain)[3]
    mat <- matrix(NA, nrow = min_iter, ncol = r3 * p)
    for (i in 1:r3) {
      for (j in 1:p) {
        col_idx <- (i - 1) * p + j
        mat[, col_idx] <- rho_w_chain[, i, j]
      }
    }
    mcmc(mat, start = 1, thin = 1)
  }))
  gelman_rho_w <- gelman.diag(mcmc_rho_w_list)
  
  # --- rho_v_mat ---
  # Manually truncate and reshape
  min_iter <- min(sapply(chains, function(chain) dim(chain$post_chains$rho_v_mat_chain)[1]))
  mcmc_rho_v_list <- mcmc.list(lapply(chains, function(chain) {
    rho_v_chain_raw <- chain$post_chains$rho_v_mat_chain
    rho_v_chain <- array(rho_v_chain_raw, dim = c(dim(rho_v_chain_raw)))
    rho_v_chain <- rho_v_chain[1:min_iter, , , drop = FALSE]
    
    r3 <- dim(rho_v_chain)[2]
    p  <- dim(rho_v_chain)[3]
    
    mat <- matrix(NA, nrow = min_iter, ncol = r3 * p)
    for (i in 1:r3) {
      for (j in 1:p) {
        col_idx <- (i - 1) * p + j
        mat[, col_idx] <- rho_v_chain[, i, j]
      }
    }
    mcmc(mat, start = 1, thin = 1)
  }))
  gelman_rho_v <- gelman.diag(mcmc_rho_v_list)
  
  gelman_rho_v <- gelman.diag(mcmc_rho_v_list)
  
  
  # --- log_lik ---
  mcmc_log_lik_list <- safe_mcmc_list("log_lik_chain", chains)
  gelman_log_lik <- gelman.diag(mcmc_log_lik_list, autoburnin = FALSE)
  
  list(
    lambda_eta = gelman_lambda_eta,
    lambda_w_vec = gelman_lambda_w,
    rho_w_mat = gelman_rho_w,
    lambda_delta = gelman_lambda_delta,
    lambda_v_vec = gelman_lambda_v,
    rho_v_mat = gelman_rho_v,
    log_lik = gelman_log_lik
  )
}



check_rhat <- function(gelman_result, threshold = 1.1, name = "Unknown") {
  if (is.null(gelman_result$psrf)) {
    cat(paste0("\t\t⚠️ R̂ diagnostic failed for ", name, "\n"))
    return(invisible(NULL))
  }
  
  rhat_values <- gelman_result$psrf[, 1]
  param_names <- rownames(gelman_result$psrf)
  
  if (is.null(param_names)) {
    param_names <- paste0(name, "_", seq_along(rhat_values))
  }
  
  non_converged_idx <- which(rhat_values > threshold)
  converged_idx <- which(rhat_values <= threshold)
  
  n_total <- length(rhat_values)
  n_non_converged <- length(non_converged_idx)
  n_converged <- length(converged_idx)
  
  if (n_non_converged == 0) {
    cat(paste0("\t\t✅ Gelman-Rubin R̂: All ", n_total,
               " parameter(s) converged for ", name, "\n"))
  } else {
    cat(paste0("\t\t⚠️ Gelman-Rubin R̂: ", n_non_converged, " of ", n_total,
               " parameter(s) did NOT converge for ", name, ":\n\t\t"))
    cat(paste(param_names[non_converged_idx], collapse = ", "), "\n")
  }
  
  invisible(rhat_values)
}

check_mpsrf <- function(gelman_result, threshold = 1.1, name = "Unknown") {
  if (is.null(gelman_result$mpsrf)) {
    cat(paste0("\t\t⚠️ Multivariate R̂ diagnostic not available for ", name, "\n"))
    return(invisible(NULL))
  }
  
  mpsrf_value <- gelman_result$mpsrf
  
  if (mpsrf_value <= threshold) {
    cat(paste0("\t\t✅ Multivariate R̂: ", round(mpsrf_value, 3),
               " — converged for ", name, "\n"))
  } else {
    cat(paste0("\t\t⚠️ Multivariate R̂: ", round(mpsrf_value, 3),
               " — not converged for ", name, "\n"))
  }
  
  invisible(mpsrf_value)
}


# ------------------------------------------------------------------------------
# @brief Plots traceplot for a specified MCMC parameter chain.
#
# This function visualizes the sampling trace of a given parameter from the
# posterior chains. It supports vector and matrix chains. For matrix chains,
# an index specifying the column must be provided. The function raises errors
# if the parameter is missing or index is out of bounds.
#
# @param post_chains  Named list of posterior chains.
# @param param_name   Character string specifying the parameter to plot.
# @param index        Optional integer index for matrix columns (default: NULL).
#
# @return None (generates a plot).
# ------------------------------------------------------------------------------
plot_trace <- function(post_chains, param_name, index = NULL) {
  if (!param_name %in% names(post_chains)) {
    stop(paste("Parameter", param_name, "not found in post_chains."))
  }
  
  chain <- post_chains[[param_name]]
  
  if (is.matrix(chain)) {
    if (is.null(index)) {
      stop(paste("Parameter", param_name, "is a matrix. Please specify 'index'."))
    }
    if (index < 1 || index > ncol(chain)) {
      stop(paste("Index out of bounds for", param_name, "."))
    }
    plot(chain[, index], type = "l",
         main = paste0("Traceplot: ", param_name, "_", index),
         xlab = "Iteration", ylab = "Value")
  } else {
    if (!is.null(index)) {
      warning("Index ignored for vector chain.")
    }
    plot(chain, type = "l",
         main = paste0("Traceplot: ", param_name),
         xlab = "Iteration", ylab = "Value")
  }
}

library(FNN)

flatten_last_dim <- function(arr) {
  d <- dim(arr)
  dim(arr) <- c(prod(d[-length(d)]), d[length(d)])
  return(arr)
}

knn_interpolate_first_dim <- function(data_array, 
                                      coords_old, 
                                      coords_new, 
                                      k = 3) {
  # Get dimensions
  dims <- dim(data_array)
  ndim <- length(dims)
  
  # Validate
  if (nrow(coords_old) != dims[1]) {
    stop("Number of rows in coords_old must match size of first dimension of array")
  }
  
  # Output shape: replace first dim with new coords length
  out_dims <- dims
  out_dims[1] <- nrow(coords_new)
  
  # Prepare output
  interpolated <- array(NA, dim = out_dims)
  
  # Indices for the remaining dimensions
  rest_dims <- dims[-1]
  rest_indices <- expand.grid(lapply(rest_dims, seq_len))
  
  pb <- progress_bar$new(
    format = "  Interp. Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
    total = nrow(rest_indices),
    width = 70
  )
  
  # Iterate over all combinations of the other modes
  for (i in seq_len(nrow(rest_indices))) {
    # Extract indices for this combination
    idx_list <- as.list(rest_indices[i, ])
    
    # Build full index for subsetting original array
    full_index <- c(list(seq_len(dims[1])), idx_list)
    
    # Extract 1D slice along first dimension
    values <- do.call(`[`, c(list(data_array), full_index, list(drop = TRUE)))
    
    # Apply kNN regression
    knn_out <- knn.reg(train = coords_old, test = coords_new, y = values, k = k)
    
    # Build target index to assign output
    target_index <- c(list(seq_len(nrow(coords_new))), idx_list)
    interpolated <- do.call(`[<-`, c(list(interpolated), target_index, list(value = knn_out$pred)))
    
    pb$tick()
  }
  
  return(interpolated)
}

fast_knn_interpolate_first_dim <- function(data_array,
                                           coords_old,
                                           coords_new,
                                           k = 5,
                                           epsilon = 1e-8) {
  library(FNN)
  
  # Dimensions
  dims <- dim(data_array)
  n_new <- nrow(coords_new)
  n_old <- nrow(coords_old)
  other_dims <- prod(dims[-1])
  
  # Flatten array to [space, ...] → [n_old, N]
  data_mat <- matrix(data_array, nrow = n_old, ncol = other_dims)
  
  # Get k nearest neighbors and distances
  knn_res <- get.knnx(data = coords_old, query = coords_new, k = k)
  
  nn_idx <- knn_res$nn.index       # [n_new, k]
  nn_dist <- knn_res$nn.dist + epsilon  # avoid div by 0
  
  # Inverse distance weights
  weights <- 1 / nn_dist
  weights <- weights / rowSums(weights)  # normalize [n_new, k]
  
  # Interpolate each column (i.e., each flattened point in the rest of the array)
  interpolated_mat <- matrix(0, nrow = n_new, ncol = other_dims)
  
  for (j in seq_len(k)) {
    wj <- weights[, j]
    idxj <- nn_idx[, j]
    interpolated_mat <- interpolated_mat + wj * data_mat[idxj, ]
  }
  
  # Reshape back to array: [new_space, ...]
  new_dims <- dims
  new_dims[1] <- n_new
  interpolated_array <- array(interpolated_mat, dim = new_dims)
  
  return(interpolated_array)
}



evaluate_sf_mcmc <- function(untested_input,
                             mcmc_setup_outputs,
                             mcmc_outputs,
                             tucker_outputs,
                             n_post_draws,
                             design_train,
                             output_dim,
                             ranks,
                             hf = FALSE,
                             transform = NULL,
                             include_noise = FALSE,
                             obs,
                             verbose = T) {
  
  n_chains <- length(mcmc_outputs)
  
  r_d <- ranks[length(ranks)]
  
  n_iter  <- mcmc_outputs[[1]]$n_iter
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
    untested_input <- matrix(untested_input, nrow = 1)
    
    gamma_output_list <- list()
    pred_output_list <- list()
    count <- 1
    attempts <- 0
    max_attempts <- 3 * n_post_draws  # To prevent infinite loops
    
    if (verbose) {
      pb <- progress_bar$new(
        format = ifelse(hf,
                        "  HF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
                        "  LF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate"),
        total = n_post_draws,
        width = 70
      ) 
    }
    
    while (count <= n_post_draws && attempts < max_attempts) {
      attempts <- attempts + 1
      
      # Sample chain and draw index
      chain_i <- sample(1:n_chains, 1)
      chain_length <- length(mcmc_outputs[[chain_i]]$post_chains$lambda_eta_chain)
      draw_i <- sample(1:chain_length, 1)
      
      chain <- mcmc_outputs[[chain_i]]$post_chains
      
      lambda_eta   <- chain$lambda_eta_chain[draw_i]
      lambda_w_vec <- chain$lambda_w_vec_chain[draw_i, ]
      rho_w_mat    <- chain$rho_mat_chain[draw_i, , ]
      
      # Ensure rho_w_mat is a proper matrix
      if (is.null(dim(rho_w_mat))) {
        rho_w_mat <- matrix(rho_w_mat,
                            nrow = dim(chain$rho_mat_chain)[2],
                            ncol = dim(chain$rho_mat_chain)[3])
      }
      
      # Construct GP components
      A <- Sigma_general(design_train, design_train, lambda_w_vec, rho_w_mat, r3 = r_d) 
      B <-  (1 / lambda_eta * solve_C_t_C)
      
      V_11 <- A + B
      
      V_12 <- Sigma_general(design_train, untested_input, lambda_w_vec, rho_w_mat, r3 = r_d)
      
      V_21 <- t(V_12)
      V_22 <- Sigma_general(untested_input, untested_input, lambda_w_vec, rho_w_mat, r3 = r_d)
      
      # Predictive mean and covariance
      eigvals <- tryCatch(
        eigen(forceSymmetric(V_11), symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA
      )
      if (any(is.na(eigvals)) || any(eigvals <= 0)) {
        print("error: V_11 not pd")
        next  # skip non-PD draw
      }
      
      L <- chol(V_11)
      y <- backsolve(L, forwardsolve(t(L), gamma_hat_1))
      Y <- backsolve(L, forwardsolve(t(L), V_12))
      
      mu <- V_21 %*% y
      sigma <- as.matrix(V_22 - V_21 %*% Y)
      
      # in case of negative diagonals due to numerical imprecision
      if (any(diag(sigma) < 0)) {
        print("update sigma")
        diag(sigma)[diag(sigma) < 0] <- 1e-12
      }
      
      # print(sigma)
      
      
      # Handle scalar sigma case
      if (length(sigma) == 1 && isTRUE(as.numeric(sigma) < 0)) {
        print("error: sigma not pd")
        next  # skip invalid draw
      }
      
      # Check if sigma is positive definite
      eigvals <- tryCatch(
        eigen(forceSymmetric(sigma), symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA
      )
      if (any(is.na(eigvals)) || any(eigvals <= 0)) {
        print("error: sigma not pd")
        next  # skip non-PD draw
      }
      
      # Draw from predictive distribution
      gamma_star <- t(rmvnorm(1, mu, sigma))
      
      # Optional debug print
      # cat(gamma_star, "\n")
      
      gamma_part <- B_eta %*% G_unfold %*% gamma_star
      noise_part <- rnorm(n, 0, sd = sqrt(1 / lambda_eta))
      
      # cat("Norms:\n")
      # cat("  gamma: ", norm(gamma_part, "2"), "\n")
      # cat("  noise: ", norm(noise_part, "2"), "\n")
      
      # # Final prediction
      # pred <- gamma_part
      # 
      # if (include_noise) {
      #   pred <- pred + noise_part
      # }
      
      # Apply transformation if needed
      if (transform == "logistic") {
        pred <- logistic(pred)
        pred <- pmin(pmax(pred, 1e-2), 1 - 1e-2)
      } else if (transform == "zero-one") {
        gamma_part <- pmin(pmax(gamma_part, 0), 1)
        if (include_noise) {
          pred <- gamma_part + noise_part
          pred <- pmin(pmax(pred, 0), 1)
        } else {
          pred <- gamma_part
        }
      } else if (transform == "logistic-soft-clip") {
        # gamma_part <- pmin(pmax(gamma_part, min(output_hf_4d_trans_scaled)), max(output_hf_4d_trans_scaled))
        if (include_noise) {
          pred <- gamma_part + noise_part
          # pred <- pmin(pmax(pred, min(output_hf_4d_trans_scaled)), max(output_hf_4d_trans_scaled))
        } else {
          pred <- gamma_part
        }
        pred <- logistic_softclip(unscaled_transformed(pred, max(output_hf_4d_trans)), epsilon = 0.0015)
      } else if (transform == "none") {
        if (include_noise) {
          pred <- gamma_part + noise_part
        } else {
          pred <- gamma_part
        }
      }
      
      # Store valid prediction
      gamma_output_list[[count]] <- gamma_star
      pred_output_list[[count]] <- pred
      if (verbose) pb$tick()
      count <- count + 1
    }
    
    # Combine into matrix
    if (length(pred_output_list) == 0) {
      stop("❌ No valid posterior draws were collected.")
    }
    
    gamma_output_mat <- do.call(cbind, gamma_output_list)
    pred_output_mat <- do.call(cbind, pred_output_list)
    return(list(
      pred_output_mat = pred_output_mat,
      gamma_output_mat = gamma_output_mat
    ))
  }
  
  out <- post_pred_sf(untested_input, n_post_draws)
  
  gamma_draws <- out$gamma_output_mat
  post_pred_draws_mat <- out$pred_output_mat
  
  post_pred_draws_arr <- array(post_pred_draws_mat, c(output_dim, n_post_draws))
  
  # Interpolate posterioir preds
  if (!hf) {
    post_pred_draws_arr <- fast_knn_interpolate_first_dim(post_pred_draws_arr, coord_lf, coord_hf, k = 1)
    post_pred_draws_mat <- flatten_last_dim(post_pred_draws_arr)
  }
  
  spatiotemp_modes <- 1:(length(dim(post_pred_draws_arr)) - 1)
  
  means <- rowMeans(post_pred_draws_mat)
  post_pred_mean <- array(means, dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  sds   <- rowSds(post_pred_draws_mat) 
  post_pred_sd <- array(sds, dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  q_vals <- rowQuantiles(post_pred_draws_mat, probs = c(0.025, 0.975))
  post_pred_q_low  <- array(q_vals[, 1], dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  post_pred_q_high <- array(q_vals[, 2], dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  mse <- (post_pred_mean - obs)^2
  mbe <- apply((post_pred_draws_arr - as.numeric(obs)), spatiotemp_modes, mean)
  mae <- apply((abs(post_pred_draws_arr - as.numeric(obs))), spatiotemp_modes, mean)
  
  cover <- (post_pred_q_low <= obs) & (obs <= post_pred_q_high)
  
  return(list(
    gamma_draws = gamma_draws,
    # post_pred_draws = post_pred_draws_arr,
    
    post_pred_mean = post_pred_mean,
    post_pred_q_low = post_pred_q_low,
    post_pred_q_high = post_pred_q_high,
    
    mse = mse,
    mbe = mse,
    mae = mae,
    
    sd = post_pred_sd,
    
    cover = cover
  ))
}

drop_first_dim <- function(x) {
  dim_x <- dim(x)
  if (!is.null(dim_x) && dim_x[1] == 1) {
    array(x[1, , , drop = FALSE], dim = dim_x[-1])
  } else {
    x
  }
}

evaluate_mf_mcmc <- function(untested_input, 
                             tucker_outputs_lf, tucker_outputs_discrep,
                             interpolate_bases_outputs,
                             mcmc_lf_setup_outputs, mcmc_mf_setup_outputs,
                             mf_mcmc_outputs, 
                             n_post_draws,
                             output_dim,
                             transform = NULL,
                             include_noise = FALSE,
                             obs,
                             verbose = T) {
  
  n_chains <- length(mf_mcmc_outputs)
  
  # B_eta_tilde <- interpolate_bases_outputs$B_eta_tilde
  
  solve_C_t_C <- mcmc_lf_setup_outputs$solve_C_t_C
  solve_K_t_K <- mcmc_mf_setup_outputs$solve_K_t_K
  
  gamma_hat_lf <- mcmc_lf_setup_outputs$gamma_hat
  gamma_hat_hf <- mcmc_mf_setup_outputs$gamma_hat_hf
  zeta_hat     <- mcmc_mf_setup_outputs$zeta_hat
  hat_vals     <- c(as.vector(gamma_hat_lf), as.vector(gamma_hat_hf), as.vector(zeta_hat))
  
  G <- tucker_outputs_lf$core_tensor
  G_discrep <- tucker_outputs_discrep$core_tensor
  
  G_unfold <- t(k_unfold(as.tensor(G), m = length(dim(G)))@data)
  G_unfold_discrep <- t(k_unfold(as.tensor(G_discrep), m = length(dim(G_discrep)))@data)
  
  # D_delta <- do.call(cbind, lapply(tucker_outputs_discrep$basis_matrices, as.vector))
  
  D_G <- matrix(0, nrow = s_hf * t_hf, ncol = ncol(G_unfold_discrep))
  counter <- 1
  if (verbose) {
    pb <- progress_bar$new(
      format = "  Matrix Calc [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
      total = r_y_discrep * r_m_discrep * r_s_discrep,
      width = 70
    )
  }
  for (y in 1:r_y_discrep) {
    for (m in 1:r_m_discrep) {
      for (s in 1:r_s_discrep) {
        basis_matrix <- calculate_basis(U_list = tuck_discrep$factor_matrices,
                                        idx = c(s, m, y))
        D_G <- D_G + as.vector(basis_matrix) %*% G_unfold_discrep[counter, , drop = FALSE]
        counter <- counter + 1
        
        if (verbose) pb$tick()
      }
      gc()
    }
  }
  
  # Prepare chain info for sampling
  chain_lengths <- sapply(mf_mcmc_outputs, function(chain) {
    length(chain$post_chains$lambda_eta_chain)
  })
  
  # Prediction function
  post_pred_mf <- function(theta_star, n_post_draws) {
    theta_star <- matrix(theta_star, nrow = 1)
    
    pred_output_list <- list()
    count <- 1
    attempts <- 0
    max_attempts <- 3 * n_post_draws  # Avoid infinite loops
    
    if (verbose) {
      pb <- progress_bar$new(
        format = "  MF Post. Preds. [:bar] :current/:total (:percent) ETA: :eta Rate: :tick_rate",
        total = n_post_draws,
        width = 70
      ) 
    }
    
    while (count <= n_post_draws && attempts < max_attempts) {
      attempts <- attempts + 1
      # print(count / attempts)
      
      # Sample chain and draw index
      chain_i <- sample(1:n_chains, 1)
      draw_i <- sample(1:chain_lengths[chain_i], 1)
      chain <- mf_mcmc_outputs[[chain_i]]$post_chains
      
      # Extract parameters
      lambda_eta   <- chain$lambda_eta_chain[draw_i]
      lambda_w_vec <- chain$lambda_w_vec_chain[draw_i, ]
      rho_w_mat <- drop_first_dim(chain$rho_w_mat_chain[draw_i, , , drop = FALSE])
      
      lambda_delta   <- chain$lambda_delta_chain[draw_i]
      lambda_v_vec   <- chain$lambda_v_vec_chain[draw_i, ]
      rho_v_mat <- drop_first_dim(chain$rho_v_mat_chain[draw_i, , , drop = FALSE])
      
      # Construct covariance blocks
      V_11 <- Sigma_general(design_train_lf, design_train_lf, lambda_w_vec, rho_w_mat, r_d_lf)
      V_12 <- Sigma_general(design_train_lf, design_train_hf, lambda_w_vec, rho_w_mat, r_d_lf)
      V_21 <- t(V_12)
      V_22 <- Sigma_general(design_train_hf, design_train_hf, lambda_w_vec, rho_w_mat, r_d_lf)
      V <- rbind(cbind(V_11, V_12), cbind(V_21, V_22))
      
      A <- bdiag(V,
                 Sigma_general(design_train_hf, design_train_hf, lambda_v_vec, rho_v_mat, r_d_discrep))
      B <- bdiag(1 / lambda_eta * solve_C_t_C, 1 / lambda_delta * solve_K_t_K)
      C <- bdiag(
        rbind(
          Sigma_general(design_train_lf, theta_star, lambda_w_vec, rho_w_mat, r_d_lf),
          Sigma_general(design_train_hf, theta_star, lambda_w_vec, rho_w_mat, r_d_lf)
        ),
        Sigma_general(design_train_hf, theta_star, lambda_v_vec, rho_v_mat, r_d_discrep)
      )
      D <- bdiag(
        Sigma_general(theta_star, theta_star, lambda_w_vec, rho_w_mat, r_d_lf),
        Sigma_general(theta_star, theta_star, lambda_v_vec, rho_v_mat, r_d_discrep)
      )
      
      
      # kappa_val <- kappa(A + B)
      # if (kappa_val > 1e-10) {
      #   print("\n\t\t⚠️ A + B is ill-conditioned")
      #   B <- B + diag(1e-10, nrow(B))
      # }
      
      
      L <- chol(A + B)
      y <- backsolve(L, forwardsolve(t(L), hat_vals))
      Y <- backsolve(L, forwardsolve(t(L), C))
      
      mu <- t(C) %*% y
      sigma <- as.matrix(D - t(C) %*% Y)
      
      
      # print(sigma)
      
      # Check PD
      eigvals <- eigen(forceSymmetric(sigma), symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals <= 0)) {
        print("\n\t\t⚠️ sigma not PD")
        next
      }
      
      out <- rmvnorm(1, mu, sigma)
      gamma_star <- out[1:r_d_lf]
      zeta_star  <- out[-(1:r_d_lf)]
      
      # cat(gamma_star, "\n")
      
      gamma_part <- interpolate_bases_outputs$B_eta_tile_x_G_unfold %*% gamma_star
      zeta_part  <- D_G %*% zeta_star
      noise_part <- rnorm(s_hf * t_hf, 0, sd = sqrt(1 / lambda_delta))
      # + rnorm(s_hf * t_hf, 0, sd = sqrt(1 / lambda_eta)) + 
      
      
      # cat("Norms:\n")
      # cat("  gamma: ", norm(gamma_part, "2"), "\n")
      # cat("  zeta:  ", norm(zeta_part, "2"), "\n")
      # cat("  noise: ", norm(noise_part, "2"), "\n")
      
      
      pred <- gamma_part + zeta_part
      
      if (include_noise) {
        pred <- pred + noise_part
      }
      
      # Apply transformation if needed
      if (transform == "logistic") {
        pred <- logistic(pred)
        pred <- pmin(pmax(pred, 1e-2), 1 - 1e-2)
      } else if (transform == "zero-one") {
        # What I had running current 01
        # gamma_plus_zeta_part <- pmin(pmax(gamma_part + zeta_part, 0), 1)
        # pred <- gamma_plus_zeta_part + noise_part
        # pred <- pmin(pmax(pred, 0), 1)
        
        gamma_plus_zeta_part <- pmin(pmax(gamma_part + zeta_part, 0), 1)
        if (include_noise) {
          pred <- gamma_plus_zeta_part + noise_part
          pred <- pmin(pmax(pred, 0), 1)
        } else {
          pred <- gamma_plus_zeta_part
        }
        
      } else if (transform == "logistic-soft-clip") {
        # gamma_plus_zeta_part <- pmin(pmax(gamma_part + zeta_part, min(output_hf_4d_trans)), max(output_hf_4d_trans))
        if (include_noise) {
          pred <- gamma_part + zeta_part + noise_part
          # pred <- gamma_plus_zeta_part + noise_part
          # pred <- pmin(pmax(pred, min(output_hf_4d_trans)), max(output_hf_4d_trans))
        } else {
          pred <- gamma_plus_zeta_part
        }
        pred <- logistic_softclip(unscaled_transformed(pred, max(output_hf_4d_trans)), epsilon = 0.0015)
      }
      
      # Store valid prediction
      pred_output_list[[count]] <- pred
      if (verbose) pb$tick()
      count <- count + 1
    }
    
    # Combine successful draws into matrix
    pred_output_mat <- do.call(cbind, pred_output_list)
    
    return(pred_output_mat)
  }
  
  # Run post pred
  post_pred_draws_mat <- post_pred_mf(untested_input, n_post_draws)
  
  post_pred_draws_arr <- array(post_pred_draws_mat, c(output_dim, n_post_draws))
  
  spatiotemp_modes <- 1:(length(dim(post_pred_draws_arr)) - 1)
  
  means <- rowMeans(post_pred_draws_mat)
  post_pred_mean <- array(means, dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  sds   <- rowSds(post_pred_draws_mat) 
  post_pred_sd <- array(sds, dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  q_vals <- rowQuantiles(post_pred_draws_mat, probs = c(0.025, 0.975))
  post_pred_q_low  <- array(q_vals[, 1], dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  post_pred_q_high <- array(q_vals[, 2], dim = dim(post_pred_draws_arr)[spatiotemp_modes])
  
  mse <- (post_pred_mean - obs)^2
  mbe <- apply((post_pred_draws_arr - as.numeric(obs)), spatiotemp_modes, mean)
  mae <- apply((abs(post_pred_draws_arr - as.numeric(obs))), spatiotemp_modes, mean)
  
  cover <- (post_pred_q_low <= obs) & (obs <= post_pred_q_high)
  
  return(list(
    # gamma_draws = gamma_draws,
    # post_pred_draws = post_pred_draws_arr,
    
    post_pred_mean = post_pred_mean,
    post_pred_q_low = post_pred_q_low,
    post_pred_q_high = post_pred_q_high,
    
    mse = mse,
    mbe = mse,
    mae = mae,
    
    sd = post_pred_sd,
    
    cover = cover
  ))
}

floor_dec <- function(x, digits = 0) {
  factor <- 10^digits
  floor(x * factor) / factor
}

ceiling_dec <- function(x, digits = 0) {
  factor <- 10^digits
  ceiling(x * factor) / factor
}

# ------------------------------------------------------------------------------
# @brief Computes cross-validation statistics for posterior predictions.
#
# This function evaluates posterior predictive performance against held-out 
# observations. It computes overall, spatial, and temporal metrics including RMSE, 
# MAE, MBE, posterior standard deviation, and coverage probability. Optionally, 
# for low-fidelity predictions, it interpolates to the high-fidelity spatial grid 
# before computing metrics.
#
# Global dependencies:
#   - Uses: coord_lf, coord_hf, s_hf, t_lf
#   - Requires: interpolate_locations(), rowMeans2(), colMeans2()
#
# @param post_pred_results  A list containing:
#                           - post_pred_mean   : [space × time] matrix of predictive means
#                           - post_pred_sd     : [space × time] matrix of predictive standard deviations
#                           - post_pred_q_low  : [space × time] matrix of lower quantiles (2.5%)
#                           - post_pred_q_high : [space × time] matrix of upper quantiles (97.5%)
#                           - post_pred_draws  : Matrix of raw posterior draws
# @param obs                Matrix of observed values at high-fidelity resolution (same dims as mean).
# @param lf                 Logical flag; if TRUE, low-fidelity predictions are interpolated to high-fidelity grid.
#
# @return A named list with three sublists:
#         - overall  : RMSE, MAE, MBE, STDEV, COVER (averaged over space × time)
#         - spatial  : Metrics averaged over time, one value per spatial location
#         - temporal : Metrics averaged over space, one value per time point
# ------------------------------------------------------------------------------
# cv_stats <- function(post_pred_results, obs, lf = FALSE, r = NA) {
#   post_pred_mean   <- post_pred_results$post_pred_mean  # [space × time]
#   post_pred_sd     <- post_pred_results$post_pred_sd
#   post_pred_q_low  <- post_pred_results$post_pred_q_low
#   post_pred_q_high <- post_pred_results$post_pred_q_high
#   
#   if (!is.na(r)) {
#     post_pred_q_low <- floor_dec(post_pred_q_low, r)
#     post_pred_q_high <- ceiling_dec(post_pred_q_high, r)
#   }
#   
#   out_dims <- dim(post_pred_mean)         # e.g. c(s_lf, d1, d2, ...)
#   other_dims <- out_dims[-1]
#   combinations <- expand.grid(lapply(other_dims, seq_len))
# 
#   if (lf == TRUE) {
#     post_pred_mean_interp   <- array(NA, dim = c(s_hf, other_dims))
#     post_pred_q_low_interp  <- array(NA, dim = c(s_hf, other_dims))
#     post_pred_q_high_interp <- array(NA, dim = c(s_hf, other_dims))
#     
#     # Interpolate each slice
#     for (i in 1:nrow(combinations)) {
#       idx <- as.integer(combinations[i, ])
#       
#       full_idx <- c(list(TRUE), as.list(idx))  # TRUE = all rows along spatial dim
#       
#       # Extract vector along spatial dim
#       mean_slice <- do.call(`[`, c(list(post_pred_mean), full_idx))
#       low_slice  <- do.call(`[`, c(list(post_pred_q_low), full_idx))
#       high_slice <- do.call(`[`, c(list(post_pred_q_high), full_idx))
#       
#       # Interpolate
#       interp_mean <- interpolate_locations(coord_lf, mean_slice, coord_hf, k = 5)
#       interp_low  <- interpolate_locations(coord_lf, low_slice,  coord_hf, k = 5)
#       interp_high <- interpolate_locations(coord_lf, high_slice, coord_hf, k = 5)
#       
#       # Set result
#       assign_idx <- c(list(TRUE), as.list(idx))
#       post_pred_mean_interp   <- do.call(`[<-`, c(list(post_pred_mean_interp), assign_idx, list(interp_mean)))
#       post_pred_q_low_interp  <- do.call(`[<-`, c(list(post_pred_q_low_interp),  assign_idx, list(interp_low)))
#       post_pred_q_high_interp <- do.call(`[<-`, c(list(post_pred_q_high_interp), assign_idx, list(interp_high)))
#     }
#     
#     post_pred_mean <- post_pred_mean_interp
#     post_pred_q_low <- post_pred_q_low_interp
#     post_pred_q_high <- post_pred_q_high_interp
#   }
#   
#   diff <- post_pred_mean - obs
#   
#   # Overall metrics
#   RMSE  <- sqrt(mean(diff^2, na.rm = TRUE))
#   MAE   <- mean(abs(diff), na.rm = TRUE)
#   MBE   <- mean(diff, na.rm = TRUE)
#   STDEV <- mean(post_pred_sd, na.rm = TRUE)
#   COVER <- mean((obs >= post_pred_q_low) & (obs <= post_pred_q_high), na.rm = TRUE)
#   
#   if (length(dim(obs)) == 2) {
#     # Spatial metrics (average over time for each space point)
#     RMSE.s  <- sqrt(rowMeans2(diff^2, na.rm = TRUE))
#     MAE.s   <- rowMeans2(abs(diff), na.rm = TRUE)
#     MBE.s   <- rowMeans2(diff, na.rm = TRUE)
#     STDEV.s <- rowMeans2(post_pred_sd, na.rm = TRUE)
#     COVER.s <- rowMeans2((obs >= post_pred_q_low) & (obs <= post_pred_q_high), na.rm = TRUE)
#     
#     # Temporal metrics (average over space for each time point)
#     RMSE.t  <- sqrt(colMeans2(diff^2, na.rm = TRUE))
#     MAE.t   <- colMeans2(abs(diff), na.rm = TRUE)
#     MBE.t   <- colMeans2(diff, na.rm = TRUE)
#     STDEV.t <- colMeans2(post_pred_sd, na.rm = TRUE)
#     COVER.t <- colMeans2((obs >= post_pred_q_low) & (obs <= post_pred_q_high), na.rm = TRUE) 
#     
#     return(list(
#       overall = list(RMSE = RMSE, MAE = MAE, MBE = MBE, STDEV = STDEV, COVER = COVER),
#       spatial = list(RMSE = RMSE.s, MAE = MAE.s, MBE = MBE.s, STDEV = STDEV.s, COVER = COVER.s),
#       temporal = list(RMSE = RMSE.t, MAE = MAE.t, MBE = MBE.t, STDEV = STDEV.t, COVER = COVER.t)
#     ))
#   } else if (length(dim(obs)) == 3) {
#     # Spatial metrics (average over time and space for each item in first dimension)
#     RMSE.s  <- sqrt(apply(diff^2, MARGIN = 1, mean, na.rm = TRUE))
#     MAE.s   <- apply(abs(diff), MARGIN = 1, mean, na.rm = TRUE)
#     MBE.s   <- apply(diff, MARGIN = 1, mean, na.rm = TRUE)
#     STDEV.s <- apply(post_pred_sd, MARGIN = 1, mean, na.rm = TRUE)
#     COVER.s <- apply((obs >= post_pred_q_low) & (obs <= post_pred_q_high), 
#                      MARGIN = 1, mean, na.rm = TRUE)
#     
#     RMSE.m  <- sqrt(apply(diff^2, MARGIN = 2, mean, na.rm = TRUE))
#     MAE.m   <- apply(abs(diff), MARGIN = 2, mean, na.rm = TRUE)
#     MBE.m   <- apply(diff, MARGIN = 2, mean, na.rm = TRUE)
#     STDEV.m <- apply(post_pred_sd, MARGIN = 2, mean, na.rm = TRUE)
#     COVER.m <- apply((obs >= post_pred_q_low) & (obs <= post_pred_q_high), 
#                      MARGIN = 2, mean, na.rm = TRUE)
#     
#     RMSE.y  <- sqrt(apply(diff^2, MARGIN = 3, mean, na.rm = TRUE))
#     MAE.y   <- apply(abs(diff), MARGIN = 3, mean, na.rm = TRUE)
#     MBE.y   <- apply(diff, MARGIN = 3, mean, na.rm = TRUE)
#     STDEV.y <- apply(post_pred_sd, MARGIN = 3, mean, na.rm = TRUE)
#     COVER.y <- apply((obs >= post_pred_q_low) & (obs <= post_pred_q_high), 
#                      MARGIN = 3, mean, na.rm = TRUE)
#     
#     return(list(
#       overall = list(RMSE = RMSE, MAE = MAE, MBE = MBE, STDEV = STDEV, COVER = COVER),
#       spatial = list(RMSE = RMSE.s, MAE = MAE.s, MBE = MBE.s, STDEV = STDEV.s, COVER = COVER.s),
#       month = list(RMSE = RMSE.m, MAE = MAE.m, MBE = MBE.m, STDEV = STDEV.m, COVER = COVER.m),
#       year = list(RMSE = RMSE.y, MAE = MAE.y, MBE = MBE.y, STDEV = STDEV.y, COVER = COVER.y)
#     ))
#   } else {
#     warning("\t\tobs has incorrect dims")
#   }
# }

cv_stats <- function(post_pred_results, obs, lf = FALSE, r = NA) {
  post_pred_mean   <- post_pred_results$post_pred_mean  # [space × time]
  post_pred_sd     <- post_pred_results$post_pred_sd
  post_pred_q_low  <- post_pred_results$post_pred_q_low
  post_pred_q_high <- post_pred_results$post_pred_q_high
  
  if (!is.na(r)) {
    post_pred_q_low <- floor_dec(post_pred_q_low, r)
    post_pred_q_high <- ceiling_dec(post_pred_q_high, r)
  }
  
  out_dims <- dim(post_pred_mean)         # e.g. c(s_lf, d1, d2, ...)
  other_dims <- out_dims[-1]
  combinations <- expand.grid(lapply(other_dims, seq_len))
  
  if (lf == TRUE) {
    post_pred_mean_interp   <- array(NA, dim = c(s_hf, other_dims))
    post_pred_sd_interp     <- array(NA, dim = c(s_hf, other_dims))
    post_pred_q_low_interp  <- array(NA, dim = c(s_hf, other_dims))
    post_pred_q_high_interp <- array(NA, dim = c(s_hf, other_dims))
    
    # Interpolate each slice
    for (i in 1:nrow(combinations)) {
      idx <- as.integer(combinations[i, ])
      
      full_idx <- c(list(TRUE), as.list(idx))  # TRUE = all rows along spatial dim
      
      # Extract vector along spatial dim
      mean_slice <- do.call(`[`, c(list(post_pred_mean), full_idx))
      sd_slice   <- do.call(`[`, c(list(post_pred_sd), full_idx))
      low_slice  <- do.call(`[`, c(list(post_pred_q_low), full_idx))
      high_slice <- do.call(`[`, c(list(post_pred_q_high), full_idx))
      
      # Interpolate
      interp_mean <- interpolate_locations(coord_lf, mean_slice, coord_hf, k = 1)
      interp_sd   <- interpolate_locations(coord_lf, sd_slice, coord_hf, k = 1)
      interp_low  <- interpolate_locations(coord_lf, low_slice,  coord_hf, k = 1)
      interp_high <- interpolate_locations(coord_lf, high_slice, coord_hf, k = 1)
      
      # Set result
      assign_idx <- c(list(TRUE), as.list(idx))
      post_pred_mean_interp   <- do.call(`[<-`, c(list(post_pred_mean_interp), assign_idx, list(interp_mean)))
      post_pred_sd_interp     <- do.call(`[<-`, c(list(post_pred_sd_interp), assign_idx, list(interp_sd)))
      post_pred_q_low_interp  <- do.call(`[<-`, c(list(post_pred_q_low_interp),  assign_idx, list(interp_low)))
      post_pred_q_high_interp <- do.call(`[<-`, c(list(post_pred_q_high_interp), assign_idx, list(interp_high)))
    }
    
    post_pred_mean <- post_pred_mean_interp
    post_pred_sd   <- post_pred_sd_interp
    post_pred_q_low <- post_pred_q_low_interp
    post_pred_q_high <- post_pred_q_high_interp
  }
  
  diff <- post_pred_mean - obs
  
  # Overall metrics
  return(list(
    MSE  = diff^2,
    MAE   = abs(diff),
    MBE   = diff,
    STDEV = post_pred_sd,
    COVER = (obs >= post_pred_q_low) & (obs <= post_pred_q_high)
    
  ))
}


# ------------------------------------------------------------------------------
# @brief Aggregates cross-validation metrics across multiple saved runs.
#
# This function loads and aggregates CV statistics (overall, spatial, temporal) 
# from a set of `.RData` files saved in a specified folder. Each file is expected 
# to contain one of: `cv_stats_lf`, `cv_stats_hf`, or `cv_stats_mf` depending on 
# the fidelity level. It compiles run-level summaries and returns a list of 
# vectors and matrices for each evaluation metric.
#
# Global dependencies:
#   - Uses: s_hf (global number of high-fidelity spatial points)
#
# @param folder_dir  String path to directory containing saved `.RData` files with CV stats.
# @param fidelity    String indicating model fidelity: "lf" (low-fidelity), 
#                    "hf" (high-fidelity only), or "mf" (multi-fidelity).
#
# @return A named list containing:
#         - run_ids   : Character vector of run identifiers (currently empty).
#         - RMSEs     : Vector of overall RMSE values for each run.
#         - MAEs      : Vector of overall MAE values for each run.
#         - MBEs      : Vector of overall MBE values for each run.
#         - STDEVs    : Vector of overall predictive standard deviations.
#         - COVERs    : Vector of overall 95% coverage proportions.
#         - RMSEs_s   : Matrix [space × runs] of spatial RMSE values.
#         - MAEs_s    : Matrix [space × runs] of spatial MAE values.
#         - MBEs_s    : Matrix [space × runs] of spatial MBE values.
#         - STDEVs_s  : Matrix [space × runs] of spatial predictive std deviations.
#         - COVERs_s  : Matrix [space × runs] of spatial coverage proportions.
#         - RMSEs_t   : Matrix [time × runs] of temporal RMSE values.
#         - MAEs_t    : Matrix [time × runs] of temporal MAE values.
#         - MBEs_t    : Matrix [time × runs] of temporal MBE values.
#         - STDEVs_t  : Matrix [time × runs] of temporal predictive std deviations.
#         - COVERs_t  : Matrix [time × runs] of temporal coverage proportions.
# ------------------------------------------------------------------------------
get_all_cv_stats <- function(folder_dir, fidelity) {
  files <- list.files(folder_dir, full.names = TRUE)
  
  # Extract the numeric part (before ".Rdata")
  file_nums <- as.numeric(gsub(".*/|\\.Rdata", "", files))
  
  # Order by number
  files <- files[order(file_nums)]
  
  print(files)
  n_runs <- length(files)
  
  load(files[1])  # loads cv_stats
  if (fidelity == "lf") {
    cv_stats <- cv_stats_lf
    ranks <- ranks_lf
  } else if (fidelity == "hf") {
    cv_stats <- cv_stats_hf
    ranks <- ranks_hf
  } else {
    cv_stats <- cv_stats_mf
    ranks <- ranks_lf
    ranks_discrep <- ranks_discrep
  }
  n_spatial <- length(cv_stats$spatial$STDEV)
  n_month <- length(cv_stats$month$RMSE)
  n_year <- length(cv_stats$year$RMSE)
  
  
  # Ranks
  ranks_mat <- matrix(NA, nrow = length(ranks), ncol = n_runs)
  
  if (fidelity == "mf") {
    ranks_discrep_mat <- matrix(NA, nrow = length(ranks_discrep), ncol = n_runs)
  } else {
    ranks_discrep_mat <- NULL
  }
  
  # Overall metrics (vectors)
  RMSEs   <- numeric(n_runs)
  MAEs    <- numeric(n_runs)
  MBEs    <- numeric(n_runs)
  STDEVs  <- numeric(n_runs)
  COVERs  <- numeric(n_runs)
  run_ids <- character(n_runs)
  
  # Spatial metrics (matrices: rows = space, cols = runs)
  RMSEs_s  <- matrix(NA, nrow = s_hf, ncol = n_runs)
  MAEs_s   <- matrix(NA, nrow = s_hf, ncol = n_runs)
  MBEs_s   <- matrix(NA, nrow = s_hf, ncol = n_runs)
  STDEVs_s <- matrix(NA, nrow = n_spatial, ncol = n_runs)
  COVERs_s <- matrix(NA, nrow = s_hf, ncol = n_runs)
  
  # Temporal metrics (matrices: rows = time, cols = runs)
  RMSEs_m  <- matrix(NA, nrow = n_month, ncol = n_runs)
  MAEs_m   <- matrix(NA, nrow = n_month, ncol = n_runs)
  MBEs_m   <- matrix(NA, nrow = n_month, ncol = n_runs)
  STDEVs_m <- matrix(NA, nrow = n_month, ncol = n_runs)
  COVERs_m <- matrix(NA, nrow = n_month, ncol = n_runs)
  
  RMSEs_y  <- matrix(NA, nrow = n_year, ncol = n_runs)
  MAEs_y   <- matrix(NA, nrow = n_year, ncol = n_runs)
  MBEs_y   <- matrix(NA, nrow = n_year, ncol = n_runs)
  STDEVs_y <- matrix(NA, nrow = n_year, ncol = n_runs)
  COVERs_y <- matrix(NA, nrow = n_year, ncol = n_runs)
  
  pb <- progress_bar$new(
    format = "  Getting CV stats [:bar] :percent ETA: :eta",
    total = length(files),
    clear = TRUE,
    width = 70
  )
  
  for (i in seq_along(files)) {
    load(files[i])
    
    if (fidelity == "lf") {
      cv_stats <- cv_stats_lf
      ranks <- ranks_lf
    } else if (fidelity == "hf") {
      cv_stats <- cv_stats_hf
      ranks <- ranks_hf
    } else {
      cv_stats <- cv_stats_mf
      ranks <- ranks_lf
      ranks_discrep <- ranks_discrep
    }
    
    ranks_mat[, i] <- ranks
    
    if (fidelity == "mf") {
      ranks_discrep_mat[, i] <- ranks_discrep
    }
    
    RMSEs[i]   <- cv_stats$overall$RMSE
    MAEs[i]    <- cv_stats$overall$MAE
    MBEs[i]    <- cv_stats$overall$MBE
    STDEVs[i]  <- cv_stats$overall$STDEV
    COVERs[i]  <- cv_stats$overall$COVER
    
    RMSEs_s[, i]  <- cv_stats$spatial$RMSE
    MAEs_s[, i]   <- cv_stats$spatial$MAE
    MBEs_s[, i]   <- cv_stats$spatial$MBE
    STDEVs_s[, i] <- cv_stats$spatial$STDEV
    COVERs_s[, i] <- cv_stats$spatial$COVER
    
    RMSEs_m[, i]  <- cv_stats$month$RMSE
    MAEs_m[, i]   <- cv_stats$month$MAE
    MBEs_m[, i]   <- cv_stats$month$MBE
    STDEVs_m[, i] <- cv_stats$month$STDEV
    COVERs_m[, i] <- cv_stats$month$COVER
    
    RMSEs_y[, i]  <- cv_stats$year$RMSE
    MAEs_y[, i]   <- cv_stats$year$MAE
    MBEs_y[, i]   <- cv_stats$year$MBE
    STDEVs_y[, i] <- cv_stats$year$STDEV
    COVERs_y[, i] <- cv_stats$year$COVER
    
    gc()
    
    pb$tick()
  }
  
  return(list(
    run_ids = run_ids,
    
    ranks_mat = ranks_mat,
    ranks_discrep_mat = ranks_discrep_mat,
    
    # Overall
    RMSEs   = RMSEs,
    MAEs    = MAEs,
    MBEs    = MBEs,
    STDEVs  = STDEVs,
    COVERs  = COVERs,
    
    # Spatial (matrices)
    RMSEs_s  = RMSEs_s,
    MAEs_s   = MAEs_s,
    MBEs_s   = MBEs_s,
    STDEVs_s = STDEVs_s,
    COVERs_s = COVERs_s,
    
    # Temporal (month)
    RMSEs_m  = RMSEs_m,
    MAEs_m   = MAEs_m,
    MBEs_m   = MBEs_m,
    STDEVs_m = STDEVs_m,
    COVERs_m = COVERs_m,
    
    # Temporal (year)
    RMSEs_y  = RMSEs_y,
    MAEs_y   = MAEs_y,
    MBEs_y   = MBEs_y,
    STDEVs_y = STDEVs_y,
    COVERs_y = COVERs_y
  ))
}

get_all_cv_stats <- function(folder_dir, fidelity) {
  files <- list.files(folder_dir, full.names = TRUE)
  
  # Extract the numeric part (before ".Rdata")
  file_nums <- as.numeric(gsub(".*/|\\.Rdata", "", files))
  
  # Order by number
  files <- files[order(file_nums)]
  
  print(files)
  n_runs <- length(files)
  
  load(files[1])  # loads post_preds
  if (fidelity == "lf") {
    post_preds <- post_preds_lf_noise
    ranks <- ranks_lf
    dims <- dim(post_preds$mse)
  } else if (fidelity == "hf") {
    post_preds <- post_preds_hf_noise
    ranks <- ranks_hf
    dims <- dim(post_preds$mse)
  } else {
    post_preds <- post_preds_mf_noise
    ranks <- ranks_lf
    ranks_discrep <- ranks_discrep
    dims <- dim(post_preds$mse)
  }
  n_spatial <- length(post_preds$spatial$sd)
  n_month <- length(post_preds$month$mse)
  n_year <- length(post_preds$year$mse)
  
  
  # Ranks
  ranks_mat <- matrix(NA, nrow = length(ranks), ncol = n_runs)
  
  if (fidelity == "mf") {
    ranks_discrep_mat <- matrix(NA, nrow = length(ranks_discrep), ncol = n_runs)
  } else {
    ranks_discrep_mat <- NULL
  }
  
  # Overall metrics (vectors)
  MSEs   <- array(NA, dim = c(dims, n_runs))
  MAEs    <- array(NA, dim = c(dims, n_runs))
  MBEs    <- array(NA, dim = c(dims, n_runs))
  STDEVs  <- array(NA, dim = c(dims, n_runs))
  COVERs  <- array(NA, dim = c(dims, n_runs))
  run_ids <- character(n_runs)
  
  
  pb <- progress_bar$new(
    format = "  Getting CV stats [:bar] :percent ETA: :eta",
    total = length(files),
    clear = TRUE,
    width = 70
  )
  
  for (i in seq_along(files)) {
    load(files[i])
    
    if (fidelity == "lf") {
      # cv_stats <- cv_stats_lf
      post_preds_noise <- post_preds_lf_noise
      ranks <- ranks_lf
    } else if (fidelity == "hf") {
      # cv_stats <- cv_stats_hf
      post_preds_noise <- post_preds_hf_noise
      ranks <- ranks_hf
    } else {
      # cv_stats <- cv_stats_mf
      post_preds_noise <- post_preds_mf_noise
      ranks <- ranks_lf
      ranks_discrep <- ranks_discrep
    }
    
    ranks_mat[, i] <- ranks
    
    if (fidelity == "mf") {
      ranks_discrep_mat[, i] <- ranks_discrep
    }
    
    MSEs[, , , i]   <- post_preds_noise$mse
    MAEs[, , , i]   <- post_preds_noise$mae
    MBEs[, , , i]   <- post_preds_noise$mbe
    STDEVs[, , , i] <- post_preds_noise$sd
    COVERs[, , , i] <- post_preds_noise$cover
    
    gc()
    
    pb$tick()
  }
  
  return(list(
    run_ids = run_ids,
    
    ranks_mat = ranks_mat,
    ranks_discrep_mat = ranks_discrep_mat,
    
    MSEs   = MSEs,
    MAEs    = MAEs,
    MBEs    = MBEs,
    STDEVs  = STDEVs,
    COVERs  = COVERs
  ))
}

