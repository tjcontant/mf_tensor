# Tucker decomposition machinery used by the tensor-GP emulators
# (R/tensor_gp_fxns.R)

library(rTensor)
library(Matrix)
library(progress)

safe_ttl <- function(tnsr, U_list, ms) {
  # Helper to check if a factor is a diagonal matrix
  is.diagonal <- function(mat) {
    inherits(mat, "ddiMatrix") # Diagonal class from Matrix package
  }

  # Keep only factors that are not diagonal
  non_diag_idx <- which(!sapply(U_list, is.diagonal))

  if (length(non_diag_idx) == 0) {
    return(tnsr)
  }

  # Only use corresponding mode indices
  ttl(tnsr, U_list[non_diag_idx], ms = ms[non_diag_idx])
}

# updated for pca
my_tucker <- function(tnsr, ranks = NULL, max_iter = 25, tol = 1e-05) {
  stopifnot(is(tnsr, "Tensor"))
  if (is.null(ranks)) {
    stop("ranks must be specified")
  }
  if (sum(ranks > tnsr@modes) != 0) {
    stop("ranks must be smaller than the corresponding mode")
  }
  if (sum(ranks <= 0) != 0) {
    stop("ranks must be positive")
  }
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    temp_mat <- rs_unfold(tnsr, m = m)@data
    if (ranks[m] == tnsr@modes[m]) {
      U_list[[m]] <- Diagonal(ranks[m]) # sparse
    } else {
      U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
    }
  }

  tnsr_norm <- rTensor::fnorm(tnsr)
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z, U_list) {
    est <- safe_ttl(Z, U_list, ms = 1:num_modes)
    curr_resid <- rTensor::fnorm(tnsr - est)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1) {
      return(FALSE)
    }
    if (abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm <
      tol) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  while ((curr_iter < max_iter) && (!converged)) {
    setTxtProgressBar(pb, curr_iter)
    modes <- tnsr@modes
    modes_seq <- 1:num_modes
    for (m in modes_seq) {
      X <- safe_ttl(tnsr, lapply(U_list[-m], t), ms = modes_seq[-m])
      if (ranks[m] < tnsr@modes[m]) {
        U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
      }
    }
    Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)
    if (CHECK_CONV(Z, U_list)) {
      converged <- TRUE
      setTxtProgressBar(pb, max_iter)
    } else {
      curr_iter <- curr_iter + 1
    }
  }
  close(pb)
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) *
    100
  est <- safe_ttl(Z, U_list, ms = 1:num_modes)
  invisible(list(
    Z = Z, U = U_list, conv = converged, est = est,
    norm_percent = norm_percent, fnorm_resid = tail(
      fnorm_resid,
      1
    ), all_resids = fnorm_resid
  ))
}

# Applies Tucker decomposition to a 3D output array
apply_tucker_decomp <- function(output, ranks, method, calculate_bases) {
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
    tucker_decomp <- my_tucker(output_tens, ranks) |> suppressMessages()
  }

  G <- tucker_decomp$Z # core tensor
  U_list <- tucker_decomp$U

  # Reconstruct the tensor from Tucker components
  output_tens_reduced <- safe_ttl(G, U_list, 1:1:n_modes)
  output_array_reduced <- output_tens_reduced@data

  # Calculate the explained variance
  norm_original <- rTensor::fnorm(output_tens)
  norm_approx <- rTensor::fnorm(output_tens_reduced)
  norm_diff <- rTensor::fnorm(output_tens - output_tens_reduced)

  explained_var_tucker <- 100 * (1 - norm_diff^2 / norm_original^2)

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
  Reduce(
    function(x, y) kronecker(y, x), # Reverse order here
    Map(function(U, i) U[, i], U_list[-length(U_list)], idx)
  )
}
