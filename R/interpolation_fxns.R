fast_knn_interpolate_first_dim <- function(data_array,
                                           coords_old,
                                           coords_new,
                                           k = 5,
                                           epsilon = 1e-8) {
  # Dimensions
  dims <- dim(data_array)
  n_new <- nrow(coords_new)
  n_old <- nrow(coords_old)
  other_dims <- prod(dims[-1])

  # Flatten array to [space, ...] → [n_old, N]
  data_mat <- matrix(data_array, nrow = n_old, ncol = other_dims)

  # Get k nearest neighbors and distances
  knn_res <- FNN::get.knnx(data = coords_old, query = coords_new, k = k)

  nn_idx <- knn_res$nn.index # [n_new, k]
  nn_dist <- knn_res$nn.dist + epsilon # avoid div by 0

  # Inverse distance weights
  weights <- 1 / nn_dist
  weights <- weights / rowSums(weights) # normalize [n_new, k]

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
