# Scoring rules: Gaussian CRPS and the (subsampled, weighted) Variogram Score

# Closed-form Gaussian CRPS (Gneiting & Raftery 2007)
crps_gaussian <- function(y, mu, sigma) {
  z <- (y - mu) / sigma
  sigma * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}

# Variogram Score (Scheuerer & Hamill 2015) for ONE spatial field
variogram_score <- function(draws, truth, p = 0.5) {
  n_loc <- length(truth)
  n_draws <- nrow(draws)

  # Observed pairwise dissimilarity matrix: |y_i - y_j|^p for every location
  # pair, computed once from the single true field.
  true_diff <- abs(outer(truth, truth, "-"))^p

  # Monte Carlo estimate of E_hat[|X_i - X_j|^p]: accumulate the same pairwise
  # dissimilarity matrix for each posterior/ensemble draw, then average
  exp_diff <- matrix(0, n_loc, n_loc)
  for (k in seq_len(n_draws)) {
    exp_diff <- exp_diff + abs(outer(draws[k, ], draws[k, ], "-"))^p
  }
  exp_diff <- exp_diff / n_draws

  # Sum of squared errors between observed and predicted dissimilarity, over
  # unique location pairs only (upper triangle -- the matrix is symmetric with a
  # zero diagonal, so the lower triangle and diagonal add no information and
  # would only double-count)
  ut <- upper.tri(true_diff)
  sum((true_diff[ut] - exp_diff[ut])^2)
}

# Variogram Score estimated from a sample of location pairs
variogram_score_subsampled <- function(draws, truth, pair_i, pair_j, p = 0.5, weights = NULL) {
  n_loc <- length(truth)
  n_draws <- nrow(draws)
  n_pairs <- length(pair_i)

  if (is.null(weights)) {
    weights <- rep(1, n_pairs)
  } else {
    weights <- weights / mean(weights)
  }

  true_diff <- abs(truth[pair_i] - truth[pair_j])^p

  exp_diff <- numeric(n_pairs)
  for (k in seq_len(n_draws)) {
    exp_diff <- exp_diff + abs(draws[k, pair_i] - draws[k, pair_j])^p
  }
  exp_diff <- exp_diff / n_draws

  total_pairs <- n_loc * (n_loc - 1) / 2
  (total_pairs / n_pairs) * sum(weights * (true_diff - exp_diff)^2)
}

# Samples location-index pairs within a given spatial distance
# and/or temporal lag
sample_pairs_restricted <- function(coord, n_time, n_pairs, max_spatial_dist = Inf, max_temporal_lag = Inf, keep_loc = NULL) {
  n_space <- nrow(coord)

  Ds <- fields::rdist(coord)
  Dt <- abs(outer(seq_len(n_time), seq_len(n_time), "-"))

  # Precomputed ONCE per possible space/time index (n_space + n_time lookups
  # total), not per sampled pair -- the candidate pool for a given index never
  # changes
  space_candidates <- lapply(seq_len(n_space), function(s) which(Ds[s, ] <= max_spatial_dist))
  time_candidates <- lapply(seq_len(n_time), function(t) which(Dt[t, ] <= max_temporal_lag))

  draw_space_j <- function(s_vec) {
    vapply(s_vec, function(s) {
      cand <- space_candidates[[s]]
      cand[sample.int(length(cand), 1)]
    }, integer(1))
  }
  draw_time_j <- function(t_vec) {
    vapply(t_vec, function(t) {
      cand <- time_candidates[[t]]
      cand[sample.int(length(cand), 1)]
    }, integer(1))
  }

  space_i <- sample.int(n_space, n_pairs, replace = TRUE)
  time_i <- sample.int(n_time, n_pairs, replace = TRUE)
  i <- space_i + n_space * (time_i - 1)

  # space_i/time_i are drawn independently, so a keep_loc rejection on the
  # COMBINED index i must redraw both together, not just one -- unlike j below,
  # whose candidate pool is already conditioned on (the now-fixed) i
  if (!is.null(keep_loc)) {
    bad_i <- !keep_loc[i]
    while (any(bad_i)) {
      space_i[bad_i] <- sample.int(n_space, sum(bad_i), replace = TRUE)
      time_i[bad_i] <- sample.int(n_time, sum(bad_i), replace = TRUE)
      i[bad_i] <- space_i[bad_i] + n_space * (time_i[bad_i] - 1)
      bad_i <- !keep_loc[i]
    }
  }

  space_j <- draw_space_j(space_i)
  time_j <- draw_time_j(time_i)
  j <- space_j + n_space * (time_j - 1)

  bad <- i == j
  if (!is.null(keep_loc)) bad <- bad | !keep_loc[j]
  while (any(bad)) {
    space_j[bad] <- draw_space_j(space_i[bad])
    time_j[bad] <- draw_time_j(time_i[bad])
    j[bad] <- space_j[bad] + n_space * (time_j[bad] - 1)
    bad <- i == j
    if (!is.null(keep_loc)) bad <- bad | !keep_loc[j]
  }

  list(i = i, j = j)
}

# Per-pair weights h_ij = 1/dist(i,j) for variogram_score_subsampled(),
# following Scheuerer & Hamill
spatiotemporal_pair_weights <- function(coord, n_time, pair_i, pair_j, spatial_scale, temporal_scale) {
  n_space <- nrow(coord)

  space_i <- ((pair_i - 1) %% n_space) + 1
  space_j <- ((pair_j - 1) %% n_space) + 1
  time_i <- ((pair_i - 1) %/% n_space) + 1
  time_j <- ((pair_j - 1) %/% n_space) + 1

  spatial_dist <- sqrt(rowSums((coord[space_i, , drop = FALSE] - coord[space_j, , drop = FALSE])^2))
  temporal_dist <- abs(time_i - time_j)

  dist <- sqrt((spatial_dist / spatial_scale)^2 + (temporal_dist / temporal_scale)^2)
  # pair_i != pair_j, so dist > 0
  1 / pmax(dist, .Machine$double.eps)
}

