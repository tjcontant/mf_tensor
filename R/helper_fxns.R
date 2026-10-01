get_date <- function() {
  format(Sys.Date(), "%d%m%y")
}

# Loads the most-recently-modified file matching `pattern` in `dir`
load_latest <- function(dir, pattern, envir = .GlobalEnv) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matching '", pattern, "' found in ", dir)
  latest <- files[which.max(file.info(files)$mtime)]
  load(latest, envir = envir)
}

scale_values <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# TRUE for rows of `query` inside the concave hull of `pts` (first two columns)
in_shrinkwrap <- function(query, pts) {
  query <- as.matrix(query)
  poly <- concaveman::concaveman(as.matrix(pts))
  n_edge <- nrow(poly) - 1
  inside <- rep(FALSE, nrow(query))
  for (i in seq_len(n_edge)) {
    p1 <- poly[i, ]
    p2 <- poly[i + 1, ]
    crosses <- ((p1[2] > query[, 2]) != (p2[2] > query[, 2])) &
      (query[, 1] < (p2[1] - p1[1]) * (query[, 2] - p1[2]) / (p2[2] - p1[2]) + p1[1])
    inside <- xor(inside, crosses)
  }
  inside
}

# Draw n candidate points uniformly in [lower, upper]^ncol(X_train), rejecting
# any candidate within min_dist of an existing training point, so the result is
# a genuinely held-out test set
generate_test_inputs <- function(n, X_train, lower = -2, upper = 0, min_dist = 0.05, weights = NULL) {
  n_col <- ncol(X_train)
  test_inputs <- matrix(NA, nrow = n, ncol = n_col)
  count <- 0

  if (is.null(weights)) weights <- rep(1, n_col)

  while (count < n) {
    candidate <- runif(n_col, lower, upper)
    dists <- sqrt(rowSums((t((t(X_train) - candidate) * weights))^2))
    if (all(dists >= min_dist)) {
      count <- count + 1
      X_train <- rbind(X_train, candidate)
      test_inputs[count, ] <- candidate
    }
  }
  return(test_inputs)
}
