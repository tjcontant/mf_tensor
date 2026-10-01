library(progress)
library(rTensor)

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

