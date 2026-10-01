# Setup -------------------------------------------------------------------

set.seed(4382)

rm(list = ls())
gc()

source("R/helper_fxns.R")
source("R/interpolation_fxns.R")

# Shubert function --------------------------------------------------------

shubert <- Vectorize(function(s1, s2, m, y, d1, d2, start, end) {
  ii <- start:end

  alpha <- 1 +  (1 / 20) * (sin(s1 * (d1 / 2 + 0.5)) + cos(s2 * (d2 / 2 + 0.5))) * sin(m) +
    sin(y) / 2
  beta <- (1 / 10) * (s1 + s2) * cos(m) + (1 / 20) * sin(m) + (1 / 5) * y

  sum1 <- sum(ii * cos((ii + 1) * d1 * alpha + ii + beta))
  sum2 <- sum(ii * cos((ii + 1) * d2 * alpha + ii + beta))

  y <- sum1 + sum2 + (1 / 5) * y
  return(y)
})

# Simulator ---------------------------------------------------------------

simulate_for_design <- function(
  design,
  coord_lf, coord_hf,
  year_ids,
  month_ids,
  start_lf, end_lf,
  k_interp = 4
) {
  n_runs <- nrow(design)
  n_months <- length(unique(month_ids))
  n_years <- length(unique(year_ids))
  n_time <- n_months * n_years

  n_lf <- nrow(coord_lf)
  n_hf <- nrow(coord_hf)

  # Allocate output arrays: [space, month, year, design]
  M_lf_all <- array(0, dim = c(n_lf, n_months, n_years, n_runs))
  M_interp_all <- array(0, dim = c(n_hf, n_months, n_years, n_runs))
  M_hf_all <- array(0, dim = c(n_hf, n_months, n_years, n_runs))
  M_disc_all <- array(0, dim = c(n_hf, n_months, n_years, n_runs))

  pb <- progress::progress_bar$new(
    format = "  generating simulations [:bar] :percent ETA: :eta",
    total = n_runs, clear = TRUE
  )

  for (i in 1:n_runs) {
    set.seed(i)

    u <- design[i, 1]
    v <- design[i, 2]

    a <- u + v
    b <- v - u

    M_lf <- array(0, dim = c(n_lf, n_months, n_years))
    M_hf <- array(0, dim = c(n_hf, n_months, n_years))
    M_disc <- array(0, dim = c(n_hf, n_months, n_years))
    M_interp <- array(0, dim = c(n_hf, n_months, n_years))

    for (t in 1:n_time) {
      yr <- year_ids[t]
      mo <- month_ids[t]
      mo_cyclic <- 2 * pi * (mo - 5) / 12

      # LF field
      M_lf[, mo, yr] <- shubert(coord_lf$x, coord_lf$y, mo_cyclic, yr, a, b,
        start = start_lf, end = end_lf
      ) |>
        scale_values() + (1 / 5) * yr

      M_lf[, mo, yr] <- M_lf[, mo, yr] + rnorm(n_lf, 0, 0.01)

      # Interpolate LF to HF grid
      M_interp[, mo, yr] <- fast_knn_interpolate_first_dim(M_lf[, mo, yr],
        coord_lf, coord_hf,
        k = k_interp
      )

      # Discrepancy
      M_disc[, mo, yr] <- mvtnorm::dmvnorm(coord_hf,
        mean = c(0, 0),
        sigma = diag(2, 2)
      ) |>
        scale_values() * (1 + abs((cos(a * mo_cyclic) + yr * b^2)))

      M_disc[, mo, yr] <- M_disc[, mo, yr] + rnorm(n_hf, 0, 0.05)

      # HF field
      M_hf[, mo, yr] <- M_disc[, mo, yr] + M_interp[, mo, yr]
    }

    M_lf_all[, , , i] <- M_lf
    M_interp_all[, , , i] <- M_interp
    M_hf_all[, , , i] <- M_hf
    M_disc_all[, , , i] <- M_disc

    pb$tick()
  }

  return(list(
    low_fi = M_lf_all,
    interpolated_low_fi = M_interp_all,
    high_fi = M_hf_all,
    discrepancy = M_disc_all
  ))
}

# Simulator specifications ------------------------------------------------

ns_1_lf <- 25
ns_2_lf <- 25
s_x_lf <- seq(-4, 4, length.out = ns_1_lf)
s_y_lf <- seq(-4, 4, length.out = ns_2_lf)
coord_lf <- expand.grid(x = s_x_lf, y = s_y_lf)
s_lf <- nrow(coord_lf)

ns_1_hf <- 50
ns_2_hf <- 50
s_x_hf <- seq(-4, 4, length.out = ns_1_hf)
s_y_hf <- seq(-4, 4, length.out = ns_2_hf)
coord_hf <- expand.grid(x = s_x_hf, y = s_y_hf)
s_hf <- nrow(coord_hf)

n_years <- 5
n_months <- 12
year_ids <- rep(1:n_years, each = n_months)
month_ids <- rep(1:n_months, times = n_years)
t_lf <- t_hf <- n_years * n_months

n_lf <- s_lf * t_lf
n_hf <- s_hf * t_hf

nx_lf <- 50
nx_hf <- 10
np <- 3

runs_lf <- 1:nx_lf

design <- lhs::maximinLHS(nx_lf, np)

# Generate simulations ----------------------------------------------------

outputs <- simulate_for_design(
  design, coord_lf, coord_hf,
  year_ids,
  month_ids,
  start_lf = 1,
  end_lf = 5,
  k_interp = 4
)

# Choose HF runs ----------------------------------------------------------

kmeans_result <- kmeans(design[runs_lf, ], centers = nx_hf)

runs_hf <- sapply(1:nx_hf, function(k) {
  idx <- which(kmeans_result$cluster == k)
  center <- kmeans_result$centers[k, ]
  dists <- rowSums((design[runs_lf[idx], ] - matrix(center, nrow = length(idx), ncol = ncol(design), byrow = TRUE))^2)
  runs_lf[idx[which.min(dists)]]
})

runs_hf <- sort(runs_hf)

# Save results ------------------------------------------------------------

output_lf <- outputs$low_fi[, , , runs_lf]
output_lf_interp <- outputs$interpolated_low_fi[, , , runs_hf]
discrepancy <- outputs$discrepancy[, , , runs_hf]
output_hf <- outputs$high_fi[, , , runs_hf] # Confirmed: interpolated LF + discrepancy

dir.create("Simulation/Data", recursive = TRUE, showWarnings = FALSE)
save.image(paste0("Simulation/Data/simulation_setup_", get_date(), ".Rdata"))

# Check results -----------------------------------------------------------

# HF should equal interpolated LF + discrepancy
stopifnot(max(abs(output_hf - (output_lf_interp + discrepancy))) < 1e-10)

# Clean up ----------------------------------------------------------------

rm(list = ls())
gc()
