# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

# Clear working environment
rm(list = ls())

# Set working directory to curretn location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Base Shubert function — frequency controlled via n_terms
shubert <- Vectorize(function(s1, s2, m, y, d1, d2, start, end) {
  ii <- start:end
  
  alpha <- 1 + (1/20) * (sin(s1) + cos(s2)) * sin(m)
  beta  <- (1/20) * (s1 + s2) + (1/40) * sin(m) + (1/5) * y
  
  sum1 <- sum(ii * cos((ii + 1) * d1 * alpha + ii + beta))
  sum2 <- sum(ii * cos((ii + 1) * d2 * alpha + ii + beta))
  
  y <- sum1 * sum2
  return(y)
})

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}


# ------------------------------------------------------------------------------
# INTERPOLATION FUNCTION
# ------------------------------------------------------------------------------
interpolate_locations <- function(coords_from, outputs, coords_to, k) {
  knn_result <- knn.reg(
    train = coords_from,
    test  = coords_to,
    y     = outputs,
    k     = k
  )
  return(knn_result$pred)
}


# ------------------------------------------------------------------------------
# UNIFIED MULTI-DESIGN SIMULATION FUNCTION
# ------------------------------------------------------------------------------
simulate_for_design <- function(design, coord_lf, coord_hf,
                                year_ids, month_ids,
                                start_lf, end_lf,
                                k_interp = 10) {
  
  n_runs   <- nrow(design)
  n_months <- length(unique(month_ids))
  n_years  <- length(unique(year_ids))
  n_time   <- n_months * n_years
  
  n_lf <- nrow(coord_lf)
  n_hf <- nrow(coord_hf)
  
  # Allocate output arrays: [space, month, year, design]
  M_lf_all      <- array(0, dim = c(n_lf, n_months, n_years, n_runs))
  M_interp_all  <- array(0, dim = c(n_hf, n_months, n_years, n_runs))
  M_hf_all      <- array(0, dim = c(n_hf, n_months, n_years, n_runs))
  M_disc_all    <- array(0, dim = c(n_hf, n_months, n_years, n_runs))
  
  pb <- progress::progress_bar$new(
    format = "  generating simulations [:bar] :percent ETA: :eta",
    total = n_runs, clear = TRUE
  )
  
  for (i in 1:n_runs) {
    a <- design[i, 1]
    b <- design[i, 2]
    
    M_lf  <- array(0, dim = c(n_lf, n_months, n_years))
    M_hf  <- array(0, dim = c(n_hf, n_months, n_years))
    M_disc <- array(0, dim = c(n_hf, n_months, n_years))
    M_interp <- array(0, dim = c(n_hf, n_months, n_years))
    
    for (t in 1:n_time) {
      yr <- year_ids[t]
      mo <- month_ids[t]
      mo_cyclic <- 2 * pi * (mo - 5) / 12
      
      # LF field
      M_lf[, mo, yr] <- shubert(coord_lf$x, coord_lf$y, mo_cyclic, yr, a, b, 
                                start = start_lf, end = end_lf) |>
        scale_values() + (1/5) * yr
      
      # Interpolate LF to HF grid
      M_interp[, mo, yr] <- interpolate_locations(coord_lf, M_lf[, mo, yr], coord_hf, k = k_interp)
      
      # Discrepancy
      sigma_val <- 1 + 0.1 * sin(mo_cyclic) + 0.1 * yr
      Sigma <- sigma_val * diag(2)
      
      M_disc[, mo, yr] <- dmvnorm(coord_hf, mean = 0.5 * c(a, b), sigma = Sigma) |>
        scale_values() * 0.25 + 1
      
      # HF field
      M_hf[, mo, yr] <- M_disc[, mo, yr] + M_interp[, mo, yr]
      
    }
    
    M_lf_all[, , , i]     <- M_lf
    M_interp_all[, , , i] <- M_interp
    M_hf_all[, , , i]     <- M_hf
    M_disc_all[, , , i]   <- M_disc
    
    pb$tick()
  }
  
  return(list(
    low_fi = M_lf_all,
    interpolated_low_fi = M_interp_all,
    high_fi = M_hf_all,
    discrepancy = M_disc_all
  ))
}


# ------------------------------------------------------------------------------
# EXAMPLE USAGE
# ------------------------------------------------------------------------------

set.seed(432)

# Spatial grids
s_x_lf <- seq(-4, 4, length.out = 25)
s_y_lf <- seq(-4, 4, length.out = 25)
coord_lf <- expand.grid(x = s_x_lf, y = s_y_lf)
s_lf <- nrow(coord_lf)

s_x_hf <- seq(-4, 4, length.out = 50)
s_y_hf <- seq(-4, 4, length.out = 50)
coord_hf <- expand.grid(x = s_x_hf, y = s_y_hf)
s_hf <- nrow(coord_hf)


# Temporal structure
n_years <- 5
n_months <- 12
year_ids <- rep(1:n_years, each = n_months)
month_ids <- rep(1:n_months, times = n_years)
t_lf <- t_hf <- n_years * n_months

n_lf <- s_lf * t_lf
n_hf <- s_hf * t_hf

# Design matrix
design <- lhs::randomLHS(100, 3)

# Run simulation for all designs
outputs <- simulate_for_design(design, coord_lf, coord_hf, year_ids, month_ids,
                               start_lf = 1, end_lf = 5)


# ------------------------------------------------------------------------------
# SELECT LF & HF RUNS
# ------------------------------------------------------------------------------
d_lf <- nrow(design)
d_hf <- round(d_lf / 10)

runs_lf <- 1:d_lf

# Use k-means to spread HF runs across design space
set.seed(428)
kmeans_result <- kmeans(design[runs_lf, ], centers = d_hf)

runs_hf <- sapply(1:d_hf, function(k) {
  idx <- which(kmeans_result$cluster == k)
  center <- kmeans_result$centers[k, ]
  dists <- rowSums((design[runs_lf[idx], ] - matrix(center, nrow = length(idx), ncol = ncol(design), byrow = TRUE))^2)
  runs_lf[idx[which.min(dists)]]
})
runs_hf <- sort(runs_hf)

# Check design selection visually
dev.off()
plot(design[runs_lf, c(1, 3)], main = "Design: LF (black) and HF (red)")
points(design[runs_hf, c(1, 3)], col = "red", pch = 16)
text(design[runs_hf, ], labels = runs_hf, pos = 3, col = "red", cex = 0.8)


# ------------------------------------------------------------------------------
# FINAL OUTPUTS
# ------------------------------------------------------------------------------

output_lf <- outputs$low_fi[, , , runs_lf]
output_lf_interp <- outputs$interpolated_low_fi[, , , runs_hf]
discrepancy <- outputs$discrepancy[, , , runs_hf]
output_hf <- outputs$high_fi[, , , runs_hf]  # Confirmed: interpolated LF + discrepancy

# rm(interpolate_locations)

# Save
save.image("save/setup.Rdata")


# ------------------------------------------------------------------------------
# CHECK: HF = LF_interp + Discrepancy
# ------------------------------------------------------------------------------

# Compute the difference
reconstructed_hf <- output_lf_interp + discrepancy
diff_array <- abs(output_hf - reconstructed_hf)

# Maximum absolute error
max_error <- max(diff_array)
mean_error <- mean(diff_array)

cat(sprintf("\nMax difference between HF and LF_interp + discrepancy: %.6e\n", max_error))
cat(sprintf("Mean difference: %.6e\n", mean_error))

# Check if within numerical tolerance
tolerance <- 1e-10
if (max_error < tolerance) {
  cat("\t✅ HF is equal to LF_interp + discrepancy within tolerance.\n")
} else {
  cat("\tWARNING: HF does not match LF_interp + discrepancy. Check simulation code.\n")
}

rm(list = ls())
