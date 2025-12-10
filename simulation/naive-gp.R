setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load functions
source("../fxns.R")

# Load previously saved setup data
load("save/setup.Rdata")

test_inputs <- matrix(c(-1.7, -1.7, -1.7), nrow = 1) # can have multiple inputs stored in rows

true_outputs <- simulate_for_design(test_inputs, 
                                    coord_lf, coord_hf, 
                                    year_ids, month_ids,
                                    start_lf = 1, end_lf = 5)
# Arrays to store results
n_space <- nrow(coord_hf)
n_months <- 12
n_years <- dim(output_hf)[3]
n_test <- nrow(test_inputs)

mean_naive_gp <- array(NA, dim = c(n_space, n_months, n_years, n_test))
mse_naive_gp <- array(NA, dim = c(n_space, n_months, n_years, n_test))
sd_naive_gp <- array(NA, dim = c(n_space, n_months, n_years, n_test))
cover_naive_gp <- array(NA, dim = c(n_space, n_months, n_years, n_test))

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta",
  total = n_space * n_months * n_years, clear = FALSE
)

for (s in 1:n_space) {
  for (m in 1:n_months) {
    for (y in 1:n_years) {
      # -----------------------
      # Fit GP silently
      # -----------------------
      gpr <- km(design = design[runs_hf, ], response = output_hf[s, m, y, ], 
                nugget.estim = F, 
                control = list(trace = F),
                lower = c(0.5, 0.5),
                upper = c())
      
      preds <- predict(gpr, newdata  = test_inputs, type = "SK")
      
      # -----------------------
      # Posterior predictive draws
      # -----------------------
      mean_naive_gp[s, m, y, ] <- preds$mean
      
      # -----------------------
      # Squared errors per test input
      # -----------------------
      mse_naive_gp[s, m, y, ] <- (preds$mean - true_outputs$high_fi[s, m, y, ])^2
      
      # -----------------------
      # Posterior SD
      # -----------------------
      sd_naive_gp[s, m, y, ] <- preds$sd
      
      # -----------------------
      # Coverage per test input (fraction of draws containing true value)
      # -----------------------
      cover_naive_gp[s, m, y, ] <- (true_outputs$high_fi[s, m, y, ] >= preds$lower95) &
        (true_outputs$high_fi[s, m, y, ] <= preds$upper95)
      # -----------------------
      # Update progress bar
      # -----------------------
      pb$tick()
    }
  }
}

save(
  mean_naive_gp,
  mse_naive_gp,
  sd_naive_gp,
  cover_naive_gp,
  test_inputs,
  file = "save/naive_gp.Rdata"
)

rm(list = ls())
gc()
