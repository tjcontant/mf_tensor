# `mf_tensor`

R code for multi-fidelity Gaussian process emulation with tensor
decompositions, accompanying [arXiv:2603.04697](https://arxiv.org/abs/2603.04697).
Supplementary material is in `supp.pdf`.

## Requirements

R >= 4.1 on macOS or Linux (MCMC chains run in forked processes), with:

```r
install.packages(c(
  "rTensor", "Matrix", "mvtnorm", "matrixStats", "truncnorm", "progress",
  "pbapply", "pbmcapply", "FNN", "fields", "lhs", "concaveman"
))
```

## Running the simulation study

From the repository root:

```r
source("Simulation/R/run_all.R")
```

This runs, in order:

1. `generate_data.R` and `generate_test_inputs.R`: simulate the LF/HF data and
   held-out test inputs.
2. `lf_tensor.R`, `hf_tensor.R`, `mf_tensor.R`: fit the LF-only, HF-only and
   multi-fidelity tensor GPs.
3. `compare_marginal_metrics.R` and `compare_variogram_score.R`: compare
   methods by absolute error, coverage, CRPS and variogram score.

**Note:** Some of the simulation setup is relaxed relative to the paper, with
fewer LF runs and MCMC iterations, so the scripts run quickly.

All scripts are in `Simulation/R/`, and results are saved to `Simulation/Data/`.
Function definitions are in `R/`.
