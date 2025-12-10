# `mf_tensor` Overview

Functions in `fxns.R` are designed for applications in both **simulation studies** and the **MPAS-Seaice** project.

---

## Workflow

### Simulation Study

The standard workflow involves three main steps: data setup, emulation, and result saving.

1.  **Data Preparation:** Run `simulation/setup.R` to create the simulation study data.
    * *Output:* Data is saved in the subfolder `simulation/save/`.

2.  **Run Emulation Functions:** Execute the scripts for the two different emulation methods:
    * **a. Tensor-based Emulation:**
        * `simulation/lf-tensor.R`
        * `simulation/hf-tensor.R`
        * `simulation/mf-tensor.R`
    * **b. Naïve-GP Emulation:**
        * `simulation/naive-gp.R`

3.  **Analysis:** Emulation results are automatically saved in `simulation/save/` for further analysis.

---

### Custom Data Replacement

You can replace the standard simulation data with your own by defining the following variables:

| Variable | Description |
| :--- | :--- |
| `output_lf`, `output_hf` | Arrays containing the **Low-Fidelity (LF)** and **High-Fidelity (HF)** outputs. |
| `design` | The full design matrix for all runs. |
| `runs_lf`, `runs_hf` | Row indices specifying which runs belong to the LF and HF fidelities. |
| `d_lf`, `d_hf` | The number of runs for each fidelity. |