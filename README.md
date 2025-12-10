# `mf_tensor` Overview

This repository provides an R implementation for **Multi-Fidelity (MF) Gaussian Process (GP) emulation** utilizing **tensor-product covariance structures**. The core functions in `fxns.R` are designed to support two main application areas: **controlled simulation studies** and the **MPAS-Seaice modeling project**.

---

## Core Functionality: Multi-Fidelity Emulation

Multi-Fidelity Emulation is a technique used to predict the output of a computationally expensive, High-Fidelity (HF) model by leveraging data from a cheaper, Low-Fidelity (LF) counterpart. This method often involves a combination of data from both models to create a more efficient and accurate emulator than using the HF data alone.

---

## Workflow

The overall process is divided into two distinct workflows: one for developing and testing the methods using simulation data, and one for applying the methods to real-world climate model output.

### 1. Simulation Study Workflow (Method Development)

This workflow is designed for testing and comparing the tensor-based MF-GP method against baseline approaches. It uses synthetic data generated within the repository.

#### Stages

1.  **Data Preparation:**
    * **Action:** Run `simulation/setup.R`.
    * **Purpose:** Generates synthetic data (LF and HF) according to a predefined design.
    * **Output:** Data is saved in the subfolder `simulation/save/`.

2.  **Run Emulation Functions:**
    * Execute the scripts for the different emulation methods, all of which use the generated data from step 1.
    
    | Emulation Method | Description | Script |
    | :--- | :--- | :--- |
    | **Tensor-LF** | Low-Fidelity Tensor GP | `simulation/lf-tensor.R` |
    | **Tensor-HF** | High-Fidelity Tensor GP | `simulation/hf-tensor.R` |
    | **Tensor-MF** | **Multi-Fidelity Tensor GP (The main method)** | `simulation/mf-tensor.R` |
    | **Naïve-GP** | Standard GP (using only HF data for baseline comparison) | `simulation/naive-gp.R` |

3.  **Analysis:**
    * **Action:** Emulation results (e.g., predicted values, variance estimates) are automatically saved in `simulation/save/` for downstream analysis and comparison.

#### Custom Data Replacement

To test the emulation methods with your own synthetic or simple dataset, you can replace the standard simulation data by defining the following variables *before* running the emulation scripts (Step 2):

| Variable | Description | Type/Shape |
| :--- | :--- | :--- |
| `output_lf`, `output_hf` | Arrays containing the **Low-Fidelity (LF)** and **High-Fidelity (HF)** outputs. | Array |
| `design` | The full design matrix for all runs (both LF and HF). | Matrix |
| `runs_lf`, `runs_hf` | Row indices specifying which rows of `design` and `output` belong to the LF and HF fidelities, respectively. | Vector |
| `d_lf`, `d_hf` | The number of runs for each fidelity. (i.e., `length(runs_lf)` and `length(runs_hf)`). | Integer |

---

### 2. MPAS-Seaice Workflow (Application)

This workflow applies the developed tensor-based MF-GP methods to real-world output from the **Model for Prediction Across Scales (MPAS-Seaice)**. The primary difference from the simulation study is the need for external data loading and the use of **Leave-One-Out (LOO) Cross-Validation** for robust performance assessment.

#### Stages

1.  **Data Preparation:**
    * **Action:** Download and save the MPAS-Seaice data.
    * **Source:** The data can be downloaded from **XXX**.
    * **File Name:** Save the data as `mpas/save/mpas-data.Rdata`.

2.  **Run Emulation Functions (LOO Cross-Validation):**
    * Execute the scripts for the different emulation methods. These scripts are set up to run the LOO procedure.
    
    | Emulation Method | Script |
    | :--- | :--- |
    | **Tensor-LF** | `mpas/mpas-lf-tensor.R` |
    | **Tensor-HF** | `mpas/mpas-hf-tensor.R` |
    | **Tensor-MF** | `mpas/mpas-mf-tensor.R` |
    | **Naïve-GP** | `mpas/mpas-naive-gp.R` |

3.  **Analysis:**
    * **Output:** The Leave-One-Out emulation results are saved in `mpas/save/`.
    * **Structure:** Results are organized into subfolders for each emulation method (e.g., `mpas/save/mf-tensor/`).
