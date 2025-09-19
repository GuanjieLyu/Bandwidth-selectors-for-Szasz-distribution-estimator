# Smoothing Parameter Selection in Szász–Mirakyan Operator–Based Distribution Estimation

This repository contains the code and supplementary material for the paper:

**“Smoothing Parameter Selection in Szász–Mirakyan Operator–Based Distribution Estimation”**

---

## Overview

The paper develops and investigates two data–driven strategies for selecting the smoothing parameter in Szász–Mirakyan operator–based distribution estimation:

1. **Bootstrap–Based Indirect Method**  
   - Employs bootstrap resampling to approximate the analytic mean squared error (MSE) expansion.  
   - Provides a regression–based approximation to the optimal smoothing level.  

2. **Estimation–Based Direct Method**  
   - Constructs plug–in estimates of bias and variance using kernel–smoothed distribution estimators.  
   - Minimizes the estimated MSE over a candidate grid of smoothing levels.  

Both approaches are analyzed theoretically and evaluated empirically through simulation studies and real–data applications.

---

## Repository Contents

- **R Code**:  
  Functions implementing the Szász–Mirakyan distribution estimator, indirect (bootstrap) selector, and direct (estimation–based) selector.  

- **Simulation Scripts**:  
  R scripts to replicate the simulation study, including bias/variance evaluation, MSE curves, and comparisons between selectors.  

- **Data Example**:  
  Code illustrating the application to the `faithful` eruption dataset.  

- **Figures**:  
  Scripts for generating the figures reported in the manuscript (CDF curves, selector comparisons, simulation results).  

---

## Requirements

- **R (≥ 4.0.0)**  
- Packages: `ggplot2`, `tidyr`, `dplyr`, `np` (for kernel distribution estimation), and other standard R libraries.

---

## How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
   cd <repo-name>
