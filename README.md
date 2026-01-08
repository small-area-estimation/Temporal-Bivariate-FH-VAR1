# TBFH-VAR1: Temporal Bivariate Fay-Herriot VAR1 Model

This repository contains R codes for the paper:

> **Cabello, E., Esteban, M.D., Morales, D. et al. Small area estimation of poverty indicators under bivariate Fay–Herriot model with correlated time effects. *AStA Adv Stat Anal (2025)*. [https://doi.org/10.1007/s10182-025-00550-5](https://doi.org/10.1007/s10182-025-00550-5)**

---

## 📄 Overview

The TBFH-VAR1 model is a generalization of the bivariate Fay-Herriot model proposed by Benavent and Morales (2021). 
It introduces a temporal dimension by modelling domain-level random effects with a first order vector autoregressive correlation structure. 
This allows the model to capture target variable correlations and domain-level random effects autocorrelation, leading to improved precision and stability in small area estimates from longitudinal surveys.

This repository includes:

- R functions to:
  - Fit the TBFH-VAR1 model
  - EBLUPs
  - MSE estimation (analytical & bootstrap)
  - Example of use


---
