# New Two-Parameter Ridge Estimator (NTPR)

## Overview
This repository contains R code for the paper:

**“New Two-Parameter Ridge Estimators for Addressing Multicollinearity in Linear Regression: Theory, Simulation, and Applications.”**

The project implements:
- classical and two-parameter estimators
- Monte Carlo simulation study
- outlier robustness analysis
- real-data application (body fat dataset)

---

## Motivation
Multicollinearity causes instability in ordinary least squares (OLS) estimation by inflating variance and reducing reliability.

This project introduces and evaluates a **New Two-Parameter Ridge Estimator (NTPR)** designed to improve estimation performance under high correlation among predictors.

---

## Main Features

- Implementation of:
  - OLS, Ridge, Liu
  - Multiple two-parameter estimators
  - Proposed NTPR estimator
- Monte Carlo simulation under various scenarios
- Outlier contamination study
- Real dataset application
- Reproducible workflow

---

## Repository Structure

```text
R/
  data_generation.R
  estimators.R
  estimators_two_parameter.R
  mse_functions.R
  simulation_core.R
  outlier_functions.R
  real_data_analysis.R

scripts/
  run_simulation_main.R
  run_simulation_outliers.R
  run_real_data_application.R
  make_tables_figures.R

data/
  raw/
  processed/
  README.md

results/
  tables/
  figures/
  saved_objects/

docs/
  methodology_notes.md
  estimator_definitions.md
  simulation_design.md

paper/
  published_paper.pdf
