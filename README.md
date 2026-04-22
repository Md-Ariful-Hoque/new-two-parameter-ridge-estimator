# New Two-Parameter Ridge Estimator for Multicollinearity

## Overview
This repository contains R code for the paper:
“New Two-Parameter Ridge Estimators for Addressing Multicollinearity in Linear Regression: Theory, Simulation, and Applications”.

## Contents
- Theoretical estimator implementation
- Monte Carlo simulation study
- Outlier contamination study
- Real-data application

## Simulation Design
- Sample sizes: 30, 50, 100, 200
- Predictors: 3, 5, 10
- Correlation: 0.80, 0.90, 0.95, 0.99
- Error SD: 1, 5, 10
- Replications: 5000

## Main Estimators
- OLS
- TPR
- NTPR
- Other two-parameter estimators

## Repository Structure
- `R/` reusable methods
- `scripts/` executable workflows
- `data/` raw and processed data
- `results/` tables and figures
- `paper/` publication files

## Author
Md Ariful Hoque
