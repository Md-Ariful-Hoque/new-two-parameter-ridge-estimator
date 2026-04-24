# Simulation Design

## Objective

To evaluate estimator performance under multicollinearity.

## Data Generation

Predictors are generated with correlation structure:

x_j = √(1 − ρ²) z_j + ρ z_p

## Parameters

### Sample Sizes
- n = 30, 50, 100, 200

### Number of Predictors
- p = 3, 5, 10

### Correlation Levels
- ρ = 0.80, 0.90, 0.95, 0.99

### Error Variance
- σ = 1, 5, 10

### Replications
- 5000 simulations

## Outlier Study

- 25% contamination added
- extreme values injected into response

## Evaluation Metric

Mean Squared Error (MSE):

MSE = E[(β̂ − β)'(β̂ − β)]

## Output

- Tables comparing estimators
- Figures for MSE trends
