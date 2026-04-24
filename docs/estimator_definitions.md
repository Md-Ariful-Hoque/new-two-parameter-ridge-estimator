# Estimator Definitions

## OLS
β̂_OLS = (X'X)⁻¹X'y

## Ridge
β̂_Ridge = (X'X + kI)⁻¹X'y

## Liu Estimator
Combines ridge shrinkage with bias parameter d.

## Two-Parameter Estimators

Includes:
- LTE
- TP
- NBE
- YC
- MRT
- DK
- GTP
- LKL
- NTP
- NRT
- TPR

## Proposed Estimator (NTPR)

The NTPR estimator is defined as:

β̂_NTPR = q_min (X'X + kI)⁻¹ X'X β̂_OLS

where:
- q_min is the minimum eigenvalue adjustment factor

## Notes

- All estimators aim to reduce variance
- Trade-off between bias and variance is key
