# Methodology Notes

## Linear Regression Model

We consider the standard linear regression model:

y = Xβ + ε

where:
- y is the response vector
- X is the design matrix
- β is the parameter vector
- ε ~ N(0, σ²I)

## Multicollinearity

Multicollinearity occurs when predictors in X are highly correlated.  
This leads to:
- unstable OLS estimates
- inflated variance
- poor inference

## Objective

To improve estimation under multicollinearity using:
- biased estimators
- shrinkage techniques
- two-parameter estimators

## Proposed Estimator

The proposed estimator (NTPR) modifies ridge-type estimators using:
- two parameters
- eigenvalue-based shrinkage

## Evaluation Criteria

Estimators are compared using:
- Mean Squared Error (MSE)
- robustness under outliers
- real data performance
