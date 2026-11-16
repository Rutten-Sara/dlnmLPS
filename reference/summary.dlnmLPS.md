# Summarize a dlnmLPS model

Summarizes the posterior of the model coefficients of the `"dlnmLPS"`
model (mean = mode under the Laplace approximation, SD, 2.5/50/97.5%
quantiles) and of the hyperparameters.

## Usage

``` r
# S3 method for class 'dlnmLPS'
summary(object, ...)
```

## Arguments

- object:

  A fitted `"dlnmLPS"` model.

- ...:

  Additional arguments to be passed.

## Value

An object of class `"summary.dlnmLPS"`.
