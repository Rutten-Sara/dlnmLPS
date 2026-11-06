# Summarize a dlnmLPS model

Computes a summary of the fixed effects of a fitted `"dlnmLPS"` model.
The summary includes posterior means (modes), standard deviations, and
approximate 95% credible intervals based on the Laplace approximation.

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
