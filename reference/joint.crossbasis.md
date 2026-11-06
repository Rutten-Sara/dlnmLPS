# Construct a joint crossbasis matrix for two interacting exposures

Constructs a distributed lag non-linear model (DLNM) crossbasis for the
interaction between two exposure variables. The exposures can be
supplied either as concurrent time series or as pre-computed lagged
exposure histories.

## Usage

``` r
joint.crossbasis(
  x = NULL,
  lag = 10,
  arglag = list(),
  argvar1 = list(),
  argvar2 = list(),
  group = NULL,
  is_identifiable = TRUE
)
```

## Arguments

- x:

  A matrix or a list containing the exposure variables. If a matrix, the
  first and second columns contain the time series of the first and
  second exposure, respectively. If a list, `x[[1]]` and `x[[2]]` should
  contain the lagged exposure histories for the first and second
  exposure.

- lag:

  Integer specifying the maximum lag.

- arglag:

  A named list of arguments defining the basis for the lag dimension
  (e.g., `fun`, `knots`, `nk`, `df`, `degree`, `Boundary.knots`).

- argvar1:

  A named list of arguments defining the basis for the first exposure
  dimension (e.g., `fun`, `knots`, `nk`, `df`, `degree`,
  `Boundary.knots`).

- argvar2:

  A named list of arguments defining the basis for the second exposure
  dimension (e.g., `fun`, `knots`, `nk`, `df`, `degree`,
  `Boundary.knots`).

- group:

  An optional factor or vector identifying groups of observations
  belonging to independent time series.

- is_identifiable:

  Logical; if `TRUE` (default), constructs an identifiable basis by
  removing confounded terms.

## Value

A joint crossbasis matrix representing the interaction between two
exposure variables across the specified lag period.
