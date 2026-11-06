# Extract predictions

Help function to extract requested predictions from an object returned
by
[`pred_crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/pred_crossbasis.md).

## Usage

``` r
extract_predictions(
  x,
  slice_x = NULL,
  slice_y = NULL,
  slice_lag = NULL,
  slice_by = NULL,
  ptype = "allfit"
)
```

## Arguments

- x:

  An object returned by
  [`pred_crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/pred_crossbasis.md).

- slice_x:

  Exposure values of the first exposure variable at which to extract
  predictions. Integer or numeric vector.

- slice_y:

  Exposure values of the second exposure variable at which to extract
  predictions (bivariate exposures only). Integer or numeric vector.

- slice_lag:

  Lag value(s) at which to extract predictions. Integer or numeric
  vector.

- slice_by:

  Optional argument to specify the levels of the by-variable at which to
  extract predictions. One value or multiple values supplied through a
  vector.

- ptype:

  Character string specifying the type of prediction. One of:

  `"matfit"`

  :   Predictions for specified combinations of exposure and lag
      variables.

  `"matRRfit"`

  :   Predictions for exponentiated associations for specified
      combinations of exposure and lag variables.

  `"allfit"`

  :   Predictions for overall cumulative associations, default value.

  `"allRRfit"`

  :   Predictions for exponentiated associations for overall cumulative
      risks.

## Value

A dataframe containing the requested exposures (and lags), predicted
value, lower and upper bounds on the credible intervals (if requested in
predictions), exceedance probabilities (if requested in predictions),
and levels of the by-variable (if applicable).

## Examples

``` r
data(barcelona_temp)
data(shapefile_bcn)

# Examples of applications

dataset <- barcelona_temp |>
  dplyr::mutate(
    dow = weekdays(date),
    year = lubridate::year(date),
    month = lubridate::month(date)
  )

group <- factor(paste0(dataset$District, "-", dataset$year))


### Example 1: Penalized and common DLNM

# Create the (unpenalized) crossbasis

knots_max <- quantile(dataset$Tmax, probs = c(0.10,0.50,0.90))
knots_lag <- dlnm::logknots(0:14, nk = 2)

basis.max_unpen <- dlnm::crossbasis(
    dataset$Tmax,
    argvar=list(fun = "ns", knots = knots_max),
    arglag=list(fun="ns", knots = knots_lag),
    lag=14,
    group=group
)

# Fit model
model_max_conf <- dlnmLPS(
    dtot ~ 1 + factor(dow) +
    sm(month, arg = list(type = "ns", df = 3)) +
    sm(year, arg = list(type = "rw1")) +
    sm(depriv, arg = list(type = "ns", df = 3)) +
    sm(basis.max_unpen, by = depriv,
            arg.by = list(type = "ns", df = 3)) +
    sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
    data = dataset,
    offset = log(pop),
    model.comparison = list(DIC = TRUE),
    family = "Poisson")

# Make predictions
pred_conf_max <- pred_crossbasis(
    basis.max_unpen,
    model_max_conf,
    at_x = 22:30,
    cen = 25,
    lag = c(0,14),
    CI = TRUE,
    by = quantile(dataset$depriv,probs= c(0.10, 0.90))
)

# Extract overall RR for 28 and 30 degrees Celsius for both deprivation indices
extract_predictions(
    pred_conf_max,
    slice_x = c(28,30),
    ptype = "allRRfit"
)
#>    x      fit    lower    upper                 by
#> 1 28 1.236967 1.148185 1.332613 -0.929699483475094
#> 2 30 1.435575 1.158485 1.778942 -0.929699483475094
#> 3 28 1.217480 1.132777 1.308515   1.21474771640633
#> 4 30 1.620973 1.462294 1.796870   1.21474771640633

# Extract lag-specific RR for 28 and 30 degrees Celsius, lag 0 and 8 and
# for 10% deprivation percentile
extract_predictions(
    pred_conf_max,
    slice_x = c(28,30),
    slice_lag = c(0,8),
    slice_by = quantile(dataset$depriv,probs= c(0.10)),
    ptype = "matRRfit"
)
#>    x lag       fit     lower    upper                 by
#> 1 28   0 1.0840469 1.0422846 1.127483 -0.929699483475094
#> 2 30   0 1.2925590 1.1909438 1.402844 -0.929699483475094
#> 3 28   8 1.0175714 1.0054077 1.029882 -0.929699483475094
#> 4 30   8 0.9985841 0.9716917 1.026221 -0.929699483475094
```
