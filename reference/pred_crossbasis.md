# Generate Predictions from a dlnmLPS Model

This function generates predictions from a fitted
[dlnmLPS](https://rutten-sara.github.io/dlnmLPS/reference/dlnmLPS.md)
model over a specified grid of exposure and lag values. Predictions are
computed relative to a user-defined reference (centering) value.

## Usage

``` r
pred_crossbasis(
  basis,
  model,
  at_x,
  at_y = NULL,
  cen,
  lag = NULL,
  by.lag = 1,
  CI.all = FALSE,
  CI = FALSE,
  ci.level = 0.95,
  exc.prob = FALSE,
  threshold = 0,
  mode = NULL,
  by = NULL
)
```

## Arguments

- basis:

  An object of class `"crossbasis"` or `"joint.crossbasis"` defining the
  exposure-lag basis used in the model.

- model:

  An object of class `"dlnmLPS"` returned by the
  [dlnmLPS](https://rutten-sara.github.io/dlnmLPS/reference/dlnmLPS.md)
  function.

- at_x:

  A numeric vector specifying constant exposure values over the selected
  lag period, or a matrix specifying exposure histories over the same
  lag period used during model estimation.

- at_y:

  A numeric vector or matrix analogous to `at_x`, specifying values or
  histories for the second exposure when `basis` is of class
  `"joint.crossbasis"`.

- cen:

  A value (or a vector of length 2 for `"joint.crossbasis"`) specifying
  the centering (reference) exposure level used for calculating relative
  risk predictions.

- lag:

  Either an integer scalar specifying the maximum lag, or a vector of
  length two specifying the lag range over which predictions are
  generated.

- by.lag:

  Numeric value specifying the increment between consecutive lag values
  used for prediction.

- CI.all:

  Logical value indicating whether credible intervals should be
  calculated only for the overall cumulative relative risk. Defaults to
  `FALSE`.

- CI:

  Logical value indicating whether credible intervals should be
  calculated for both individual lag-specific predictions and overall
  cumulative relative risks. Defaults to `FALSE`.

- ci.level:

  Numeric value specifying the credible interval level.

- exc.prob:

  Logical value indicating whether exceedance probabilities should be
  calculated. Defaults to `FALSE`.

- threshold:

  Numeric value defining the threshold used for exceedance probability
  calculations (at scale of linear predictor). Defaults to `0`.

- mode:

  Optional user-specified posterior mode.

- by:

  Level or value of the `by` variable at which the DLNM relationship is
  evaluated. For a continuous `by` variable, this is a (vector of)
  numeric value(s) defining the level at which the effect is evaluated.
  For a factor `by` variable, this is a (vector of) factor level(s)
  defining the corresponding stratum.

## Value

An object of class `"dlnmLPS_pred"` containing the following components
(structured by value of the `by`-variable, if specified):

- matfit:

  Matrix of predictions for the specified combinations of exposure and
  lag values.

- Qlower_matfit:

  Matrix of lower credible interval bounds for the values in `matfit`.

- Qupper_matfit:

  Matrix of upper credible interval bounds for the values in `matfit`.

- matRRfit:

  Matrix of exponentiated associations from `matfit`.

- Qlower_matRRfit:

  Matrix of exponentiated bounds from `Qlower_matfit`.

- Qupper_matRRfit:

  Matrix of exponentiated bounds from `Qupper_matfit`.

- allfit:

  Vector of estimated overall cumulative predictions.

- Qlower_allfit:

  Vector of lower credible interval bounds for the overall cumulative
  predictions in `allfit`.

- Qupper_allfit:

  Vector of upper credible interval bounds for the overall cumulative
  predictions in `allfit`.

- allRRfit:

  Matrix of exponentiated associations from `allfit`.

- Qlower_allRRfit:

  Matrix of exponentiated bounds from `Qlower_allfit`.

- Qupper_allRRfit:

  Matrix of exponentiated bounds from `Qupper_allfit`.

- exc_prob_matfit:

  Matrix of exceedance probabilities for matfit, if requested.

- exc_prob_allfit:

  Vector of exceedance probabilities for allfit, if requested.

## Examples

``` r
if (FALSE) { # \dontrun{
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

# Create the (penalized) crossbasis
basis.max <- dlnm::crossbasis(
    dataset$Tmax,
    argvar = list(df = 5, fun = "ps"),
    arglag = list(df = 5, fun = "ps"),
    lag = 14,
    group = group
)

# Fit model
model_max_common <- dlnmLPS(
    dtot ~ 1 + factor(dow) +
    sm(month, arg = list(type = "ns", df = 3)) +
    sm(year, arg = list(type = "rw1")) +
    sm(basis.max, arg = list(pen= 2)) +
    sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
    data = dataset,
    offset = log(pop),
    model.comparison = list(DIC = TRUE),
    family = "Poisson"
)

# Make predictions
pred_max <- pred_crossbasis(
     basis = basis.max,
     model = model_max_common,
     at_x = 22:30,
     cen = 25,
     lag = c(0,14),
     CI = TRUE
)

# Plot results (lag-specific and overall RR)
plot(pred_max, slice_lag = c(0,3,8), ptype = "slices", CI = TRUE, xlab = "temp", ylab = "RR")
plot(pred_max, ptype = "overall", CI = TRUE, xlab = "temp", ylab = "RR")


### Example 2: Penalized and common joint DLNM

# Create the (penalized) joint crossbasis
basis.max_min <- joint.crossbasis(
    cbind(dataset$Tmax, dataset$Tmin),
    argvar1 = list(df = 5, fun = "ps"),
    argvar2 = list(df = 5, fun = "ps"),
    arglag = list(df = 5, fun = "ps"),
    lag = 14,
    group = group
)

# Fit model
model_joint_common <- dlnmLPS(
    dtot ~ 1 + factor(dow) +
    sm(month, arg = list(type = "ns", df = 3)) +
    sm(year, arg = list(type = "rw1")) +
    sm(basis.max_min, arg = list(pen= 2)) +
    sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
    data = dataset,
    offset = log(pop),
    model.comparison = list(DIC = TRUE),
    family = "Poisson"
)

# Make predictions
pred_max_min <- pred_crossbasis(
     basis = basis.max_min,
     model = model_joint_common,
     at_x = c(29,29,29,30,30,30,31,31,31),
     at_y = c(20,21,22,20,21,22,20,21,22),
     cen = c(25,15),
     lag = c(0,14),
     CI = TRUE
)

# Plot results (contour plot overall RR)
contour(pred_max_min, ptype = "overall", cols = function(n) hcl.colors(n, "Inferno"))


### Example 3: Unpenalized and spatially varying DLNM

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
model_max_spat <- dlnmLPS(
    dtot ~ 1 + factor(dow) +
    sm(month, arg = list(type = "ns", df = 3)) +
    sm(year, arg = list(type = "rw1")) +
    sm(basis.max_unpen, by = District,
            arg.by = list(type = "Leroux", map = shapefile_bcn)) +
    sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
    data = dataset,
    offset = log(pop),
    model.comparison = list(DIC = TRUE),
    family = "Poisson")

summary(model_max_spat)

# Make predictions
pred_spat_max <- pred_crossbasis(
    basis.max_unpen,
    model_max_spat,
    at_x = 22:30,
    cen = 25,
    lag = c(0,14),
    CI = TRUE,
    by = unique(dataset$District)
)

# Plot results (region-specific overall RR)
plot(pred_spat_max, ptype = "overall", CI = FALSE, xlab = "temp", ylab = "RR", ylim = c(0.9,3.0),
    multiple.plots = FALSE)


### Example 4: Unpenalized and deprivation-varying DLNM

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

# Extract lag-specific RR for 28 and 30 degrees Celsius, lag 0 and 8 and
# for 10% deprivation percentile
extract_predictions(
    pred_conf_max,
    slice_x = c(28,30),
    slice_lag = c(0,8),
    slice_by = quantile(dataset$depriv,probs= c(0.10)),
    ptype = "matRRfit"
)

# Plot results (contour plot exposure-lag-response)
contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))
} # }
```
