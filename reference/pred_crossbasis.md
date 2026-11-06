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
#> 
#> Call:
#> 
#> dtot ~ 1 + sm(month, arg = list(type = "ns", df = 3)) + sm(year, arg = list(type = "rw1")) + sm(basis.max_unpen, by = District, arg.by = list(type = "Leroux", map = shapefile_bcn)) + sm(District, arg = list(type = "Leroux", map = shapefile_bcn)) + factor(dow)
#> 
#> Fixed effects:
#> 
#>                                   mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept)                     -9.772 0.172    -10.110   -9.772     -9.434
#> month1                          -0.181 0.013     -0.207   -0.181     -0.155
#> month2                          -0.196 0.019     -0.233   -0.196     -0.158
#> month3                           0.000 0.011     -0.021    0.000      0.020
#> year2007                        -0.013 0.006     -0.025   -0.013     -0.001
#> year2008                         0.002 0.005     -0.008    0.002      0.013
#> year2009                         0.007 0.005     -0.004    0.007      0.018
#> year2010                         0.006 0.006     -0.005    0.006      0.017
#> year2011                        -0.003 0.007     -0.016   -0.003      0.010
#> basis.max_unpenv1.l1            -0.006 0.047     -0.098   -0.006      0.086
#> basis.max_unpenv1.l2            -0.070 0.046     -0.161   -0.070      0.020
#> basis.max_unpenv1.l3             0.046 0.047     -0.046    0.046      0.138
#> basis.max_unpenv1.l4            -0.102 0.046     -0.193   -0.102     -0.012
#> basis.max_unpenv2.l1             0.004 0.041     -0.076    0.004      0.084
#> basis.max_unpenv2.l2            -0.095 0.040     -0.173   -0.095     -0.016
#> basis.max_unpenv2.l3             0.101 0.042      0.020    0.101      0.183
#> basis.max_unpenv2.l4            -0.098 0.040     -0.176   -0.098     -0.020
#> basis.max_unpenv3.l1             0.007 0.090     -0.169    0.007      0.183
#> basis.max_unpenv3.l2            -0.268 0.090     -0.445   -0.268     -0.092
#> basis.max_unpenv3.l3             0.241 0.091      0.062    0.241      0.420
#> basis.max_unpenv3.l4            -0.257 0.090     -0.433   -0.257     -0.080
#> basis.max_unpenv4.l1            -0.009 0.046     -0.100   -0.009      0.082
#> basis.max_unpenv4.l2            -0.142 0.042     -0.224   -0.142     -0.060
#> basis.max_unpenv4.l3             0.346 0.049      0.251    0.346      0.441
#> basis.max_unpenv4.l4            -0.229 0.042     -0.312   -0.229     -0.146
#> basis.max_unpenv1.l1:District01 -0.041 0.033     -0.105   -0.041      0.024
#> basis.max_unpenv1.l1:District02 -0.005 0.029     -0.063   -0.005      0.052
#> basis.max_unpenv1.l1:District03  0.001 0.033     -0.064    0.001      0.065
#> basis.max_unpenv1.l1:District04 -0.023 0.035     -0.092   -0.023      0.046
#> basis.max_unpenv1.l1:District05  0.031 0.033     -0.033    0.031      0.096
#> basis.max_unpenv1.l1:District06  0.045 0.033     -0.020    0.045      0.109
#> basis.max_unpenv1.l1:District07  0.009 0.031     -0.051    0.009      0.069
#> basis.max_unpenv1.l1:District08 -0.013 0.033     -0.078   -0.013      0.053
#> basis.max_unpenv1.l1:District09  0.019 0.032     -0.044    0.019      0.081
#> basis.max_unpenv1.l1:District10 -0.023 0.033     -0.088   -0.023      0.041
#> basis.max_unpenv1.l2:District01 -0.043 0.034     -0.109   -0.043      0.023
#> basis.max_unpenv1.l2:District02  0.003 0.030     -0.056    0.003      0.062
#> basis.max_unpenv1.l2:District03  0.016 0.034     -0.050    0.016      0.082
#> basis.max_unpenv1.l2:District04  0.013 0.036     -0.058    0.013      0.084
#> basis.max_unpenv1.l2:District05  0.009 0.034     -0.057    0.009      0.075
#> basis.max_unpenv1.l2:District06 -0.022 0.034     -0.088   -0.022      0.045
#> basis.max_unpenv1.l2:District07 -0.002 0.031     -0.063   -0.002      0.059
#> basis.max_unpenv1.l2:District08  0.007 0.034     -0.059    0.007      0.073
#> basis.max_unpenv1.l2:District09  0.022 0.032     -0.042    0.022      0.085
#> basis.max_unpenv1.l2:District10 -0.003 0.034     -0.069   -0.003      0.063
#> basis.max_unpenv1.l3:District01  0.004 0.033     -0.061    0.004      0.069
#> basis.max_unpenv1.l3:District02  0.007 0.029     -0.051    0.007      0.065
#> basis.max_unpenv1.l3:District03  0.005 0.033     -0.059    0.005      0.070
#> basis.max_unpenv1.l3:District04 -0.002 0.035     -0.071   -0.002      0.067
#> basis.max_unpenv1.l3:District05  0.016 0.033     -0.048    0.016      0.081
#> basis.max_unpenv1.l3:District06  0.020 0.033     -0.045    0.020      0.085
#> basis.max_unpenv1.l3:District07 -0.006 0.030     -0.065   -0.006      0.054
#> basis.max_unpenv1.l3:District08 -0.019 0.034     -0.086   -0.019      0.047
#> basis.max_unpenv1.l3:District09 -0.010 0.032     -0.073   -0.010      0.054
#> basis.max_unpenv1.l3:District10 -0.016 0.033     -0.080   -0.016      0.048
#> basis.max_unpenv1.l4:District01  0.013 0.032     -0.050    0.013      0.076
#> basis.max_unpenv1.l4:District02 -0.004 0.029     -0.062   -0.004      0.054
#> basis.max_unpenv1.l4:District03  0.001 0.032     -0.062    0.001      0.064
#> basis.max_unpenv1.l4:District04  0.006 0.035     -0.063    0.006      0.075
#> basis.max_unpenv1.l4:District05  0.001 0.033     -0.063    0.001      0.066
#> basis.max_unpenv1.l4:District06 -0.002 0.032     -0.065   -0.002      0.061
#> basis.max_unpenv1.l4:District07 -0.016 0.030     -0.076   -0.016      0.043
#> basis.max_unpenv1.l4:District08  0.007 0.032     -0.057    0.007      0.071
#> basis.max_unpenv1.l4:District09  0.007 0.031     -0.054    0.007      0.068
#> basis.max_unpenv1.l4:District10 -0.013 0.033     -0.077   -0.013      0.051
#> basis.max_unpenv2.l1:District01  0.018 0.033     -0.046    0.018      0.082
#> basis.max_unpenv2.l1:District02 -0.004 0.030     -0.062   -0.004      0.053
#> basis.max_unpenv2.l1:District03  0.006 0.033     -0.059    0.006      0.071
#> basis.max_unpenv2.l1:District04 -0.001 0.035     -0.071   -0.001      0.068
#> basis.max_unpenv2.l1:District05 -0.020 0.033     -0.086   -0.020      0.045
#> basis.max_unpenv2.l1:District06 -0.018 0.033     -0.082   -0.018      0.046
#> basis.max_unpenv2.l1:District07 -0.018 0.031     -0.078   -0.018      0.042
#> basis.max_unpenv2.l1:District08  0.009 0.033     -0.056    0.009      0.074
#> basis.max_unpenv2.l1:District09  0.022 0.032     -0.042    0.022      0.085
#> basis.max_unpenv2.l1:District10  0.008 0.033     -0.056    0.008      0.072
#> basis.max_unpenv2.l2:District01 -0.019 0.033     -0.084   -0.019      0.046
#> basis.max_unpenv2.l2:District02  0.000 0.030     -0.059    0.000      0.059
#> basis.max_unpenv2.l2:District03 -0.003 0.034     -0.068   -0.003      0.063
#> basis.max_unpenv2.l2:District04 -0.004 0.036     -0.075   -0.004      0.067
#> basis.max_unpenv2.l2:District05 -0.007 0.034     -0.073   -0.007      0.060
#> basis.max_unpenv2.l2:District06  0.048 0.033     -0.017    0.048      0.113
#> basis.max_unpenv2.l2:District07  0.013 0.031     -0.048    0.013      0.074
#> basis.max_unpenv2.l2:District08 -0.013 0.033     -0.078   -0.013      0.051
#> basis.max_unpenv2.l2:District09  0.007 0.033     -0.057    0.007      0.072
#> basis.max_unpenv2.l2:District10 -0.022 0.033     -0.088   -0.022      0.044
#> basis.max_unpenv2.l3:District01  0.003 0.033     -0.061    0.003      0.067
#> basis.max_unpenv2.l3:District02 -0.011 0.029     -0.069   -0.011      0.046
#> basis.max_unpenv2.l3:District03 -0.004 0.033     -0.069   -0.004      0.061
#> basis.max_unpenv2.l3:District04 -0.014 0.035     -0.083   -0.014      0.054
#> basis.max_unpenv2.l3:District05 -0.025 0.033     -0.089   -0.025      0.040
#> basis.max_unpenv2.l3:District06 -0.008 0.033     -0.073   -0.008      0.056
#> basis.max_unpenv2.l3:District07  0.000 0.030     -0.060    0.000      0.059
#> basis.max_unpenv2.l3:District08  0.022 0.033     -0.043    0.022      0.088
#> basis.max_unpenv2.l3:District09  0.043 0.033     -0.021    0.043      0.107
#> basis.max_unpenv2.l3:District10 -0.005 0.033     -0.069   -0.005      0.059
#> basis.max_unpenv2.l4:District01  0.002 0.032     -0.060    0.002      0.064
#> basis.max_unpenv2.l4:District02  0.008 0.029     -0.049    0.008      0.065
#> basis.max_unpenv2.l4:District03  0.023 0.032     -0.039    0.023      0.085
#> basis.max_unpenv2.l4:District04  0.018 0.035     -0.050    0.018      0.086
#> basis.max_unpenv2.l4:District05  0.010 0.033     -0.054    0.010      0.074
#> basis.max_unpenv2.l4:District06 -0.007 0.031     -0.068   -0.007      0.055
#> basis.max_unpenv2.l4:District07 -0.014 0.030     -0.073   -0.014      0.045
#> basis.max_unpenv2.l4:District08 -0.017 0.032     -0.079   -0.017      0.045
#> basis.max_unpenv2.l4:District09 -0.014 0.031     -0.076   -0.014      0.047
#> basis.max_unpenv2.l4:District10 -0.009 0.032     -0.072   -0.009      0.053
#> basis.max_unpenv3.l1:District01 -0.003 0.038     -0.079   -0.003      0.072
#> basis.max_unpenv3.l1:District02  0.007 0.032     -0.055    0.007      0.069
#> basis.max_unpenv3.l1:District03 -0.014 0.038     -0.090   -0.014      0.061
#> basis.max_unpenv3.l1:District04  0.000 0.039     -0.077    0.000      0.077
#> basis.max_unpenv3.l1:District05 -0.012 0.036     -0.083   -0.012      0.059
#> basis.max_unpenv3.l1:District06  0.005 0.038     -0.070    0.005      0.080
#> basis.max_unpenv3.l1:District07  0.006 0.033     -0.058    0.006      0.071
#> basis.max_unpenv3.l1:District08 -0.004 0.042     -0.085   -0.004      0.078
#> basis.max_unpenv3.l1:District09  0.020 0.038     -0.055    0.020      0.094
#> basis.max_unpenv3.l1:District10 -0.004 0.036     -0.075   -0.004      0.066
#> basis.max_unpenv3.l2:District01  0.011 0.039     -0.065    0.011      0.086
#> basis.max_unpenv3.l2:District02  0.004 0.032     -0.059    0.004      0.066
#> basis.max_unpenv3.l2:District03 -0.007 0.039     -0.083   -0.007      0.069
#> basis.max_unpenv3.l2:District04 -0.014 0.040     -0.092   -0.014      0.063
#> basis.max_unpenv3.l2:District05 -0.006 0.037     -0.078   -0.006      0.065
#> basis.max_unpenv3.l2:District06  0.010 0.039     -0.065    0.010      0.086
#> basis.max_unpenv3.l2:District07 -0.002 0.033     -0.067   -0.002      0.063
#> basis.max_unpenv3.l2:District08  0.003 0.042     -0.079    0.003      0.086
#> basis.max_unpenv3.l2:District09  0.003 0.039     -0.073    0.003      0.078
#> basis.max_unpenv3.l2:District10 -0.002 0.036     -0.073   -0.002      0.070
#> basis.max_unpenv3.l3:District01 -0.012 0.038     -0.087   -0.012      0.062
#> basis.max_unpenv3.l3:District02  0.006 0.031     -0.055    0.006      0.068
#> basis.max_unpenv3.l3:District03 -0.013 0.038     -0.087   -0.013      0.062
#> basis.max_unpenv3.l3:District04 -0.008 0.039     -0.084   -0.008      0.068
#> basis.max_unpenv3.l3:District05 -0.005 0.036     -0.076   -0.005      0.065
#> basis.max_unpenv3.l3:District06  0.023 0.038     -0.051    0.023      0.097
#> basis.max_unpenv3.l3:District07  0.005 0.033     -0.059    0.005      0.069
#> basis.max_unpenv3.l3:District08 -0.007 0.041     -0.088   -0.007      0.074
#> basis.max_unpenv3.l3:District09  0.009 0.038     -0.065    0.009      0.083
#> basis.max_unpenv3.l3:District10  0.001 0.036     -0.069    0.001      0.071
#> basis.max_unpenv3.l4:District01 -0.014 0.039     -0.089   -0.014      0.062
#> basis.max_unpenv3.l4:District02  0.002 0.032     -0.060    0.002      0.064
#> basis.max_unpenv3.l4:District03 -0.011 0.039     -0.087   -0.011      0.065
#> basis.max_unpenv3.l4:District04 -0.012 0.039     -0.089   -0.012      0.065
#> basis.max_unpenv3.l4:District05 -0.007 0.036     -0.079   -0.007      0.064
#> basis.max_unpenv3.l4:District06  0.017 0.038     -0.059    0.017      0.092
#> basis.max_unpenv3.l4:District07  0.016 0.033     -0.049    0.016      0.080
#> basis.max_unpenv3.l4:District08  0.013 0.042     -0.070    0.013      0.095
#> basis.max_unpenv3.l4:District09 -0.008 0.038     -0.083   -0.008      0.067
#> basis.max_unpenv3.l4:District10  0.005 0.036     -0.066    0.005      0.076
#> basis.max_unpenv4.l1:District01  0.012 0.037     -0.061    0.012      0.084
#> basis.max_unpenv4.l1:District02 -0.009 0.032     -0.072   -0.009      0.054
#> basis.max_unpenv4.l1:District03  0.002 0.038     -0.072    0.002      0.075
#> basis.max_unpenv4.l1:District04 -0.007 0.039     -0.084   -0.007      0.070
#> basis.max_unpenv4.l1:District05 -0.007 0.037     -0.080   -0.007      0.066
#> basis.max_unpenv4.l1:District06 -0.018 0.037     -0.090   -0.018      0.053
#> basis.max_unpenv4.l1:District07 -0.002 0.033     -0.068   -0.002      0.063
#> basis.max_unpenv4.l1:District08  0.003 0.038     -0.071    0.003      0.077
#> basis.max_unpenv4.l1:District09  0.006 0.039     -0.071    0.006      0.083
#> basis.max_unpenv4.l1:District10  0.021 0.036     -0.049    0.021      0.091
#> basis.max_unpenv4.l2:District01  0.008 0.038     -0.066    0.008      0.082
#> basis.max_unpenv4.l2:District02 -0.003 0.032     -0.066   -0.003      0.060
#> basis.max_unpenv4.l2:District03 -0.017 0.038     -0.092   -0.017      0.057
#> basis.max_unpenv4.l2:District04 -0.006 0.040     -0.084   -0.006      0.072
#> basis.max_unpenv4.l2:District05 -0.001 0.037     -0.074   -0.001      0.071
#> basis.max_unpenv4.l2:District06  0.012 0.037     -0.061    0.012      0.085
#> basis.max_unpenv4.l2:District07  0.010 0.033     -0.055    0.010      0.076
#> basis.max_unpenv4.l2:District08  0.009 0.038     -0.066    0.009      0.084
#> basis.max_unpenv4.l2:District09 -0.015 0.039     -0.092   -0.015      0.062
#> basis.max_unpenv4.l2:District10  0.003 0.036     -0.068    0.003      0.074
#> basis.max_unpenv4.l3:District01 -0.002 0.037     -0.075   -0.002      0.070
#> basis.max_unpenv4.l3:District02 -0.009 0.032     -0.071   -0.009      0.054
#> basis.max_unpenv4.l3:District03  0.001 0.038     -0.074    0.001      0.076
#> basis.max_unpenv4.l3:District04 -0.001 0.039     -0.077   -0.001      0.076
#> basis.max_unpenv4.l3:District05 -0.002 0.037     -0.075   -0.002      0.070
#> basis.max_unpenv4.l3:District06 -0.012 0.037     -0.085   -0.012      0.060
#> basis.max_unpenv4.l3:District07 -0.001 0.033     -0.066   -0.001      0.064
#> basis.max_unpenv4.l3:District08  0.010 0.039     -0.066    0.010      0.086
#> basis.max_unpenv4.l3:District09  0.023 0.039     -0.054    0.023      0.101
#> basis.max_unpenv4.l3:District10 -0.006 0.036     -0.076   -0.006      0.064
#> basis.max_unpenv4.l4:District01 -0.011 0.037     -0.083   -0.011      0.061
#> basis.max_unpenv4.l4:District02  0.008 0.032     -0.054    0.008      0.071
#> basis.max_unpenv4.l4:District03  0.012 0.037     -0.061    0.012      0.085
#> basis.max_unpenv4.l4:District04  0.021 0.039     -0.056    0.021      0.097
#> basis.max_unpenv4.l4:District05  0.020 0.037     -0.053    0.020      0.092
#> basis.max_unpenv4.l4:District06  0.018 0.036     -0.052    0.018      0.088
#> basis.max_unpenv4.l4:District07 -0.008 0.033     -0.073   -0.008      0.057
#> basis.max_unpenv4.l4:District08 -0.013 0.037     -0.086   -0.013      0.059
#> basis.max_unpenv4.l4:District09 -0.025 0.039     -0.101   -0.025      0.051
#> basis.max_unpenv4.l4:District10 -0.020 0.036     -0.090   -0.020      0.050
#> District01                      -0.115 0.127     -0.363   -0.115      0.133
#> District02                       0.077 0.112     -0.142    0.077      0.297
#> District03                      -0.011 0.126     -0.258   -0.011      0.236
#> District04                      -0.049 0.128     -0.300   -0.049      0.202
#> District05                      -0.033 0.122     -0.273   -0.033      0.207
#> District06                       0.202 0.126     -0.045    0.202      0.450
#> District07                       0.009 0.115     -0.216    0.009      0.235
#> District08                       0.016 0.134     -0.247    0.016      0.278
#> District09                      -0.003 0.125     -0.248   -0.003      0.243
#> District10                      -0.094 0.121     -0.332   -0.094      0.144
#> factor(dow)Monday                0.015 0.012     -0.009    0.015      0.038
#> factor(dow)Saturday              0.018 0.012     -0.006    0.018      0.042
#> factor(dow)Sunday                0.020 0.012     -0.004    0.020      0.044
#> factor(dow)Thursday              0.002 0.012     -0.022    0.002      0.026
#> factor(dow)Tuesday               0.009 0.012     -0.015    0.009      0.032
#> factor(dow)Wednesday             0.012 0.012     -0.011    0.012      0.036
#>                                   mode
#> (Intercept)                     -9.772
#> month1                          -0.181
#> month2                          -0.196
#> month3                           0.000
#> year2007                        -0.013
#> year2008                         0.002
#> year2009                         0.007
#> year2010                         0.006
#> year2011                        -0.003
#> basis.max_unpenv1.l1            -0.006
#> basis.max_unpenv1.l2            -0.070
#> basis.max_unpenv1.l3             0.046
#> basis.max_unpenv1.l4            -0.102
#> basis.max_unpenv2.l1             0.004
#> basis.max_unpenv2.l2            -0.095
#> basis.max_unpenv2.l3             0.101
#> basis.max_unpenv2.l4            -0.098
#> basis.max_unpenv3.l1             0.007
#> basis.max_unpenv3.l2            -0.268
#> basis.max_unpenv3.l3             0.241
#> basis.max_unpenv3.l4            -0.257
#> basis.max_unpenv4.l1            -0.009
#> basis.max_unpenv4.l2            -0.142
#> basis.max_unpenv4.l3             0.346
#> basis.max_unpenv4.l4            -0.229
#> basis.max_unpenv1.l1:District01 -0.041
#> basis.max_unpenv1.l1:District02 -0.005
#> basis.max_unpenv1.l1:District03  0.001
#> basis.max_unpenv1.l1:District04 -0.023
#> basis.max_unpenv1.l1:District05  0.031
#> basis.max_unpenv1.l1:District06  0.045
#> basis.max_unpenv1.l1:District07  0.009
#> basis.max_unpenv1.l1:District08 -0.013
#> basis.max_unpenv1.l1:District09  0.019
#> basis.max_unpenv1.l1:District10 -0.023
#> basis.max_unpenv1.l2:District01 -0.043
#> basis.max_unpenv1.l2:District02  0.003
#> basis.max_unpenv1.l2:District03  0.016
#> basis.max_unpenv1.l2:District04  0.013
#> basis.max_unpenv1.l2:District05  0.009
#> basis.max_unpenv1.l2:District06 -0.022
#> basis.max_unpenv1.l2:District07 -0.002
#> basis.max_unpenv1.l2:District08  0.007
#> basis.max_unpenv1.l2:District09  0.022
#> basis.max_unpenv1.l2:District10 -0.003
#> basis.max_unpenv1.l3:District01  0.004
#> basis.max_unpenv1.l3:District02  0.007
#> basis.max_unpenv1.l3:District03  0.005
#> basis.max_unpenv1.l3:District04 -0.002
#> basis.max_unpenv1.l3:District05  0.016
#> basis.max_unpenv1.l3:District06  0.020
#> basis.max_unpenv1.l3:District07 -0.006
#> basis.max_unpenv1.l3:District08 -0.019
#> basis.max_unpenv1.l3:District09 -0.010
#> basis.max_unpenv1.l3:District10 -0.016
#> basis.max_unpenv1.l4:District01  0.013
#> basis.max_unpenv1.l4:District02 -0.004
#> basis.max_unpenv1.l4:District03  0.001
#> basis.max_unpenv1.l4:District04  0.006
#> basis.max_unpenv1.l4:District05  0.001
#> basis.max_unpenv1.l4:District06 -0.002
#> basis.max_unpenv1.l4:District07 -0.016
#> basis.max_unpenv1.l4:District08  0.007
#> basis.max_unpenv1.l4:District09  0.007
#> basis.max_unpenv1.l4:District10 -0.013
#> basis.max_unpenv2.l1:District01  0.018
#> basis.max_unpenv2.l1:District02 -0.004
#> basis.max_unpenv2.l1:District03  0.006
#> basis.max_unpenv2.l1:District04 -0.001
#> basis.max_unpenv2.l1:District05 -0.020
#> basis.max_unpenv2.l1:District06 -0.018
#> basis.max_unpenv2.l1:District07 -0.018
#> basis.max_unpenv2.l1:District08  0.009
#> basis.max_unpenv2.l1:District09  0.022
#> basis.max_unpenv2.l1:District10  0.008
#> basis.max_unpenv2.l2:District01 -0.019
#> basis.max_unpenv2.l2:District02  0.000
#> basis.max_unpenv2.l2:District03 -0.003
#> basis.max_unpenv2.l2:District04 -0.004
#> basis.max_unpenv2.l2:District05 -0.007
#> basis.max_unpenv2.l2:District06  0.048
#> basis.max_unpenv2.l2:District07  0.013
#> basis.max_unpenv2.l2:District08 -0.013
#> basis.max_unpenv2.l2:District09  0.007
#> basis.max_unpenv2.l2:District10 -0.022
#> basis.max_unpenv2.l3:District01  0.003
#> basis.max_unpenv2.l3:District02 -0.011
#> basis.max_unpenv2.l3:District03 -0.004
#> basis.max_unpenv2.l3:District04 -0.014
#> basis.max_unpenv2.l3:District05 -0.025
#> basis.max_unpenv2.l3:District06 -0.008
#> basis.max_unpenv2.l3:District07  0.000
#> basis.max_unpenv2.l3:District08  0.022
#> basis.max_unpenv2.l3:District09  0.043
#> basis.max_unpenv2.l3:District10 -0.005
#> basis.max_unpenv2.l4:District01  0.002
#> basis.max_unpenv2.l4:District02  0.008
#> basis.max_unpenv2.l4:District03  0.023
#> basis.max_unpenv2.l4:District04  0.018
#> basis.max_unpenv2.l4:District05  0.010
#> basis.max_unpenv2.l4:District06 -0.007
#> basis.max_unpenv2.l4:District07 -0.014
#> basis.max_unpenv2.l4:District08 -0.017
#> basis.max_unpenv2.l4:District09 -0.014
#> basis.max_unpenv2.l4:District10 -0.009
#> basis.max_unpenv3.l1:District01 -0.003
#> basis.max_unpenv3.l1:District02  0.007
#> basis.max_unpenv3.l1:District03 -0.014
#> basis.max_unpenv3.l1:District04  0.000
#> basis.max_unpenv3.l1:District05 -0.012
#> basis.max_unpenv3.l1:District06  0.005
#> basis.max_unpenv3.l1:District07  0.006
#> basis.max_unpenv3.l1:District08 -0.004
#> basis.max_unpenv3.l1:District09  0.020
#> basis.max_unpenv3.l1:District10 -0.004
#> basis.max_unpenv3.l2:District01  0.011
#> basis.max_unpenv3.l2:District02  0.004
#> basis.max_unpenv3.l2:District03 -0.007
#> basis.max_unpenv3.l2:District04 -0.014
#> basis.max_unpenv3.l2:District05 -0.006
#> basis.max_unpenv3.l2:District06  0.010
#> basis.max_unpenv3.l2:District07 -0.002
#> basis.max_unpenv3.l2:District08  0.003
#> basis.max_unpenv3.l2:District09  0.003
#> basis.max_unpenv3.l2:District10 -0.002
#> basis.max_unpenv3.l3:District01 -0.012
#> basis.max_unpenv3.l3:District02  0.006
#> basis.max_unpenv3.l3:District03 -0.013
#> basis.max_unpenv3.l3:District04 -0.008
#> basis.max_unpenv3.l3:District05 -0.005
#> basis.max_unpenv3.l3:District06  0.023
#> basis.max_unpenv3.l3:District07  0.005
#> basis.max_unpenv3.l3:District08 -0.007
#> basis.max_unpenv3.l3:District09  0.009
#> basis.max_unpenv3.l3:District10  0.001
#> basis.max_unpenv3.l4:District01 -0.014
#> basis.max_unpenv3.l4:District02  0.002
#> basis.max_unpenv3.l4:District03 -0.011
#> basis.max_unpenv3.l4:District04 -0.012
#> basis.max_unpenv3.l4:District05 -0.007
#> basis.max_unpenv3.l4:District06  0.017
#> basis.max_unpenv3.l4:District07  0.016
#> basis.max_unpenv3.l4:District08  0.013
#> basis.max_unpenv3.l4:District09 -0.008
#> basis.max_unpenv3.l4:District10  0.005
#> basis.max_unpenv4.l1:District01  0.012
#> basis.max_unpenv4.l1:District02 -0.009
#> basis.max_unpenv4.l1:District03  0.002
#> basis.max_unpenv4.l1:District04 -0.007
#> basis.max_unpenv4.l1:District05 -0.007
#> basis.max_unpenv4.l1:District06 -0.018
#> basis.max_unpenv4.l1:District07 -0.002
#> basis.max_unpenv4.l1:District08  0.003
#> basis.max_unpenv4.l1:District09  0.006
#> basis.max_unpenv4.l1:District10  0.021
#> basis.max_unpenv4.l2:District01  0.008
#> basis.max_unpenv4.l2:District02 -0.003
#> basis.max_unpenv4.l2:District03 -0.017
#> basis.max_unpenv4.l2:District04 -0.006
#> basis.max_unpenv4.l2:District05 -0.001
#> basis.max_unpenv4.l2:District06  0.012
#> basis.max_unpenv4.l2:District07  0.010
#> basis.max_unpenv4.l2:District08  0.009
#> basis.max_unpenv4.l2:District09 -0.015
#> basis.max_unpenv4.l2:District10  0.003
#> basis.max_unpenv4.l3:District01 -0.002
#> basis.max_unpenv4.l3:District02 -0.009
#> basis.max_unpenv4.l3:District03  0.001
#> basis.max_unpenv4.l3:District04 -0.001
#> basis.max_unpenv4.l3:District05 -0.002
#> basis.max_unpenv4.l3:District06 -0.012
#> basis.max_unpenv4.l3:District07 -0.001
#> basis.max_unpenv4.l3:District08  0.010
#> basis.max_unpenv4.l3:District09  0.023
#> basis.max_unpenv4.l3:District10 -0.006
#> basis.max_unpenv4.l4:District01 -0.011
#> basis.max_unpenv4.l4:District02  0.008
#> basis.max_unpenv4.l4:District03  0.012
#> basis.max_unpenv4.l4:District04  0.021
#> basis.max_unpenv4.l4:District05  0.020
#> basis.max_unpenv4.l4:District06  0.018
#> basis.max_unpenv4.l4:District07 -0.008
#> basis.max_unpenv4.l4:District08 -0.013
#> basis.max_unpenv4.l4:District09 -0.025
#> basis.max_unpenv4.l4:District10 -0.020
#> District01                      -0.115
#> District02                       0.077
#> District03                      -0.011
#> District04                      -0.049
#> District05                      -0.033
#> District06                       0.202
#> District07                       0.009
#> District08                       0.016
#> District09                      -0.003
#> District10                      -0.094
#> factor(dow)Monday                0.015
#> factor(dow)Saturday              0.018
#> factor(dow)Sunday                0.020
#> factor(dow)Thursday              0.002
#> factor(dow)Tuesday               0.009
#> factor(dow)Wednesday             0.012
#> 
#> Hyperparameters:
#> 
#> Precision Leroux re District  21.12
#> Precision RW year  4842
#> Precision spatial DLNM basis.max_unpen  388.7
#> Rho DLNM basis.max_unpen  0.5129
#> Rho Leroux re District  0.4606
#> 
#> Model comparison:
#> 
#> DIC: 30375.69 ; Effective Parameters: 82.35752 

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

# Plot results (region-specific overal RR)
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
#>    x      fit    lower    upper                 by
#> 1 28 1.236969 1.148188 1.332615 -0.929699483475094
#> 2 30 1.435580 1.158489 1.778947 -0.929699483475094
#> 3 28 1.217482 1.132780 1.308517   1.21474771640633
#> 4 30 1.620977 1.462298 1.796875   1.21474771640633

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
#> 1 28   0 1.0840473 1.0422849 1.127483 -0.929699483475094
#> 2 30   0 1.2925594 1.1909441 1.402845 -0.929699483475094
#> 3 28   8 1.0175714 1.0054077 1.029882 -0.929699483475094
#> 4 30   8 0.9985843 0.9716919 1.026221 -0.929699483475094

# Plot results (contour plot exposure-lag-response)
contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))



```
