# Contour plots of predicted outcomes

Produce contour plots of predicted outcomes obtained from
[`pred_crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/pred_crossbasis.md).
In contrast to
[`plot.dlnmLPS_pred()`](https://rutten-sara.github.io/dlnmLPS/reference/plot.DLNMLPS_pred.md),
this function only produces contour plots of surfaces and therefore is
only available for plot configurations that result in a
three-dimensional surface. For Poisson or negative binomial models with
a log link, relative risks (RRs) are plotted; for binomial models with a
logit link, odds ratios (ORs) are plotted; and for Gaussian models with
an identity link, additive differences in the outcome are plotted.

## Usage

``` r
# S3 method for class 'dlnmLPS_pred'
contour(
  x,
  slice_lag = NULL,
  slice_x = NULL,
  slice_y = NULL,
  slice_by = NULL,
  ptype = "overall",
  xlim = NULL,
  ylim = NULL,
  zlim = NULL,
  cols = NULL,
  xlab = NULL,
  ylab = NULL,
  zlab = NULL,
  mar = NULL,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`pred_crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/pred_crossbasis.md).

- slice_lag:

  Lag value(s) at which to obtain slices. Integer or numeric vector.

- slice_x:

  Exposure values of the first exposure variable at which to obtain
  slices. Integer or numeric vector.

- slice_y:

  Exposure values of the second exposure variable at which to obtain
  slices (bivariate exposures only). Integer or numeric vector.

- slice_by:

  Optional argument to specify the levels of the by-variable to be
  plotted. One value or multiple values supplied through a vector.

- ptype:

  Character string specifying the type of plot. One of:

  `"overall"`

  :   Overall cumulative exposure-response plot.

  `"cross"`

  :   Complete exposure-lag-response surface (only available for
      univariate exposures).

  `"slices"`

  :   Slices through the exposure-lag-response surface.

- xlim, ylim, zlim:

  Numeric vectors of length two giving axis limits.

- cols:

  Colours used for lines, surfaces or contours.

- xlab, ylab, zlab:

  Axis labels.

- mar:

  Figure margins.

- ...:

  Additional graphical arguments passed to the underlying plotting
  functions. For 3D plots, additional arguments of
  [`graphics::persp`](https://rdrr.io/r/graphics/persp.html) can be
  passed.

## Value

Produces a contour plot and returns `NULL` invisibly.

## Details

This function provides contour representations of the three-dimensional
surfaces available in
[`plot.dlnmLPS_pred()`](https://rutten-sara.github.io/dlnmLPS/reference/plot.DLNMLPS_pred.md).

The available plots are:

- the exposure-lag-response surface (`ptype = "cross"`) for univariate
  exposures;

- the overall cumulative exposure-response surface (`ptype = "overall"`)
  for bivariate exposures when no exposure slice is specified;

- slice surfaces (`ptype = "slices"`) whenever the selected slicing
  variables define a three-dimensional surface.

Contour plots are not available for two-dimensional plots (e.g.,
cumulative exposure-response curves or exposure/lag slices shown as
multiple lines). Confidence intervals are therefore not supported.

## See also

[`plot.dlnmLPS_pred()`](https://rutten-sara.github.io/dlnmLPS/reference/plot.DLNMLPS_pred.md),
[`pred_crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/pred_crossbasis.md)

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


### Example 1: Unpenalized and spatially varying DLNM

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
plot(pred_spat_max, ptype = "overall", CI = TRUE, xlab = "temp", ylab = "RR", ylim = c(0.9,3.0))


### Example 2: Unpenalized and deprivation-varying DLNM

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

# Plot results (contour plot exposure-lag-response)
contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))


### Example 3: Penalized and common DLNM

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


### Example 4: Penalized and common joint DLNM

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
} # }
```
