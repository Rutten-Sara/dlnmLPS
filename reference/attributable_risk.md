# Attributable risk from dlnmLPS model

This function calculates the attributable risk (attributable fraction
and attributable number) from a fitted
[dlnmLPS](https://rutten-sara.github.io/dlnmLPS/reference/dlnmLPS.md)
model.

## Usage

``` r
attributable_risk(
  basis,
  cases,
  model,
  cen,
  x,
  y = NULL,
  by = NULL,
  group = NULL,
  dir = "back",
  tot = TRUE,
  range_x = NULL,
  range_y = NULL,
  samples = FALSE,
  nsamples = 1000,
  exc.threshold = NULL
)
```

## Arguments

- basis:

  An object of class `"crossbasis"` or `"joint.crossbasis"` defining the
  exposure-lag basis used in the model.

- cases:

  A vector of cases or, only if `dir = "forw"`, a matrix of future
  cases.

- model:

  An object of class `"dlnmLPS"` returned by the
  [dlnmLPS](https://rutten-sara.github.io/dlnmLPS/reference/dlnmLPS.md)
  function.

- cen:

  A value (or a vector of length 2 for `"joint.crossbasis"`) specifying
  the centering (reference) exposure level used for calculating relative
  risk predictions.

- x:

  Vector of the first exposure or, only if `dir = "back"` a matrix of
  lagged exposures, for which the attributable fraction needs to be
  calculated.

- y:

  Vector of the second exposure or, only if `dir = "back"` a matrix of
  lagged exposures, for which the attributable fraction needs to be
  calculated. Only for `basis` of class`"joint.crossbasis"`.

- by:

  Vector of the `by` variable-levels (same length as cases) at which the
  DLNM relationship is evaluated. For a continuous `by` variable, this
  is a vector of numeric values defining the level at which the effect
  is evaluated. For a factor `by` variable, this is a vector of factor
  levels defining the corresponding stratum.

- group:

  Vector of region indicators to specify an optional grouping structure
  of the output. For spatially varying DLNMs, this argument
  automatically corresponds to the ID vector in `by` (i.e.,
  region-specific attributable risks are calculated). For non-spatially
  varying DLNMs, a vector can also be specified (same length as cases)
  to group the output per region.

- dir:

  Either "back" (default) or "forw" for backward and forward perspective
  on attributable risk respectively.

- tot:

  If `True`, the total attributable risk over the entire period is
  calculated.

- range_x:

  The range of the first exposure for which the attributable risk is
  calculated. If `NULL`(default), the whole range is used.

- range_y:

  The range of the second exposure for which the attributable risk is
  calculated. If `NULL`(default), the whole range is used. Only for
  `basis` of class`"joint.crossbasis"`.

- samples:

  Set to `TRUE` if samples from the posterior distribution should be
  returned. Samples can only be calculated for `tot = TRUE`.

- nsamples:

  Number of samples from the posterior distribution to be returned.
  (only when `samples = TRUE`).

- exc.threshold:

  Numeric value defining the threshold used for exceedance probability
  calculations of the attributable fraction. If `NULL` (default), no
  exceedance probabilities are calculated.

## Value

A list containing the attributable fractions and numbers.

If `samples = TRUE`, the list is structured in two main elements. Note
that this option is only possible with `tot = TRUE`. The two elements
are:

- samples element:

  List with samples from the posterior distribution of the total
  attributable fraction `af` and total attributable number `an`. The
  group `ID` (if applicable) and exceedance probabilities (if requested)
  are also reported. For non-NULL grouping structures, this output is
  structured per group.

- fit element:

  The fitted total `af`, `an`, and the group `ID` if applicable. For
  non-NULL grouping structures, this output is structured per group.

If `samples = FALSE`, the list contains the fitted total `af`, `an`, and
the group `ID` if applicable, when `tot = TRUE`. When `tot = FALSE`,
this list contains the day-specific `af`, `an` and `cases`, as well as
the group `ID` if applicable. For non-NULL grouping structures, this
output is structured per group.

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

# Calculate attributable fraction (forward perspective) for 2007

data2007_district <- dataset |>
  dplyr::filter(year == 2007)

# Calculate attributable risk
af2007_common <- attributable_risk(
  basis = basis.max,
  cases = data2007_district$dtot,
  model = model_max_common,
  cen = 25,
  x = data2007_district$Tmax,
  group = data2007_district$District, # Calculate af per district
  dir = "forw",
  tot = TRUE,
  samples = TRUE  # Samples from posterior distribution
)

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

# Calculate attributable fraction (backward perspective) for 2007

af2007_joint <- attributable_risk(
  basis = basis.max_min,
  cases = data2007_district$dtot,
  model = model_joint_common,
  cen = c(25,15),
  x = data2007_district$Tmax,
  y = data2007_district$Tmin,
  group = data2007_district$District, # Calculate af per district
  dir = "back",
  tot = TRUE,
  samples = TRUE  # Samples from posterior distribution
)


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


# Calculate attributable fraction (backward perspective) for 2007

af2007_spat <- attributable_risk(
  basis = basis.max_unpen,
  cases = data2007_district$dtot,
  model = model_max_spat,
  cen = 25,
  x = data2007_district$Tmax,
  by = data2007_district$District, # Spatially varying DLNM by district
  dir = "back",
  tot = TRUE,
  samples = TRUE  # Samples from posterior distribution
)


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

# Calculate attributable fraction (forward perspective) for 2007

af2007_conf <- attributable_risk(
  basis = basis.max_unpen,
  cases = data2007_district$dtot,
  model = model_max_conf,
  cen = 25,
  x = data2007_district$Tmax,
  by = data2007_district$depriv, # By-value for every district
  group = data2007_district$District, # Group output by district
  dir = "forw",
  tot = TRUE,
  samples = TRUE  # Samples from posterior distribution
)
} # }
```
