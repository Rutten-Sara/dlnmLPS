# Effect-modified DLNMs

``` r

library(dlnmLPS)
```

## Introduction

Classical distributed lag non-linear models (DLNMs) assume a common
exposure-lag-response relationship across geographical areas. However,
this relationship might be altered by an effect modifier that differs
between spatial units. For example, the relationship between temperature
and mortality may differ between regions according to the deprivation
index.

`dlnmLPS` provides methods for fitting effect-modified distributed lag
non-linear models, allowing the exposure-lag-response relationship to
vary according to a continuous or categorical covariate (Rutten et al.
2026). This framework allows for flexible interactions between the DLNM
and a possibly non-linear smooth covariate.

In this vignette, we demonstrate how to fit effect-modified DLNMs using
`dlnmLPS` and show how to predict and visualize the
exposure-lag-response relationship.

The basic functionality of `dlnmLPS` is introduced in
[`vignette("introduction_dlnmlps")`](https://rutten-sara.github.io/dlnmLPS/articles/introduction_dlnmlps.md).

## Effect-modified DLNMs

### Example data

For illustration, we use the following simulated example data:

``` r

# Example data
data(barcelona_temp)
```

This simulated dataset contains daily minimum, maximum and mean
temperature in all 10 districts of Barcelona during the summer period (1
June–30 September) between 2007 and 2011. The population size,
deprivation index and daily number of deaths are included for every
district. For this example, we will investigate how the effect of
maximum temperature on mortality varies spatially according to the
deprivation index. The outcome and exposure variable can then be
inspected:

``` r

# Inspect the data
head(barcelona_temp)
#>         date     Tmax     Tmin    Tmean District      pop   depriv dtot
#> 1 2007-06-01 25.85997 17.93684 21.89841       01 659845.8 1.188957   26
#> 2 2007-06-02 25.71266 17.21590 21.46428       01 659845.8 1.188957   27
#> 3 2007-06-03 25.33624 16.95453 21.14538       01 659845.8 1.188957   23
#> 4 2007-06-04 25.47256 15.42090 20.44673       01 659845.8 1.188957   23
#> 5 2007-06-05 23.89031 16.21367 20.05199       01 659845.8 1.188957   26
#> 6 2007-06-06 24.08012 16.47465 20.27739       01 659845.8 1.188957   18
```

We can extract some additional information from the date such as the day
of the week, year and month.

``` r

# Extract date information
dataset <- barcelona_temp |>
  dplyr::mutate(
    dow = weekdays(date),
    year = lubridate::year(date),
    month = lubridate::month(date)
  )
```

### Fitting an effect-modified DLNM

The `dlnmLPS` package provides a lot of flexibility in modelling effect
modification. In this example, we will fit four different models:

- **Categorical effect modification**, a penalized DLNM modified by a
  categorical deprivation index;

- **Linear effect modification**, a penalized DLNM modified by a
  continuous deprivation index, using a linear interaction;

- **Smooth unpenalized effect modification**, a penalized DLNM modified
  by a continuous deprivation index, expressed by natural splines using
  2 knots (on the 10th and 90th percentiles);

- **Smooth penalized effect modification**, a penalized DLNM modified by
  a continuous deprivation index, expressed by penalized splines using 5
  df.

We include month as a natural spline with 3 degrees of freedom, year as
a random walk of order one, day of the week as a factor and district as
an independent random effect. We use a penalized DLNM with 8 degrees of
freedom and penalty of order 2.

For the first model, we split the deprivation index into a categorical
variable with three levels (below the 25th percentile, 25th-75th
percentile, above the 75th percentile). We specify the modifying
covariate through the option `by` and include `type = "factor"` in the
argument `arg.by`. We also include the categorical deprivation index as
a main effect in the model. The model can be fitted as:

``` r

# Fit a categorical modified DLNM
group <- factor(paste0(dataset$District, "-", dataset$year))

basis.max_pen <- dlnm::crossbasis(
  dataset$Tmax,
  argvar = list(fun = "ps", df = 8),
  arglag = list(fun = "ps", df = 8),
  lag = 14,
  group = group
)

dataset$depriv_cat <- ifelse(dataset$depriv <= quantile(dataset$depriv,0.25), "low",
                             ifelse(dataset$depriv >= quantile(dataset$depriv,0.75), "high", "medium"))

model1 <- dlnmLPS(dtot ~ factor(dow) + factor(depriv_cat) +
                        sm(basis.max_pen, arg = list(pen = 2),
                           by = depriv_cat, arg.by = list(type = "factor")) +
                        sm(month, arg = list(type = "ns", df = 3)) +
                        sm(year, arg = list(type = "rw1")) +
                        sm(District, arg = list(type = "ind")),
                   data = dataset,
                   offset = log(pop),
                   model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                   family = "Poisson")
```

The second model can similarly be fitted by specifying `type = "linear"`
instead of `type = "factor"`. Deprivation index is now also included as
a linear main effect:

``` r

# Fit a linear modified DLNM
model2 <- dlnmLPS(dtot ~ factor(dow) + depriv +
                        sm(basis.max_pen, arg = list(pen = 2),
                           by = depriv, arg.by = list(type = "linear")) +
                        sm(month, arg = list(type = "ns", df = 3)) +
                        sm(year, arg = list(type = "rw1")) +
                        sm(District, arg = list(type = "ind")),
                   data = dataset,
                   offset = log(pop),
                   model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                   family = "Poisson")
```

For the third model, we include deprivation index as a smooth effect,
using natural splines with two knots. The model can be fitted as:

``` r

# Fit an unpenalized smooth modified DLNM
knots_depriv <- quantile(dataset$depriv, probs = c(0.10,0.90))

model3 <- dlnmLPS(dtot ~ factor(dow) + sm(depriv, arg = list(type = "ns", knots = knots_depriv)) +
                        sm(basis.max_pen, arg = list(pen = 2),
                           by = depriv, arg.by = list(type = "ns", knots = knots_depriv)) +
                        sm(month, arg = list(type = "ns", df = 3)) +
                        sm(year, arg = list(type = "rw1")) +
                        sm(District, arg = list(type = "ind")),
                   data = dataset,
                   offset = log(pop),
                   model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                   family = "Poisson")
```

Lastly, in the fourth model, we include an interaction with deprivation
index using penalized splines:

``` r

# Fit a penalized smooth modified DLNM
model4 <- dlnmLPS(dtot ~ factor(dow) + sm(depriv, arg = list(type = "ps", df = 5, pen = 2)) +
                        sm(basis.max_pen, arg = list(pen = 2),
                           by = depriv, arg.by = list(type = "ps", df = 5, pen = 2)) +
                        sm(month, arg = list(type = "ns", df = 3)) +
                        sm(year, arg = list(type = "rw1")) +
                        sm(District, arg = list(type = "ind")),
                   data = dataset,
                   offset = log(pop),
                   model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                   family = "Poisson")
```

We can compare the performance of the four models using model comparison
criteria DIC, WAIC and CPO:

``` r

# Compare the estimated associations
model1$DIC; model1$WAIC; -sum(log(model1$CPO))
#> [1] 30411.85
#> [1] 30415.3
#> [1] 15207.56
model2$DIC; model2$WAIC; -sum(log(model2$CPO))
#> [1] 30389.48
#> [1] 30392.5
#> [1] 15196.13
model3$DIC; model3$WAIC; -sum(log(model3$CPO))
#> [1] 30389.52
#> [1] 30392.24
#> [1] 15196.01
model4$DIC; model4$WAIC; -sum(log(model4$CPO))
#> [1] 30378.99
#> [1] 30382.03
#> [1] 15190.93
```

The fourth model is preferred based on all model selection criteria.
Using this model, we predict the exposure-lag-response curve for three
different levels of the deprivation index (2.5th, 50th and 97.5th
percentiles) at temperatures between 22 and 30 degrees Celsius, relative
to a centering value of 25 degrees. We also request the probability that
the logarithm of the relative risk exceeds 0:

``` r

pred_penalized <- pred_crossbasis(
  basis = basis.max_pen,
  model = model4,
  at_x = seq(22,30,length = 50),
  cen = 25,
  lag = c(0,14),
  CI = TRUE,
  by = round(quantile(dataset$depriv, probs = c(0.025,0.50,0.975)),2),
  exc.prob = TRUE
)
```

Using the
[`extract_predictions()`](https://rutten-sara.github.io/dlnmLPS/reference/extract_predictions.md)
function, we can easily extract (part of) the requested predictions. For
example, if interest is in the lag-specific RR for lag 0 and lag 2 at 22
and 30 degrees and the 2.5th and 97.5th percentiles of deprivation, we
can use:

``` r

extract_predictions(pred_penalized, slice_x = c(22,30), slice_lag = c(0,2), 
                    slice_by = round(quantile(dataset$depriv, probs = c(0.025,0.975)),2), 
                    ptype = "matRRfit")
#>    x lag       fit     lower    upper   exc.prob    by
#> 1 22   0 0.9442494 0.8760118 1.017802 0.06695034 -2.09
#> 2 30   0 1.1271132 1.0665366 1.191130 0.99998909 -2.09
#> 3 22   2 0.9529744 0.9034242 1.005242 0.03852666 -2.09
#> 4 30   2 1.0263029 0.9943068 1.059329 0.94593421 -2.09
#> 5 22   0 0.9615895 0.8950778 1.033043 0.14208096  1.45
#> 6 30   0 1.1422374 1.0804183 1.207594 0.99999860  1.45
#> 7 22   2 1.0057789 0.9596735 1.054099 0.59509832  1.45
#> 8 30   2 1.0925640 1.0596552 1.126495 0.99999999  1.45
```

To plot the exposure-response relationship at lag 0, 3 and 8 for the
2.5th, 50th and 97.5th percentiles of the deprivation index, we can use:

``` r

plot(pred_penalized, slice_lag = c(0,3,8), ptype = "slices", CI = TRUE, 
     xlab = "temperature", ylab = "RR")
```

![](effect-modified-dlnms_files/figure-html/unnamed-chunk-12-1.png)

Based on this plot, the differences in the lag-specific RR between
deprivation levels seem small. The plot title indicates the level of
deprivation. The option `multiple.plots = FALSE` allows visualizing the
results for different levels of the deprivation index on one plot. Note
that this option is only available for a single (scalar) slice value or
overall RR-curve. Furthermore, the option `slice_by` can be used to
specify only a subset of the predicted `by`-values. For example, if we
only want to plot the 2.5th and 97.5th deprivation percentiles, we can
use:

``` r

plot(pred_penalized, slice_x = 30, 
     slice_by = round(quantile(dataset$depriv,probs = c(0.025,0.975)),2), 
     ptype = "slices", CI = TRUE, xlab = "lag", ylab = "RR", multiple.plots = FALSE)
```

![](effect-modified-dlnms_files/figure-html/unnamed-chunk-13-1.png)

``` r

plot(pred_penalized, slice_by = round(quantile(dataset$depriv, probs = c(0.025,0.975)),2), 
     ptype = "overall", CI = TRUE, xlab = "temperature", ylab = "RR",  
     multiple.plots = FALSE)
```

![](effect-modified-dlnms_files/figure-html/unnamed-chunk-13-2.png)

At high temperatures, there is some evidence that the overall RR is
lower for districts with a low deprivation index. A contour plot of the
exposure-lag-response relationship can also be made for different
deprivation levels:

``` r

contour(pred_penalized, ptype = "cross", xlab = "lag", ylab = "temperature")
```

![](effect-modified-dlnms_files/figure-html/unnamed-chunk-14-1.png)![](effect-modified-dlnms_files/figure-html/unnamed-chunk-14-2.png)![](effect-modified-dlnms_files/figure-html/unnamed-chunk-14-3.png)

We could also calculate the attributable fraction in 2007 for example.
Using the district-specific maximum temperatures, deprivation index and
deaths, and the backward perspective, we get:

``` r

# Calculate exposure and outcome time series
data2007_district <- dataset |>
  dplyr::filter(year == 2007) |>
  dplyr::arrange(District, date)

# Calculate attributable risk
af2007_district <- attributable_risk(
  basis = basis.max_pen, 
  cases = data2007_district$dtot, 
  model = model4,
  cen = 25,
  x = data2007_district$Tmax,
  by = data2007_district$depriv,
  group = data2007_district$District,
  dir = "back",
  tot = TRUE,
  samples = TRUE)
```

We could for example also calculate a counterfactual attributable
fraction, in which the deprivation index in every area is equal to the
2.5th percentile of the deprivation index.

``` r

# Calculate counterfactual attributable risk
depriv_lower <- round(quantile(data2007_district$depriv, probs = 0.025, na.rm = TRUE),2)
af2007_counter <- attributable_risk(
  basis = basis.max_pen, 
  cases = data2007_district$dtot, 
  model = model4,
  cen = 25,
  x = data2007_district$Tmax,
  by = rep(depriv_lower, length(data2007_district$depriv)),
  group = data2007_district$District,
  dir = "back",
  tot = TRUE,
  samples = TRUE)
```

Both results could be visualized as:

``` r

# Extract results
af_ci <- dplyr::bind_rows(lapply(seq_along(af2007_district$samples), function(i) {
  x <- af2007_district$samples[[i]]
  data.frame(
    ID = x$ID[1],
    estimate = af2007_district$fit[[i]]$af,
    lower = quantile(x$af, 0.025, na.rm = TRUE),
    upper = quantile(x$af, 0.975, na.rm = TRUE)
  )
}))

counter_ci <- dplyr::bind_rows(lapply(seq_along(af2007_counter$samples), function(i) {
  x <- af2007_counter$samples[[i]]
  data.frame(
    ID = x$ID[1],
    estimate = af2007_counter$fit[[i]]$af,
    lower = quantile(x$af, 0.025, na.rm = TRUE),
    upper = quantile(x$af, 0.975, na.rm = TRUE)
  )
}))

af_plot_df <- dplyr::bind_rows(
  dplyr::mutate(af_ci, scenario = "Observed"),
  dplyr::mutate(counter_ci, scenario = "Counterfactual")
)


# Visualize results
ggplot2::ggplot(af_plot_df, ggplot2::aes(x = estimate, y = ID, colour = scenario)) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), orientation = "y", 
                         position = ggplot2::position_dodge(0.5)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::scale_colour_manual(values = c("Observed" = "#0072B2","Counterfactual" = "#D55E00")) +
  ggplot2::labs(x = "Attributable fraction", y = "District",colour = "Scenario") +
  ggplot2::theme_minimal()
```

![](effect-modified-dlnms_files/figure-html/unnamed-chunk-17-1.png)

It can be seen that a lower deprivation index would result in a lower
attributable fraction for some districts.

## Summary

In this vignette, we introduced the effect-modification functionality of
`dlnmLPS`. We demonstrated how to fit four different models, differing
in the flexibility of the effect modification. We used penalized splines
to showcase the functionality of the package, but similar models can be
fitted with unpenalized splines. See also
[`vignette("introduction_dlnmlps")`](https://rutten-sara.github.io/dlnmLPS/articles/introduction_dlnmlps.md)
for more information.

## References

Rutten, Sara, Thomas Neyens, Elisa Duarte, Antonio Gasparrini, and
Christel Faes. 2026. “Distributed Lag Non-Linear Models with Spatial
Effect Modification Using Laplacian p-Splines.” *Environmental and
Ecological Statistics*, ahead of print, June.
<https://doi.org/10.1007/s10651-026-00739-0>.
