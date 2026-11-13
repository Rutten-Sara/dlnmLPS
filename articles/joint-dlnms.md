# Bivariate (joint) DLNMs

``` r

library(dlnmLPS)
```

## Introduction

Traditionally, distributed lag non-linear models (DLNMs) are used to
model the effect of one, possibly delayed and non-linear exposure on an
outcome. However, interest might be in the joint effect of two
(bivariate) exposures. For example, the effect of maximum and minimum
temperature on mortality might be synergistic, requiring a joint
modelling strategy.

`dlnmLPS` provides methods for modelling the delayed and non-linear
effect of two exposures jointly (Rutten et al. 2024). In this vignette,
we demonstrate how to fit and visualize these bivariate DLNMs using the
`dlnmLPS` framework.

The basic functionality of `dlnmLPS` is introduced in
[`vignette("introduction_dlnmlps")`](https://rutten-sara.github.io/dlnmLPS/articles/introduction_dlnmlps.md).

## Joint DLNMs

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
district. For this example, we will investigate how maximum and minimum
temperatures jointly impact mortality. The outcome and exposure variable
can then be inspected:

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

### Fitting an additive DLNM

While fitting a joint DLNM requires an extension of the basic
functionality of DLNMs, fitting a model with an additive DLNM of maximum
and minimum temperatures is straightforward.

We include month as natural splines with 3 degrees of freedom, year as a
random walk of order one, day of the week as a factor and District as an
(independent) random effect. Fitting a penalized additive DLNM, using 8
degrees of freedom and a penalty of order 2, can be done as follows:

``` r

# Fit a penalized additive DLNM
group <- factor(paste0(dataset$District, "-", dataset$year))

basis.max_pen <- dlnm::crossbasis(
  dataset$Tmax,
  argvar = list(fun = "ps", df = 8),
  arglag = list(fun = "ps", df = 8),
  lag = 14,
  group = group
)

basis.min_pen <- dlnm::crossbasis(
  dataset$Tmin,
  argvar = list(fun = "ps", df = 8),
  arglag = list(fun = "ps", df = 8),
  lag = 14,
  group = group
)

model_additive <- dlnmLPS(dtot ~ factor(dow) + 
                                sm(basis.max_pen, arg = list(pen = 2)) +
                                sm(basis.min_pen, arg = list(pen = 2)) +
                                sm(month, arg = list(type = "ns", df = 3)) +
                                sm(year, arg = list(type = "rw1")) +
                                sm(District, arg = list(type = "ind")),
                           data = dataset,
                           offset = log(pop),
                           model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                           family = "Poisson")
```

### Fitting a joint DLNM

Fitting a joint DLNM is very similar to fitting a traditional DLNM.
However, since the `dlnm` package does not directly provide the
functionality to specify a joint crossbasis, the function
[`joint.crossbasis()`](https://rutten-sara.github.io/dlnmLPS/reference/joint.crossbasis.md)
of the `dlnmLPS` package should be used instead. This function is very
similar to the
[`dlnm::crossbasis()`](https://rdrr.io/pkg/dlnm/man/crossbasis.html)
function, requiring arguments for both exposure variables as well as the
lag. The data can be supplied as a matrix for which the first and second
column contain the timeseries of the first and second exposure
respectively.

Fitting an unpenalized version of the model can be done through:

``` r

# Fit an unpenalized joint DLNM
knots_max <- quantile(dataset$Tmax, probs = c(0.10,0.50,0.90))
knots_min <- quantile(dataset$Tmin, probs = c(0.10,0.50,0.90))
knots_lag <- dlnm::logknots(0:14, nk = 2)

joint.basis <- joint.crossbasis(
                  cbind(dataset$Tmax, dataset$Tmin),
                  argvar1 = list(fun = "ns", knots = knots_max),
                  argvar2 = list(fun = "ns", knots = knots_min),
                  arglag = list(fun = "ns", knots = knots_lag),
                  lag = 14,
                  group = group
                  )

model_unpenalized <- dlnmLPS(dtot ~ factor(dow) + 
                                   sm(joint.basis) +
                                   sm(month, arg = list(type = "ns", df = 3)) +
                                   sm(year, arg = list(type = "rw1")) +
                                   sm(District, arg = list(type = "ind")),
                              data = dataset,
                              offset = log(pop),
                              model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                              family = "Poisson")
```

Similarly, a penalized version of the DLNM, using 8 degrees of freedom
and a second order penalty, can be fitted using:

``` r

# Fit a penalized joint DLNM
joint.basis_pen <- joint.crossbasis(
                    cbind(dataset$Tmax, dataset$Tmin),
                    argvar1 = list(fun = "ps", df = 8),
                    argvar2 = list(fun = "ps", df = 8),
                    arglag = list(fun = "ps", df = 8),
                    lag = 14,
                    group = group
                    )

model_penalized <- dlnmLPS(dtot ~ factor(dow) + 
                                 sm(joint.basis_pen, arg = list(pen = 2)) +
                                 sm(month, arg = list(type = "ns", df = 3)) +
                                 sm(year, arg = list(type = "rw1")) +
                                 sm(District, arg = list(type = "ind")),
                              data = dataset,
                              offset = log(pop),
                              model.comparison = list(DIC = TRUE, WAIC = TRUE, CPO = TRUE),
                              family = "Poisson")
```

To visualize this model, we need a grid of realistic maximum and minimum
temperature combinations. A prediction is made for every `at_x`-`at_y`
pair, rather than for cross-pairs. This also implies that the length of
both vectors should be equal. A centering value now consists of a vector
containing the reference value for both the first and the second
exposure. Predictions can then be made and exceedance probabilities can
be requested, using:

``` r

# Realistic maximum - minimum temperature data frame
df = data.frame()
k=1
for (i in seq(23,30,by=0.5)){
  lower_i = min(barcelona_temp$Tmin[barcelona_temp$Tmax<i+0.25 & barcelona_temp$Tmax>i-0.25])
  upper_i = max(barcelona_temp$Tmin[barcelona_temp$Tmax<i+0.25 & barcelona_temp$Tmax>i-0.25])
  for (j in seq(14,21,by=0.5)){
    if (j<upper_i & j>lower_i){
      df[k,1] = i
      df[k,2] = j
      k=k+1
    }
  }
}
colnames(df)=c("x","y")

# Predict RR
pred_penalized <- pred_crossbasis(
  basis = joint.basis_pen,
  model = model_penalized,
  at_x = df$x,
  at_y = df$y,
  cen = c(25,15),
  lag = c(0,14),
  CI = TRUE,
  exc.prob = TRUE
)
```

We can again extract specific predictions of interest using
[`extract_predictions()`](https://rutten-sara.github.io/dlnmLPS/reference/extract_predictions.md).
For example, if we want to extract the lag-specific RRs at the first
three `x` and `y` combinations, we can use:

``` r

extract_predictions(pred_penalized, slice_x = df$x[1:3], slice_y = df$y[1:3], ptype = "matRRfit")
#>     x    y lag       fit     lower     upper   exc.prob
#> 1  23 14.0   0 0.9852232 0.9672790 1.0035002 0.05621185
#> 2  23 14.5   0 0.9855114 0.9677840 1.0035635 0.05752740
#> 3  23 15.0   0 0.9877308 0.9688876 1.0069405 0.10452622
#> 4  23 14.0   1 0.9861360 0.9726874 0.9997705 0.02314584
#> 5  23 14.5   1 0.9874432 0.9737271 1.0013525 0.03831566
#> 6  23 15.0   1 0.9904178 0.9754364 1.0056293 0.10783565
#> 7  23 14.0   2 0.9881595 0.9759904 1.0004805 0.02978325
#> 8  23 14.5   2 0.9902222 0.9777724 1.0028306 0.06399125
#> 9  23 15.0   2 0.9936023 0.9800292 1.0073634 0.18021068
#> 10 23 14.0   3 0.9904948 0.9779017 1.0032500 0.07174049
#> 11 23 14.5   3 0.9928084 0.9801674 1.0056125 0.13481190
#> 12 23 15.0   3 0.9960596 0.9825564 1.0097483 0.28537646
#> 13 23 14.0   4 0.9926184 0.9808486 1.0045294 0.11172518
#> 14 23 14.5   4 0.9945830 0.9828344 1.0064721 0.18515179
#> 15 23 15.0   4 0.9971241 0.9846119 1.0097953 0.32743017
#> 16 23 14.0   5 0.9946184 0.9826611 1.0067213 0.19094221
#> 17 23 14.5   5 0.9958647 0.9840622 1.0078087 0.24786217
#> 18 23 15.0   5 0.9973745 0.9849670 1.0099384 0.34031316
#> 19 23 14.0   6 0.9966630 0.9845423 1.0089331 0.29618309
#> 20 23 14.5   6 0.9970953 0.9852384 1.0090949 0.31682387
#> 21 23 15.0   6 0.9975509 0.9852258 1.0100302 0.34953467
#> 22 23 14.0   7 0.9988407 0.9875366 1.0102742 0.42083816
#> 23 23 14.5   7 0.9985686 0.9875410 1.0097193 0.40020198
#> 24 23 15.0   7 0.9981546 0.9867159 1.0097260 0.37672788
#> 25 23 14.0   8 1.0011361 0.9895516 1.0128563 0.57582035
#> 26 23 14.5   8 1.0003858 0.9891639 1.0117350 0.52671343
#> 27 23 15.0   8 0.9993792 0.9878672 1.0110253 0.45816472
#> 28 23 14.0   9 1.0035039 0.9921742 1.0149631 0.72700594
#> 29 23 14.5   9 1.0026054 0.9916616 1.0136700 0.67891286
#> 30 23 15.0   9 1.0013622 0.9901785 1.0126722 0.59388802
#> 31 23 14.0  10 1.0055494 0.9949785 1.0162327 0.84763374
#> 32 23 14.5  10 1.0048183 0.9946035 1.0151380 0.82174095
#> 33 23 15.0  10 1.0036826 0.9932333 1.0142417 0.75439916
#> 34 23 14.0  11 1.0066039 0.9957681 1.0175576 0.88336392
#> 35 23 14.5  11 1.0062480 0.9958078 1.0167977 0.87909883
#> 36 23 15.0  11 1.0054794 0.9948762 1.0161957 0.84381310
#> 37 23 14.0  12 1.0061197 0.9956823 1.0166664 0.87424567
#> 38 23 14.5  12 1.0062321 0.9961856 1.0163800 0.88753147
#> 39 23 15.0  12 1.0059944 0.9958180 1.0162749 0.87535995
#> 40 23 14.0  13 1.0045321 0.9942363 1.0149346 0.80518142
#> 41 23 14.5  13 1.0050748 0.9951923 1.0150554 0.84232206
#> 42 23 15.0  13 1.0053751 0.9954504 1.0153988 0.85521776
#> 43 23 14.0  14 1.0027260 0.9917013 1.0138733 0.68531598
#> 44 23 14.5  14 1.0035239 0.9929970 1.0141624 0.74338216
#> 45 23 15.0  14 1.0041873 0.9937729 1.0147109 0.78394616
```

A contour plot and 3D plot of the overall cumulative RR can be made:

``` r

contour(pred_penalized, ptype = "overall", xlab = "maximum temperature", ylab = "minimum temperature")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-10-1.png)

``` r

plot(pred_penalized, ptype = "overall", ticktype="detailed",theta=230,
        xlab = "maximum", ylab = "minimum", zlab = "RR", 
        ltheta=200,phi=30,lphi=30,
        nticks = 4, shade = 0.75,r=sqrt(3),d=5,cex.axis=0.5, cex.lab=0.8,border=NA)
```

![](joint-dlnms_files/figure-html/unnamed-chunk-10-2.png)

The plots seem to suggest that both maximum and minimum temperature are
important when estimating RRs. We can for example also visualize the lag
relationship for specific exposure levels:

``` r

plot(pred_penalized, slice_x = rep(30, sum(pred_penalized$x == 30)), 
     slice_y = pred_penalized$y[pred_penalized$x == 30], ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-11-1.png)

or the exposures-outcome relationship for the immediate risk (lag = 0):

``` r

contour(pred_penalized, slice_lag = 0, xlab = "maximum temperature", ylab = "minimum temperature", 
        ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-12-1.png)

and risk after five days:

``` r

contour(pred_penalized, slice_lag = 5, xlab = "maximum temperature", ylab = "minimum temperature",
        ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-13-1.png)

The attributable fraction in 2007 can also be calculated for example. We
make use of the district-specific maximum and minimum temperature
timeseries and of the backward perspective:

``` r

# Calculate exposure and outcome timeseries
data2007_district <- dataset |>
  dplyr::filter(year == 2007) |>
  dplyr::arrange(District, date)

# Calculate attributable risk
af2007_district <- attributable_risk(
  basis = joint.basis_pen, 
  cases = data2007_district$dtot, 
  model = model_penalized,
  cen = c(25,15),
  x = data2007_district$Tmax,
  y = data2007_district$Tmin,
  group = data2007_district$District,
  dir = "back",
  tot = TRUE,
  samples = TRUE)
```

These results could then be visualized as:

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

# Visualize results
ggplot2::ggplot(af_ci, ggplot2::aes(x = estimate, y = ID)) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), orientation = "y") +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Attributable fraction", y = "District") +
  ggplot2::theme_minimal()
```

![](joint-dlnms_files/figure-html/unnamed-chunk-15-1.png)

### Fitting a spatial and effect modified joint DLNM

Just as for traditional DLNMs, a joint DLNM can be spatially varying or
subject to effect modification (see
[`vignette("spatially-varying-dlnms")`](https://rutten-sara.github.io/dlnmLPS/articles/spatially-varying-dlnms.md)
and
[`vignette("effect-modified-dlnms")`](https://rutten-sara.github.io/dlnmLPS/articles/effect-modified-dlnms.md)
for more information). The syntax remains exactly the same as for the
univariate DLNM.

## Summary

In this vignette, we introduced the joint modelling functionality of
`dlnmLPS`. We demonstrated how to fit and visualize these models. For
more information, see also
[`vignette("introduction_dlnmlps")`](https://rutten-sara.github.io/dlnmLPS/articles/introduction_dlnmlps.md),
[`vignette("spatially-varying-dlnms")`](https://rutten-sara.github.io/dlnmLPS/articles/spatially-varying-dlnms.md)
and
[`vignette("effect-modified-dlnms")`](https://rutten-sara.github.io/dlnmLPS/articles/effect-modified-dlnms.md).

## References

Rutten, S., E. Duarte, T. Neyens, D. Lauwaet, and C. Faes. 2024. “A
Penalized Distributed-Lag Non-Linear Model for Modeling the Joint
Delayed Effect of Two Predictors: Impact of Minimum and Maximum
Temperature on Mortality.” *medRxiv*, ahead of print.
<https://doi.org/10.1101/2024.11.29.24318041>.
