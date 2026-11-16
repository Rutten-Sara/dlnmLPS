# Bivariate (joint) DLNMs

``` r

library(dlnmLPS)
set.seed(123)
```

## Introduction

Traditionally, distributed lag non-linear models (DLNMs) are used to
model the possibly delayed and non-linear effect of one exposure on an
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
    dow = lubridate::wday(date),
    year = lubridate::year(date),
    month = lubridate::month(date)
  )
```

### Fitting an additive DLNM

While fitting a joint DLNM requires an extension of the basic
functionality of DLNMs, fitting a model with an additive DLNM of maximum
and minimum temperatures is straightforward.

We include month as a natural spline with 3 degrees of freedom, year as
a random walk of order one, day of the week as a factor and district as
an (independent) random effect. Fitting a penalized additive DLNM, using
8 degrees of freedom and a penalty of order 2, can be done as follows:

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
columns contain the time series of the first and second exposure,
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
and a second-order penalty, can be fitted using:

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
#> 1  23 14.0   0 0.9852169 0.9672663 1.0035005 0.05620042
#> 2  23 14.5   0 0.9855044 0.9677709 1.0035628 0.05750335
#> 3  23 15.0   0 0.9877244 0.9688750 1.0069405 0.10448385
#> 4  23 14.0   1 0.9861312 0.9726787 0.9997699 0.02314140
#> 5  23 14.5   1 0.9874389 0.9737188 1.0013523 0.03830794
#> 6  23 15.0   1 0.9904150 0.9754296 1.0056307 0.10783262
#> 7  23 14.0   2 0.9881577 0.9759850 1.0004822 0.02980018
#> 8  23 14.5   2 0.9902215 0.9777683 1.0028334 0.06403024
#> 9  23 15.0   2 0.9936038 0.9800272 1.0073684 0.18032519
#> 10 23 14.0   3 0.9904954 0.9778984 1.0032547 0.07181751
#> 11 23 14.5   3 0.9928105 0.9801657 1.0056184 0.13495216
#> 12 23 15.0   3 0.9960638 0.9825570 1.0097563 0.28563789
#> 13 23 14.0   4 0.9926201 0.9808474 1.0045341 0.11183759
#> 14 23 14.5   4 0.9945859 0.9828345 1.0064778 0.18533758
#> 15 23 15.0   4 0.9971288 0.9846139 1.0098028 0.32773112
#> 16 23 14.0   5 0.9946202 0.9826605 1.0067254 0.19106609
#> 17 23 14.5   5 0.9958671 0.9840624 1.0078135 0.24803088
#> 18 23 15.0   5 0.9973782 0.9849685 1.0099442 0.34055017
#> 19 23 14.0   6 0.9966642 0.9845411 1.0089365 0.29628089
#> 20 23 14.5   6 0.9970967 0.9852377 1.0090985 0.31693983
#> 21 23 15.0   6 0.9975531 0.9852260 1.0100344 0.34968463
#> 22 23 14.0   7 0.9988410 0.9875355 1.0102759 0.42086726
#> 23 23 14.5   7 0.9985690 0.9875401 1.0097212 0.40024526
#> 24 23 15.0   7 0.9981555 0.9867156 1.0097281 0.37679847
#> 25 23 14.0   8 1.0011356 0.9895500 1.0128569 0.57578141
#> 26 23 14.5   8 1.0003855 0.9891625 1.0117358 0.52668956
#> 27 23 15.0   8 0.9993792 0.9878663 1.0110262 0.45816866
#> 28 23 14.0   9 1.0035031 0.9921724 1.0149632 0.72693959
#> 29 23 14.5   9 1.0026048 0.9916600 1.0136703 0.67885776
#> 30 23 15.0   9 1.0013618 0.9901774 1.0126726 0.59385587
#> 31 23 14.0  10 1.0055486 0.9949773 1.0162322 0.84758812
#> 32 23 14.5  10 1.0048177 0.9946026 1.0151377 0.82170186
#> 33 23 15.0  10 1.0036821 0.9932327 1.0142415 0.75437138
#> 34 23 14.0  11 1.0066032 0.9957672 1.0175571 0.88333570
#> 35 23 14.5  11 1.0062475 0.9958071 1.0167973 0.87907670
#> 36 23 15.0  11 1.0054791 0.9948758 1.0161953 0.84379814
#> 37 23 14.0  12 1.0061191 0.9956818 1.0166657 0.87422480
#> 38 23 14.5  12 1.0062316 0.9961852 1.0163794 0.88751509
#> 39 23 15.0  12 1.0059940 0.9958178 1.0162741 0.87534828
#> 40 23 14.0  13 1.0045316 0.9942364 1.0149333 0.80516727
#> 41 23 14.5  13 1.0050742 0.9951923 1.0150542 0.84230799
#> 42 23 15.0  13 1.0053744 0.9954505 1.0153973 0.85520545
#> 43 23 14.0  14 1.0027255 0.9917018 1.0138719 0.68530036
#> 44 23 14.5  14 1.0035233 0.9929972 1.0141609 0.74336013
#> 45 23 15.0  14 1.0041865 0.9937731 1.0147091 0.78392187
```

A contour plot and 3D plot of the overall cumulative RR can be made:

``` r

contour(pred_penalized, ptype = "overall", xlab = "maximum temperature", ylab = "minimum temperature")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-10-1.png)

``` r

plot(pred_penalized, ptype = "overall", ticktype="detailed",theta=230,
        xlab = "maximum temperature", ylab = "minimum temperature", zlab = "RR", 
        ltheta=200,phi=30,lphi=30,
        nticks = 4, shade = 0.75,r=sqrt(3),d=5,cex.axis=0.5, cex.lab=0.8,border=NA)
```

![](joint-dlnms_files/figure-html/unnamed-chunk-10-2.png)

The plots seem to suggest that both maximum and minimum temperatures are
important when estimating RRs. We can for example also visualize the lag
relationship for specific exposure levels:

``` r

plot(pred_penalized, slice_x = rep(30, sum(pred_penalized$x == 30)), 
     slice_y = pred_penalized$y[pred_penalized$x == 30], ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-11-1.png)

or the exposure-response relationship for the immediate risk (lag = 0):

``` r

contour(pred_penalized, slice_lag = 0, xlab = "maximum temperature", ylab = "minimum temperature", 
        ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-12-1.png)

and the risk at lag 5:

``` r

contour(pred_penalized, slice_lag = 5, xlab = "maximum temperature", ylab = "minimum temperature",
        ptype = "slices")
```

![](joint-dlnms_files/figure-html/unnamed-chunk-13-1.png)

The attributable fraction in 2007 can also be calculated for example. We
make use of the district-specific maximum and minimum temperature time
series and of the backward perspective:

``` r

# Calculate exposure and outcome time series
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

### Spatially varying and effect-modified joint DLNMs

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

Rutten, Sara, Elisa Duarte, Thomas Neyens, Dirk Lauwaet, and Christel
Faes. 2024. “A Penalized Distributed-Lag Non-Linear Model for Modeling
the Joint Delayed Effect of Two Predictors: Impact of Minimum and
Maximum Temperature on Mortality.” In *medRxiv*.
<https://doi.org/10.1101/2024.11.29.24318041>.
