# dlnmLPS

![Version](https://img.shields.io/badge/Version-0.1.0-red)![Languages](https://img.shields.io/badge/Languages-R%2C%20C%2B%2B-informational)![Lifecycle](https://img.shields.io/badge/lifecycle-Experimental-orange)![CodeSize](https://img.shields.io/github/languages/code-size/Rutten-Sara/dlnmLPS?color=orange&label=Code%20size&style=plastic)![License](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue)![UserStars](https://img.shields.io/github/stars/Rutten-Sara/dlnmLPS)

**dlnmLPS** is an R package for fitting **d**istributed **l**ag
**n**on-linear **m**odels within a Bayesian framework using
**L**aplacian-**P**-**S**plines. It provides a flexible approach for
modelling complex exposure-lag-response relationships, including
relationships that vary spatially or according to pre-specified
continuous or categorical modifiers, as well as exposure-lag-response
relationships for two joint exposures. In addition to
exposure-lag-response relationships, `dlnmLPS` allows other sources of
heterogeneity to be incorporated, including smooth covariate effects,
random effects and spatial and temporal structures.

More specifically, `dlnmLPS` can be used to model:

- **distributed lag non-linear associations**, allowing nonlinear
  exposure-response relationships over multiple time lags. DLNMs can be
  unpenalized (natural splines or B-splines) or penalized (P-splines);

- **random effects**, including independent random effects and spatially
  structured random effects (ICAR, BYM or Leroux). These random effects
  can possibly vary over time (according to an independent or random
  walk structure);

- **temporal random effects**, including first- and second-order
  random-walk structures. These temporal effects can possibly depend on
  location (according to an independent, ICAR or Leroux structure);

- **smooth exposure effects**, using natural splines (`ns`), B-splines
  (`bs`), or P-splines (`ps`); and

- **more complex DLNM extensions**, including effect modification,
  spatially varying DLNMs, and models with bivariate exposures.

The different components can be combined to construct models that
account simultaneously for nonlinear exposure-response relationships,
distributed lag effects, and spatial or temporal dependence. Supported
outcome distributions are Poisson, negative binomial, binomial and
Gaussian.

## Documentation

Detailed documentation for `dlnmLPS` is available through the tabs on
this website:

- The **Reference** tab provides an overview of the functions available
  in `dlnmLPS`, including their arguments and return values.

- The **Articles** tab provides more detailed explanations of the
  package functionality, with examples illustrating how to specify and
  fit different models.

## Installation

The development version of `dlnmLPS` can be installed from the GitHub
repository using:

``` r

# install.packages("remotes")
remotes::install_github("Rutten-Sara/dlnmLPS")
```

Access to the repository is required to install the package.

To install `dlnmLPS` together with its vignettes, containing short
tutorials, use:

``` r

# install.packages("remotes")
remotes::install_github("Rutten-Sara/dlnmLPS", build_vignettes = TRUE)
```

Afterwards, the vignettes can be requested using:

``` r

browseVignettes("dlnmLPS")
```

## Citation

Associated literature:

1.  Rutten, S., Duarte, E., Neyens, T., Lauwaet, D., & Faes, C. (2024).
    A penalized distributed lag non-linear model for modeling the joint
    delayed effect of two predictors: impact of minimum and maximum
    temperature on mortality. *medRxiv*.
    <https://doi.org/10.1101/2024.11.29.24318041>

2.  Rutten, S., Sumalinab, B., Gressani, O., Neyens, T., Duarte, E.,
    Hens, N., & Faes, C. (2026). Penalized distributed lag non-linear
    models for small area data using Laplacian-P-splines. *Statistics
    and Computing*, 36(1), 38.
    <https://doi.org/10.1007/s11222-025-10790-9>

3.  Rutten, S., Neyens, T., Duarte, E., Gasparrini, A., & Faes, C.
    (2026). Distributed lag non-linear models with spatial effect
    modification using Laplacian P-splines. *Environmental and
    Ecological Statistics*. <https://doi.org/10.1007/s10651-026-00739-0>

4.  Rutten, S., Neyens, T., Duarte, E., Gasparrini, A., & Faes, C.
    (2026). Spatially varying distributed lag non-linear models using
    Laplacian P-splines. *arXiv preprint arXiv:2604.09012*.
    <https://doi.org/10.48550/arXiv.2604.09012>

## Package version

This is version 0.1.0 (2026-09-01).

`dlnmLPS` is currently under active development. This is an early
development version, and the package interface may change as
functionality is developed and tested.

## License

`dlnmLPS`: a fast and flexible Bayesian tool for Distributed Lag
Non-Linear Models. Copyright (C) 2026 Sara Rutten. This program is free
software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any
later version.

`dlnmLPS` is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY. See the `LICENSE` file for details.
