
# DLNMLPS

<!-- badges: start -->
![Version](https://img.shields.io/badge/Version-0.1.0-red)
![Languages](https://img.shields.io/badge/Languages-R%2C%20C%2B%2B-informational)
![Lifecycle](https://img.shields.io/badge/lifecycle-Experimental-orange)
![CodeSize](https://img.shields.io/github/languages/code-size/Rutten-Sara/DLNMLPS?color=orange&label=Code%20size&style=plastic)
![License](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue)
![UserStars](https://img.shields.io/github/stars/Rutten-Sara/DLNMLPS?style=social)
<!-- badges: end -->

**DLNMLPS** is an R package for fitting **D**istributed **L**ag **N**on-linear
**M**odels within a Bayesian framework using **L**aplacian-**P**-**S**plines. It
provides a flexible approach for modelling complex
exposure-lag-response relationships, including relationships that
vary spatially or according to other continuous or categorical
modifiers. The package also allows for the inclusion of structured
spatial and temporal effects.

## Installation

The development version of DLNMLPS can be installed from the
private GitHub repository using:

```r
# install.packages("remotes")
remotes::install_github("Rutten-Sara/DLNMLPS")
```
Access to the repository is required to install the package.

## Citation

Associated literature:
1. Rutten, S., Duarte, E., Neyens, T., Lauwaet, D., & Faes, C.
   (2024). A penalized distributed-lag non-linear model for modeling
   the joint delayed effect of two predictors: impact of minimum and
   maximum temperature on mortality. *medRxiv*.
   <https://doi.org/10.1101/2024.11.29.24318041>
2. Rutten, S., Sumalinab, B., Gressani, O., Neyens, T., Duarte, E., 
   Hens, N., & Faes, C. (2026). Penalized distributed lag non-linear 
   models for small area data using Laplacian-P-splines. 
   *Statistics and Computing*, 36(1), 38.
   <https://doi.org/10.1007/s11222-025-10790-9>
3. Rutten, S., Neyens, T., Duarte, E., Gasparrini, A., & Faes, C. (2026). 
   Distributed lag non-linear models with spatial effect modification 
   using Laplacian P-splines. *Environmental and Ecological Statistics*.
  <https://doi.org/10.1007/s10651-026-00739-0>
4. Rutten, S., Neyens, T., Duarte, E., Gasparrini, A., & Faes, C.
   (2026). Spatially varying distributed lag non-linear models using
   Laplacian P-splines. *arXiv preprint arXiv:2604.09012*.
   <https://doi.org/10.48550/arXiv.2604.09012>

To cite the DLNMLPS package in publications, use:
    Rutten, S., Neyens, T., Duarte, E. & Faes, C. (2026).
    DLNMLPS: Flexible Bayesian Distributed Lag Non-Linear Models.[Computer Software].
   
## Package version

This is version 0.1.0 (2026-09-01). 

DLNMLPS is currently under active development. This is an early
development version, and the package interface may change
as functionality is developed and tested.

## License
DLNMLPS: a fast and flexible Bayesian tool for Distributed Lag Non-Linear Models.
Copyright (C) 2026 Sara Rutten. 
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

DLNMLPS is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY. See the `LICENSE` file for details.


