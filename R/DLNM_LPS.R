#' Fit DLNM using Laplacian P-splines
#'
#' `DLNM_LPS` fits distributed lag non-linear models (DLNMs) using a
#' Laplace approximation framework. The method is primarily designed for
#' penalized spline representations, but also accommodates unpenalized spline
#' bases within the DLNM framework. The model allows inclusion of structured
#' components, such as spatial random effects and random walks. Smooth
#' covariate effects can be specified using natural splines, P-splines, or
#' B-splines.
#'
#' @param model A model formula. Smooth terms specified using [sm] can be
#'   included on the right-hand side to define random effects, random walks,
#'   DLNM cross-basis terms, or other smooth functions.
#' @param data A data frame containing the response variable and covariates
#'   specified in the model formula.
#' @param offset An optional model offset term.
#' @param model.comparison A list specifying optional model evaluation criteria to calculate.
#'   Possible elements are `DIC`, `WAIC`, and `CPO`, each taking a logical value
#'   (`TRUE` or `FALSE`). All criteria default to `FALSE`.
#' @param family The distribution family used for the response model. Currently
#'   supported options are `poisson` (default) and `NB` for the negative
#'   binomial distribution.
#' @param refit Optional number of additional model refinement iterations. Each
#'   iteration updates the coefficient and hyperparameter estimates based on the
#'   previous iteration. Default is `0`.
#'
#' @returns A list of class `"DLNMLPS"`. See [DLNMLPS_Object] for a description of
#' the returned object and its components.
#'
#' @export
#'
DLNM_LPS <- function(model,
                         data = NULL,
                         offset = NULL,
                         model.comparison = list(DIC = FALSE, WAIC = FALSE, CPO = FALSE),
                         family = "poisson", refit = 0){

  model.comparison <- utils::modifyList(
    list(DIC = F, WAIC = F, CPO = F),
    model.comparison
  )

  offset <- eval(substitute(offset), data, parent.frame())



  if(family == "poisson"){
    output = DLNM_LPS_pois(model,
                              data,
                              offset,
                              model.comparison,
                              refit)
  }else if (family =="NB"){
    output = DLNM_LPS_NB(model,
                             data,
                             offset,
                             model.comparison,
                             refit)

  }
}
