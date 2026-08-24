#' DLNMLPS objects
#'
#' Objects of class `"DLNMLPS"` as returned by [DLNM_LPS()].
#'
#' @section Components:
#' \describe{
#'   \item{model}{The fitted model.}
#'   \item{xi_mode}{Posterior mode of the latent field.}
#'   \item{Prec}{Posterior precision matrix.}
#'   \item{Xpred}{Model design matrix.}
#'   \item{y}{Response vector.}
#'   \item{offset}{Model offset.}
#'   \item{smooth}{Smooth term information.}
#'   \item{v_mode}{Posterior mode of log-hyperparameters.}
#'   \item{DIC}{Deviance Information Criterion.}
#'   \item{`pd DIC`}{Effective number of parameters for DIC.}
#'   \item{WAIC}{Widely Applicable Information Criterion.}
#'   \item{`pd WAIC`}{Effective number of parameters for WAIC.}
#'   \item{CPO}{Conditional Predictive Ordinates.}
#'   \item{constraints}{Constraint matrix to prevent identifiability issues.}
#'   \item{types}{A vector with the type of each model component.}
#' }
#'
#' @name DLNMLPS_Object
NULL
