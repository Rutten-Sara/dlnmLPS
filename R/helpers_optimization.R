# Helpers optimization function   ---------------------------------------------

#' Specify a smooth or structured model component
#'
#' `sm()` specifies smooth functions and structured effects to be included
#' in the model formula supplied to [DLNM_LPS()].
#'
#' @param x The variable or basis to which the model component is applied.
#'   This can be:
#'   \describe{
#'     \item{Crossbasis object}{An object of class [`dlnm::crossbasis()`],
#'     representing a distributed lag non-linear model
#'     term.}
#'     \item{Numeric vector}{A continuous covariate for constructing a
#'     smooth function.}
#'     \item{Factor vector}{A categorical covariate for defining random
#'     effects or other structured effects.}
#'   }
#' @param arg A list specifying additional arguments. The argument
#'   `arg$type` determines the type of component:
#'   \describe{
#'     \item{`crossbasis`}{A penalized or unpenalized DLNM cross-basis.
#'     Objects of class `crossbasis` are recognized
#'     automatically and therefore do not need to be specified explicitly in
#'     `arg$type`. The additional argument `pen` can be specified, where
#'     `NULL` (default) specifies an unpenalized basis and an integer
#'     specifies the penalty structure.}
#'     \item{`ps`, `ns`, `bs`}{Smooth functions based on P-splines,
#'     natural splines, or B-splines. Additional arguments can include
#'     those passed to the corresponding spline constructors, such as
#'     `df`, `knots`, `intercept`, `Boundary.knots`, `degree`, and
#'     penalty specifications.}
#'     \item{`ind`, `ICAR`, `Leroux`, `Convolution`}{Random
#'     effects, including spatial random effects. Additional arguments
#'     include spatial connectivity information through `map` and
#'     `connect` (logical value indicating whether the spatial adjacancy
#'     graph should be connected, defaults to `TRUE`).}
#'     \item{`rw1`, `rw2`}{First- and second-order random walks.
#'     Additional arguments include `cyclic` for specifying cyclic random
#'     walks (Defaults to `FALSE`).}
#'   }
#' @param by An optional variable used to modify or stratify the model
#'   component specified by `x`. This can be:
#'   \describe{
#'     \item{Numeric vector}{A continuous variable used to define a
#'     varying or modified effect of the component. For a crossbasis
#'     specified in `x`, this can be used to define a continuously varying
#'     DLNM effect. For a (spatial) random effect specified in `x`, this can be
#'     used to define a temporally varying spatial effect, when the vector
#'     is a time indicator.}
#'     \item{Factor vector}{A categorical variable used to define
#'     stratum-specific effects. For a crossbasis/random walk specified in `x`,
#'     this allows the DLNM association/random walk to spatially vary when the
#'     factor is a region indicator.
#'     }
#'   }
#' @param arg.by A list specifying additional arguments for the `by`
#'   variable. The argument `arg.by$type` determines how the `by` variable
#'   is modelled:
#'   \describe{
#'     \item{`ps`, `ns`, `bs`}{Smooth functions based on P-splines,
#'     natural splines, or B-splines. Additional arguments can include
#'     those passed to the corresponding spline constructors, such as
#'     `df`, `knots`, `intercept`, `Boundary.knots`, `degree`, and
#'     penalty specifications.}
#'     \item{`ind`, `ICAR`, `Leroux`}{Random effects, including spatial
#'     random effects. Additional arguments include spatial connectivity
#'     information through `map` and `connect` (logical value indicating
#'     whether the spatial adjacency graph should be connected, defaults
#'     to `TRUE`).}
#'     \item{`rw1`, `rw2`}{First- and second-order random walks.
#'     Additional arguments include `cyclic` for specifying cyclic random
#'     walks (Defaults to `FALSE`). Only for time-structured (spatial)
#'     random effects (main type is `ind`, `ICAR` or `Leroux`).}
#'   }
#'
#' @returns
#' An object representing a smooth or structured model component. This object
#' is used internally by [DLNM_LPS()] when constructing the model from the
#' supplied formula.
#'
#' @examples
#' # `sm()` is intended to be used inside a DLNM_LPS() formula:
#' \dontrun{
#' # Penalized DLNM cross-basis
#' # Cross-basis objects are recognized automatically
#' sm(crossbasis, arg = list(pen = 2))
#'
#' # Leroux spatial random effect
#' sm(region_ID, arg = list(type = "Leroux", map = map))
#'
#' # Cyclic second-order random walk
#' sm(time, arg = list(type = "rw2", cyclic = TRUE))
#'
#' # Penalized P-spline smooth
#' sm(x, arg = list(type = "ps", pen = 2, df = 5))
#'
#' # Penalized DLNM cross-basis by location (region_ID)
#' sm(
#'   crossbasis,
#'   arg = list(pen = 2),
#'   by = region_ID,
#'   arg.by = list(type = "Leroux", map = map)
#' )
#'
#' # Penalized DLNM cross-basis by continuous modifier x (using ps)
#' sm(
#'   crossbasis,
#'   arg = list(pen = 2),
#'   by = x,
#'   arg.by = list(type = "ps", pen = 2, df = 5)
#' )
#'
#' # Spatially varying second-order random walk
#' sm(
#'   time,
#'   arg = list(type = "rw2", cyclic = TRUE),
#'   by = region_ID,
#'   arg.by = list(type = "Leroux", map = map)
#' )
#'
#' # Time varying spatial effect
#' sm(
#'   region_ID,
#'   arg = list(type = "Leroux", map = map),
#'   by = time,
#'   arg.by = list(type = "rw1", cyclic = TRUE)
#' )
#' }
sm <- function(x , arg = NULL, by = NULL, arg.by = NULL){

}

#' Splitting the model formula in different terms
#' @noRd
split_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    c(split_terms(expr[[2]]), split_terms(expr[[3]]))
  } else {
    list(expr)
  }
}

#' Verifying whether term is a smooth function ('sm')
#' @noRd
is_sm <- function(term) {
  is.call(term) && identical(term[[1]], as.name("sm"))
}


#' Evaluating the expression
#' @noRd
eval_expression <- function(x_expr, data = NULL, env = globalenv()){
  if(is.symbol(x_expr) | is.language(x_expr)){
    x_name <- deparse(x_expr)

    if(!is.null(data) & x_name %in% names(data)){
      data[[x_name]]
    }else{
      eval(x_expr, envir = env)
    }
  }else{
    x_expr
  }
}


#' Initializing the setup
#' @noRd
initialize_laplace_setup_constraint <- function(model, data, include_dispersion = FALSE) {
  rhs <- model[[3]]
  all_terms <- split_terms(rhs)


  has_no_intercept <- attr(stats::terms(model), "intercept") == 0

  is_intercept_term <- vapply(all_terms, function(x) {
    identical(x, 0) ||
      identical(x, 1) ||
      (is.call(x) &&
         identical(x[[1]], as.name("-")) &&
         identical(x[[2]], 1))
  }, logical(1))

  all_terms <- all_terms[!is_intercept_term]

  is_sm_res <- vapply(all_terms, is_sm, logical(1))
  all_fixed_terms <- all_terms[!is_sm_res]
  all_sm_terms    <- all_terms[is_sm_res]

  all_terms <- c(all_sm_terms, all_fixed_terms )


  if(!is.null(data)){
    mod <- stats::model.frame( paste(deparse1(model[[2]], backtick = TRUE),model[[1]],1),
                               data = data, na.action = stats::na.pass)
  }else{
    mod <- stats::model.frame( paste(deparse1(model[[2]], backtick = TRUE),model[[1]],1),
                               na.action = stats::na.pass)
  }

  n_terms <- length(all_terms)
  y <- as.numeric(stats::model.extract(mod, "response"))

  v_init <- if (include_dispersion) 0 else NULL
  names(v_init) <- if (include_dispersion) "scale" else NULL
  output <- list(
    all_terms = all_terms,
    y = y, # response,
    Wcb = NULL,
    pen_dlnm = if (include_dispersion) 1 else 0,
    Pv.list = list(),
    logpv.fixed = list(),
    v_init = v_init,
    sm_list = vector("list", n_terms),
    constraint = NULL,
    has_no_intercept = has_no_intercept,
    has_no_intercept_smooth = has_no_intercept
  )

  output


}

#' Checking validitiy of the input
#' @noRd
check_sm_term_inputs <- function(sm_args_iter, data, env) {

  if (is.null(sm_args_iter$x)) {
    stop("No x vector/matrix specified in function sm")
  }
  sm_args_iter$x <- eval_expression(sm_args_iter$x, data, env)

  sm_args_iter$arg <- lapply(sm_args_iter$arg, eval, envir = env)

  if (("joint.crossbasis" %in% class(sm_args_iter$x)) ||
      ("crossbasis" %in% class(sm_args_iter$x))) {
    sm_args_iter$arg$type <- "crossbasis"
  }

  if (!"type" %in% names(sm_args_iter$arg) || is.null(sm_args_iter$arg$type)) {
    stop("Please specify the type of your sm function in arg.")
  }

  if (!sm_args_iter$arg$type %in% c("crossbasis", "ps", "ns", "bs", "ind", "ICAR","Leroux", "Convolution", "rw1", "rw2")) {
    stop("sm type must be either crossbasis, ps, ns, bs, rw1, rw2, ind, ICAR, Leroux or Convolution")
  }

  if (sm_args_iter$arg$type %in% c("ind","ICAR","Leroux","Convolution","rw1","rw2")) {
    if(sm_args_iter$arg$type %in% c("rw1","rw2")){
      lev <- sort(unique(sm_args_iter$x))
      sm_args_iter$x <- factor(sm_args_iter$x, levels = lev)
    }else{
      sm_args_iter$x <- as.factor(sm_args_iter$x)
    }

    sm_args_iter$arg <- utils::modifyList(
      list( map = NULL, connect = TRUE, cyclic = FALSE),
      sm_args_iter$arg
    )
  } else if (sm_args_iter$arg$type == "crossbasis") {
    sm_args_iter$arg <- utils::modifyList(
      list(type = "crossbasis", pen = NULL, map = NULL),
      sm_args_iter$arg
    )
  } else {
    sm_args_iter$arg <- utils::modifyList(
      list(type = "ps", pen = 2, df = 10,
           knots = NULL, intercept = FALSE,
           Boundary.knots = range(sm_args_iter$x, na.rm = TRUE),
           degree = 3),
      sm_args_iter$arg
    )
  }


  if(sm_args_iter$arg$type == "crossbasis" & !is.null(sm_args_iter$by) & is.null(sm_args_iter$arg.by)){
    stop("Specify the arguments of your by variable in the sm function")
  }

  if(!is.null((sm_args_iter$arg.by))){

    sm_args_iter$arg.by <- lapply(sm_args_iter$arg.by, eval, envir = env)

    if (! "type" %in% names(sm_args_iter$arg.by)){
      stop("Specify a type argument in arg.by that is equal to ICAR, Leroux, ind, rw1, rw2, ps, ns, bs, linear or factor")
    }
    if(!sm_args_iter$arg.by$type %in% c("ind", "Leroux","ICAR","rw1","rw2", "ps","ns","bs", "linear","factor")) {
      stop("Specify a type argument in arg.by that is equal to ICAR, Leroux, ind, rw1, rw2, ps, ns, bs, linear or factor")
    }
    sm_args_iter$by  = eval_expression(sm_args_iter$by, data, env)

    if(sm_args_iter$arg.by$type%in% c("ind","ICAR","Leroux","rw1","rw2")){
      if(sm_args_iter$arg.by$type %in% c("rw1","rw2")){
        lev <- sort(unique(sm_args_iter$by))
        sm_args_iter$by <- factor(sm_args_iter$by, levels = lev)
      }else{
        sm_args_iter$by = as.factor(sm_args_iter$by )
      }

      sm_args_iter$arg.by <- utils::modifyList(
        list(map = NULL, connect = TRUE, cyclic = FALSE),
        sm_args_iter$arg.by
      )
    }else if (sm_args_iter$arg.by$type == "linear"){
      sm_args_iter$arg.by <- list(type = "linear")

    }else if (sm_args_iter$arg.by$type == "factor"){
      sm_args_iter$arg.by <- list(type = "factor")
    }else{
      sm_args_iter$arg.by <- utils::modifyList(
        list(type = "ps", pen = 2, df = 10,
             knots = NULL, intercept = TRUE, Boundary.knots = range(sm_args_iter$by, na.rm = TRUE),
             degree = 3),
        sm_args_iter$arg.by
      )
    }
  }

  sm_args_iter
}

