#' Construct a joint crossbasis matrix for two interacting exposures
#'
#' Constructs a distributed lag non-linear model (DLNM) crossbasis for the
#' interaction between two exposure variables. The exposures can be supplied
#' either as concurrent time series or as pre-computed lagged exposure histories.
#' @param x A matrix or a list containing the exposure variables.
#'   If a matrix, the first and second columns contain the time series of the
#'   first and second exposure, respectively. If a list, `x[[1]]` and `x[[2]]`
#'   should contain the lagged exposure histories for the first and second
#'   exposure.
#' @param lag Integer specifying the maximum lag.
#' @param arglag A named list of arguments defining the basis for the lag
#'   dimension (e.g., `fun`, `knots`, `nk`, `df`, `degree`,
#'   `Boundary.knots`).
#' @param argvar1 A named list of arguments defining the basis for the first
#'   exposure dimension (e.g., `fun`, `knots`, `nk`, `df`, `degree`,
#'   `Boundary.knots`).
#' @param argvar2 A named list of arguments defining the basis for the second
#'   exposure dimension (e.g., `fun`, `knots`, `nk`, `df`, `degree`,
#'   `Boundary.knots`).
#' @param group An optional factor or vector identifying groups of observations
#'   belonging to independent time series.
#' @param is_identifiable Logical; if `TRUE` (default), constructs an
#'   identifiable basis by removing confounded terms.
#'
#' @return A joint crossbasis matrix representing the interaction between two
#'   exposure variables across the specified lag period.
#' @export
#'
joint.crossbasis <- function(x = NULL, lag=10, arglag = list(), argvar1 = list(), argvar2 = list(),
                            group = NULL, is_identifiable = TRUE){

  if (is.list(x)){
    x_fit = x[[1]]
    y_fit = x[[2]]
  }else{
    x_fit = x[,1]
    y_fit = x[,2]
  }

  cn <- NULL


  basisvar1 <- do.call("joint.onebasis",utils::modifyList(argvar1, list(x = as.numeric(x_fit), type = "var", intercept = T)))
  basisvar2 <- do.call("joint.onebasis",utils::modifyList(argvar2, list(x = as.numeric(y_fit), type = "var", intercept = T)))
  basislag <- do.call("joint.onebasis",utils::modifyList(arglag, list(lag = lag, type = "lag", intercept = T)))


  nlag_basis <- ncol(basislag)
  make_lagged_basis <- function(basis_col) {
    if (!is.list(x)) {
      as.matrix(tsModel::Lag(basis_col, 0:lag, group = group))
    } else {
      matrix(basis_col, ncol = length(0:lag))
    }
  }

  basisvar1_lagged <- lapply(seq_len(ncol(basisvar1)), function(v) make_lagged_basis(basisvar1[, v]))
  basisvar2_lagged <- lapply(seq_len(ncol(basisvar2)), function(v) make_lagged_basis(basisvar2[, v]))



  if (is_identifiable){
    if (!is.list(x)){
      crossbasis <- matrix(NA, nrow = length(x_fit), ncol =  (ncol(basisvar1)-1)*ncol(basislag)+(ncol(basisvar2)-1)*ncol(basislag)+(ncol(basisvar1)-1)*(ncol(basisvar2)-1)*ncol(basislag))
    }else{
      crossbasis <- matrix(NA, nrow = dim(x_fit)[1], ncol =  (ncol(basisvar1)-1)*ncol(basislag)+(ncol(basisvar2)-1)*ncol(basislag)+(ncol(basisvar1)-1)*(ncol(basisvar2)-1)*ncol(basislag))
    }
    col_start <- 1L
    for (v in seq_len(ncol(basisvar1))){
      mat1 <- basisvar1_lagged[[v]]
      for (w in seq_len(ncol(basisvar2))){
        if (!(v==w & v==1)){
          cols <- col_start:(col_start + nlag_basis - 1L)
          crossbasis[, cols] <- (mat1 * basisvar2_lagged[[w]]) %*% basislag
          cn[cols] <- paste0("v1_", v, ".v2_", w, ".l_", seq_len(nlag_basis))
          col_start <- col_start + nlag_basis
        }
      }
    }
  } else {
    if (!is.list(x)){
      crossbasis <- matrix(0, nrow = length(x_fit), ncol = ncol(basisvar1)*ncol(basisvar2)*ncol(basislag))
    }else{
      crossbasis <- matrix(0, nrow = dim(x_fit)[1], ncol = ncol(basisvar1)*ncol(basisvar2)*ncol(basislag))
    }

    col_start <- 1L
    for (v in seq_len(ncol(basisvar1))){
      mat1 <- basisvar1_lagged[[v]]
      for (w in seq_len(ncol(basisvar2))){
        cols <- col_start:(col_start + nlag_basis - 1L)
        crossbasis[, cols] <- (mat1 * basisvar2_lagged[[w]]) %*% basislag
        cn[cols] <- paste0("v1_", v, ".v2_", w, ".l_", seq_len(nlag_basis))
        col_start <- col_start + nlag_basis
      }
    }
  }
  ind <- match(names(formals(f_formals(attributes(basisvar1)$fun))), names(attributes(basisvar1)),
               nomatch = 0)

  argvar1 <- c(attributes(basisvar1)["fun"], attributes(basisvar1)[ind])


  ind <- match(names(formals(f_formals(attributes(basisvar2)$fun))), names(attributes(basisvar2)),
               nomatch = 0)

  argvar2 <- c(attributes(basisvar2)["fun"], attributes(basisvar2)[ind])


  ind <- match(names(formals(f_formals(attributes(basislag)$fun))), names(attributes(basislag)),
               nomatch = 0)

  arglag <- c(attributes(basislag)["fun"], attributes(basislag)[ind])

  attributes(crossbasis)<-  c(attributes(crossbasis), list(argvar1 = argvar1,
                                                           argvar2 = argvar2,
                                                           arglag = arglag,
                                                           lag = c(0,lag),
                                                           range_x = c(min(x_fit, na.rm = T),max(x_fit, na.rm = T)),
                                                           range_y = c(min(y_fit, na.rm = T),max(y_fit, na.rm = T)),
                                                           is_identifiable = is_identifiable ))
  class(crossbasis) <- c("joint.crossbasis","matrix")


  colnames(crossbasis) <- cn

  return(crossbasis)

}





joint.onebasis <- function(x = NULL, type = "lag", fun = "ns", lag, ...){

  # Check input specifications
  if (type == "lag"){
    checklagbasis(fun, lag)
    if (is.null(x)) x <- 0:lag
    nl <- paste("lag",x)
    names(x) <- nl
  }
  else if (type == "var"){checkvarbasis(fun, x)}

  basis <- dlnm::onebasis(x, fun, ...)

  attributes(basis)<-  c(attributes(basis), list(fun=fun))

  class(basis) <- c("joint.onebasis","matrix")
  return(basis)
}



# Function to check input of lagbasis
checklagbasis <- function(fun, lag){
  if (!(fun %in% c("strata","ns","bs", "ps"))){
    stop("Invalid fun value.")
  }
  if (!(lag-floor(lag)==0)){
    stop("lag should be an integer.")
  }
  else return(TRUE)
}


# Function to check input of variable basis
checkvarbasis <- function(fun, x){
  if (!(fun %in% c("strata","ns","bs", "ps"))){
    stop("Invalid fun value.")
  }
  if (is.null(x)){
    stop("A vector x should be provided.")
  }
  else return(TRUE)
}

# Check formals
f_formals <- function(fun) {tryCatch(
  get(fun, mode = "function"),
  error = function(e) utils::getFromNamespace(fun, "dlnm")
)
}


