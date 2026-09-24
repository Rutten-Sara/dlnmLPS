#' Summarize a DLNMLPS model
#'
#' Computes a summary of the fixed effects of a fitted `"DLNMLPS"` model.
#' The summary includes posterior means (modes), standard deviations, and
#' approximate 95% credible intervals based on the Laplace approximation.
#' @param object A fitted `"DLNMLPS"` model.
#' @param ... Additional arguments to be passed.
#' @returns An object of class \code{"summary.DLNMLPS"}.
#' @export
#'
summary.DLNMLPS <- function(object,...){

  constraint_matrix <- object$constraints

  mode <- object$xi_mode
  Prec <- object$Prec

  Xpred <- Matrix::diag(1, nrow = length(mode))

  if(is.null(constraint_matrix)){
    sd_logpredX <- Predict_se(Prec,Xpred,seq_along(mode))
  }else{
    sd_logpredX <- Predict_se_con(Prec,Xpred,seq_along(mode), constraint_matrix)
  }

  quantiles <- mapply(function(mean_val, sd_val) {
    stats::qnorm(p = c(0.025,0.5,0.975), mean = mean_val, sd = sd_val)
  },  mode, sd_logpredX)


  summary_fixed <- data.frame(
    mean       = mode,
    sd         = sd_logpredX,
    `0.025quant` = quantiles[1,],
    `0.5quant`   =mode,
    `0.975quant` = quantiles[3,],
    mode       = mode,
    check.names  = FALSE
  )

  rownames(summary_fixed) <- names(mode)

  if(!is.null(object$v_mode)){
    hyperparameters <- make_summary_hyper(object$v_mode)
  }else{
    hyperparameters <- NULL
  }

  structure(
    list(
      call = object$model,
      summary.coef = summary_fixed,
      hyperparameters = hyperparameters,
      dic = object$DIC,
      pd.dic = object$`pd DIC`,
      waic = object$WAIC,
      pd.waic = object$`pd WAIC`,
      cpo = object$CPO
    ),
    class = "summary.DLNMLPS"
  )

}

#' Print a summary of a DLNMLPS model
#'
#' @param x An object of class \code{"summary.DLNMLPS"}.
#' @param ... Additional arguments to be passed.
#' @export
print.summary.DLNMLPS <- function(x, ...) {

  cat("\nCall:\n\n")
  cat(deparse1(x$call))

  cat("\n\nFixed effects:\n\n")
  print(round(x$summary.coef, 3))

  if (!is.null(x$hyperparameters)) {

    cat("\nHyperparameters:\n\n")

    hp <- x$hyperparameters

    for (i in seq_len(nrow(hp))) {

      if (hp$hyperparameter[i] == "scale") {
        cat("scale  ", hp$value[i], "\n", sep = "")
      } else {
          cat(
            hp$hyperparameter[i], " ",
            hp$type[i], " ",
            hp$variable[i], "  ",
            hp$value[i],
            "\n",
            sep = ""
          )
      }
    }
  }

  if (!is.null(x$dic) ||
      !is.null(x$waic) ||
      !is.null(x$cpo)) {

    cat("\nModel comparison:\n\n")

    if (!is.null(x$dic)) {
      cat(
        "DIC:", x$dic,
        "; Effective Parameters:", x$pd.dic, "\n"
      )
    }

    if (!is.null(x$waic)) {
      cat(
        "WAIC:", x$waic,
        "; Effective Parameters:", x$pd.waic, "\n"
      )
    }

    if (!is.null(x$cpo)) {
      cat(
        "sum(log(CPO)):",
        sum(log(x$cpo), na.rm = TRUE),
        "\n"
      )
    }
  }

  invisible(x)
}


#' Construction of hyperparameter summary
#' @noRd
make_summary_hyper <- function(v_mode) {

  nm <- names(v_mode)

  # Identify type
  type_hyper <- ifelse(grepl("^Rho", nm),"rho","exp")

  # Extract the hyperparameter type
  hyperparameter <-  ifelse(
    nm == "scale",
    "scale",
    sub("^((Penalty|Precision|Rho)[^ ]*).*", "\\1", nm)
  )

  # Extract the model type
  type <-  ifelse(
    nm == "scale",
    NA_character_,
    sub(
      "^(Penalty|Precision|Rho)\\s+(.*?)\\s+sm\\(.*$",
      "\\2",
      nm
    )
  )

  # Extract the variable from sm(...)
  variable <- ifelse(
    nm == "scale",
    NA_character_,
    sub(".*sm\\(([^,\\)]+).*", "\\1", nm)
  )

  # Clean whitespace
  variable <- trimws(variable)

  # Transform values
  value <- v_mode

  value[type_hyper == "exp"] <- exp(value[type_hyper == "exp"])

  value[type_hyper == "rho"] <- leroux_rho(value[type_hyper == "rho"])

  key <- paste(hyperparameter, type, variable, sep = "|||")
  groups <- split(seq_along(value), key)



  out <- lapply(groups, function(i) {
    data.frame(
      hyperparameter = hyperparameter[i[1]],
      type = type[i[1]],
      variable = variable[i[1]],
      value = paste(
        format(value[i], digits = 4, trim = TRUE),
        collapse = " "
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}
