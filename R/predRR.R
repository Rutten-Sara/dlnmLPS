#' Generate Predictions from a DLNMLPS Model
#'
#' This function generates relative risk (RR) predictions from a fitted
#' [DLNM_LPS] model over a specified grid of exposure and lag values. Predictions
#' are computed relative to a user-defined reference (centering) value.
#'
#' @param basis An object of class `"crossbasis"` or `"joint.crossbasis"`
#'   defining the exposure-lag basis used in the model.
#' @param model An object of class `"DLNMLPS"` returned by the [DLNM_LPS] function.
#' @param at_x A numeric vector specifying constant exposure values over the
#'   selected lag period, or a matrix specifying exposure histories over the
#'   same lag period used during model estimation.
#' @param at_y A numeric vector or matrix analogous to `at_x`, specifying values
#'   or histories for the second exposure when `basis` is of class
#'   `"joint.crossbasis"`.
#' @param cen A value  (or a vector of length 2 for `"joint.crossbasis"`)
#'   specifying the centering (reference) exposure level used
#'   for calculating relative risk predictions.
#' @param lag Either an integer scalar specifying the maximum lag, or a vector
#'   of length two specifying the lag range over which predictions are generated.
#' @param by.lag Numeric value specifying the increment between consecutive lag
#'   values used for prediction.
#' @param CI.all Logical value indicating whether credible intervals should be
#'   calculated only for the overall cumulative relative risk. Defaults to
#'   `FALSE`.
#' @param CI Logical value indicating whether credible intervals should be
#'   calculated for both individual lag-specific predictions and overall
#'   cumulative relative risks. Defaults to `FALSE`.
#' @param ci.level Numeric value specifying the credible interval level.
#' @param exc.prob Logical value indicating whether exceedance probabilities
#'   should be calculated. Defaults to `FALSE`.
#' @param threshold Numeric value defining the threshold used for exceedance
#'   probability calculations. Defaults to `1`.
#' @param mode Optional user-specified posterior mode.
#' @param by Level or value of the `by` variable at which the DLNM relationship
#'   is evaluated. For a continuous `by` variable, this is a (vector of) numeric
#'   value(s) defining the level at which the effect is evaluated. For a factor
#'   `by` variable, this is a (vector of) factor level(s) defining the
#'   corresponding stratum.
#'
#' @returns An object of class `"DLNMLPS_pred"` containing the following
#' components (structured by value of the `by`-variable, if specified):
#'
#' \describe{
#'   \item{logpredX}{Matrix of estimated log relative risks (logRR) for the
#'   specified combinations of exposure and lag values.}
#'   \item{Qlower_logpredX}{Matrix of lower credible interval bounds for the
#'   logRR values in `logpredX`.}
#'   \item{Qupper_logpredX}{Matrix of upper credible interval bounds for the
#'   logRR values in `logpredX`.}
#'   \item{pred_all}{Vector of estimated overall cumulative relative risks.}
#'   \item{Qlower_all}{Vector of lower credible interval bounds for the overall
#'   cumulative relative risks in `pred_all`.}
#'   \item{Qupper_all}{Vector of upper credible interval bounds for the overall
#'   cumulative relative risks in `pred_all`.}
#'   \item{exc_prob}{Vector of exceedance probabilities, if requested.}
#' }
#'
#' @export

predRR <- function(basis, model, at_x, at_y = NULL, cen, lag, by.lag = 1, CI.all = FALSE, CI = FALSE,
                   ci.level = 0.95, exc.prob = FALSE, threshold = 1, mode = NULL, by = NULL){


  constraint_matrix <- model$constraints

  index <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                 gsub("\\s+", "", colnames(model$Xpred)),
                 fixed = TRUE)


  rhs <-  model$model[[3]]
  all_terms <- split_terms(rhs)

  index_term <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                 gsub("\\s+", "", all_terms),
                 fixed = TRUE)

  kind_term <- model$types[which(index_term)-1]

  dim_basis <- ncol(basis)

  if(sum(index)==dim_basis & (!is.null(by))){
    stop("No by variable for fixed crossbasis")
  }


  if(is.null(mode)){
    mode <- model$xi_mode
  }

  mode_cb <- mode[index]
  names(mode_cb) <- colnames(model$Xpred)[index]




  ind_cb <- which(index)[1:dim_basis]
  xi_mode_cb = mode[ind_cb]


  if(exc.prob){
    CI.all = TRUE
  }


  origlag <- attr(basis,"lag")
  if (is.null(lag)) {lag = origlag}


  if(length(lag) == 1){
    lag <- seq(0,lag, by = by.lag)
  }else if(length(lag)==2){
    lag <- seq(lag[1],lag[2], by = by.lag)
  }else{
    stop("lag should be an integer scalar or vector of length 2")
  }

  prec_xi <- model$Prec


  if(attr(basis, "class")[1] == "joint.crossbasis"){



    if (is.matrix(at_x)){
      if (dim(at_x)[2] !=  sum(attr(basis,"lag"))+1)
        stop("matrix in 'at' must have ncol=length(origlag)")
      if (is.null(rownames(at_x))){
        rownames(at_x) <- seq(nrow(at_x))}
      if (is.null(rownames(at_y))){
        rownames(at_y) <- seq(nrow(at_y))}
    }

    predvar <- if(is.matrix(at_x)) {c(rownames(at_x),"cen")
    } else c(at_x,cen[1])


    # Data matrix for prediction
    Xpred <- data_pred_2d(basis, at_x, at_y, lag, predvar, cen)

    predvar <- predvar[-length(predvar)]

  }else{

    if (is.matrix(at_x)){
      if (dim(at_x)[2] !=  sum(attr(basis,"lag"))+1)
        stop("matrix in 'at' must have ncol=length(origlag)")
      if (is.null(rownames(at_x))){
        rownames(at_x) <- seq(nrow(at_x))}
    }

    predvar <- if(is.matrix(at_x)) {rownames(at_x)
    } else at_x


    # Data matrix for prediction
    Xpred <- data_pred_1d(basis, at_x, lag, cen)

  }

  if(is.null(by)){

  logpredX <- as.numeric(Xpred %*% xi_mode_cb)

  if(CI){
    if(is.null(constraint_matrix)){
      sd_logpredX <- Predict_se(prec_xi,Xpred,ind_cb)
    }else{
      sd_logpredX <- Predict_se_con(prec_xi,Xpred,ind_cb, constraint_matrix)
    }

    quantiles <- mapply(function(mean_val, sd_val) {
      stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
    }, logpredX, sd_logpredX)
    Qlower_logpredX <- quantiles[1,]
    Qupper_logpredX <- quantiles[2,]

    Qlower_logpredX <- matrix(Qlower_logpredX, nrow = length(predvar), byrow = FALSE)
    Qupper_logpredX <- matrix(Qupper_logpredX, nrow = length(predvar), byrow = FALSE)

  }else{
    sd_logpredX <- Qlower_logpredX <- Qupper_logpredX <- NULL

  }


  # Overall RR
  Xpredall <- 0
  for (j in 1:length(lag)) {
    if(lag[j]==round(lag[j])){
      ind_all <- seq(length(predvar)) + length(predvar) * (j - 1)
      Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }
  }


  pred_all <- exp(as.numeric(Xpredall %*% xi_mode_cb))


  if(CI | CI.all){
    if(is.null(constraint_matrix)){
      sd_all <- Predict_se(prec_xi,Xpredall,ind_cb)
    }else{
      sd_all <- Predict_se_con(prec_xi,Xpredall,ind_cb, constraint_matrix)
    }

    quantiles_all <- mapply(function(mean_val, sd_val) {
      stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
    }, log(pred_all), sd_all)
    Qlower_all = exp(quantiles_all[1,])
    Qupper_all = exp(quantiles_all[2,])
  }else{
    sd_all <- Qlower_all <- Qupper_all <- NULL
  }


  if(exc.prob){
    prob_all <- mapply(function(mean_val, sd_val) {
      stats::pnorm(log(threshold), mean = mean_val, sd = sd_val, lower.tail = FALSE)
    }, log(pred_all), sd_all)
    exc.prob_res = prob_all
  }else{
    exc.prob_res = NULL
  }



  output <- list("logpredX" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                 "Qlower_logpredX" =  Qlower_logpredX,
                 "Qupper_logpredX" =  Qupper_logpredX,
                 "Qlower_all" = Qlower_all,
                 "Qupper_all" = Qupper_all,
                 "pred_all" = pred_all,
                 "exc.prob" = exc.prob_res,
                 "x" = at_x,
                 "y" = at_y,
                 "lag" = lag)
  } else if (kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
    if(is.null(by)){stop("Specify at least one region ID")}

    output = list()

    Xpred_full = cbind(Xpred, Xpred)

    # Overall RR
    Xpredall <- 0
    for (j in 1:length(lag)) {
      if(lag[j]==round(lag[j])){
        ind_all <- seq(length(predvar)) + length(predvar) * (j - 1)
        Xpredall <- Xpredall + Xpred_full[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
      }
    }

    ind.cb_spat_list = list()

    for(i in factor(by)){

      ind_spat_ID <- grepl(i, colnames(model$Xpred),
                           fixed = TRUE)

      ind_spat <- index & ind_spat_ID

      xi_mode_cb_spat_i = mode[ind_spat]
      ind.cb_spat_list[[i]] = which(ind_spat)

      ##########################
      logpredX <- as.numeric(Xpred_full %*% c(xi_mode_cb, xi_mode_cb_spat_i))

      pred_all <- exp(as.numeric(Xpredall %*% c(xi_mode_cb, xi_mode_cb_spat_i)))

      output[[i]] <- list("logpredX" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                          "pred_all" = pred_all,
                          "by" = i,
                          "x" = at_x,
                          "y" = at_y,
                          "lag" = lag)


    }


    if(CI | CI.all){
      if(is.null(constraint_matrix)){
        sd_all <- Predict_se_all(prec_xi,Xpredall,ind_cb, ind.cb_spat_list)
      }else{
        sd_all <- Predict_se_all_con(prec_xi,Xpredall,ind_cb, ind.cb_spat_list, constraint_matrix)
      }


      for(i in 1:length(by)){
        quantiles_all <- mapply(function(mean_val, sd_val) {
          stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
        }, log(output[[i]]$pred_all), sd_all[i,])
        output[[i]]$Qlower_all = exp(quantiles_all[1,])
        output[[i]]$Qupper_all = exp(quantiles_all[2,])

        if(exc.prob){
          prob_all <- mapply(function(mean_val, sd_val) {
            stats::pnorm(log(threshold), mean = mean_val, sd = sd_val, lower.tail = FALSE)
          }, log(output[[i]]$pred_all), sd_all[i,])
          output[[i]]$exc.prob = prob_all
        }
      }
    }


    if(CI){
      if(is.null(constraint_matrix)){
        sd_logpredX <- Predict_se_all(prec_xi,Xpred_full,ind_cb, ind.cb_spat_list)
      }else{
        sd_logpredX <- Predict_se_all_con(prec_xi,Xpred_full,ind_cb, ind.cb_spat_list, constraint_matrix)
      }


      for(i in 1:length(by)){
        quantiles <- mapply(function(mean_val, sd_val) {
          stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
        }, output[[i]]$logpredX, sd_logpredX[i,])
        output[[i]]$Qlower_logpredX <- matrix(quantiles[1,], nrow = length(predvar), byrow = FALSE)
        output[[i]]$Qupper_logpredX <- matrix(quantiles[2,], nrow = length(predvar), byrow = FALSE)
      }

    }




  } else if (kind_term %in% c("modified DLNM", "joint modified DLNM")){
    type_smooth <- model$smooth[[which(index_term)-1]]
    output = list()

    for(by.value in 1:length(by)){

    if("ns" %in% class(type_smooth)){
      arg.by = as.matrix(splines::ns(by[by.value], knots = attributes(type_smooth)$knots,
                                     Boundary.knots = attributes(type_smooth)$Boundary.knots,
                                     intercept = attributes(type_smooth)$intercept))
    }else if("bs" %in% class(type_smooth)){
      arg.by = as.matrix(splines::bs(by[by.value], knots = attributes(type_smooth)$knots,
                                     Boundary.knots = attributes(type_smooth)$Boundary.knots,
                                     intercept = attributes(type_smooth)$intercept,
                                     degree = attributes(type_smooth)$degree))
    } else if ("ps" %in% class(type_smooth)){
      arg.by = as.matrix(dlnm::ps(by[by.value], knots = attributes(type_smooth)$knots,
                                  df = attributes(type_smooth)$df,
                                  degree = attributes(type_smooth)$degree,
                                  intercept = attributes(type_smooth)$intercept))
    } else if (!is.null(attributes(type_smooth)$Dimnames[[1]])){
      index_by <- grepl(by[by.value], attributes(type_smooth)$Dimnames[[2]],
                        fixed = TRUE)
      arg.by <-rep(0,length(attributes(type_smooth)$Dimnames[[2]]))
      arg.by[index_by] =1
      arg.by <- as.matrix(arg.by, nrow = 1)
    } else{
      arg.by = as.matrix(by[by.value], ncol = 1)
    }

    arg.by <- matrix(rep(arg.by, each = length(predvar)*length(lag)), nrow = length(predvar)*length(lag), byrow = FALSE) # Replicate values for every lag


    Mat_crossbasis_spat <- Xpred


    interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
      arg.by * Mat_crossbasis_spat[, j]
    })

    interaction_matrix <- do.call(cbind, interaction_list)

    # Global design matrix
    Xpred_full <- cbind(Xpred,interaction_matrix)

    # Overall RR
    Xpredall <- 0
    for (j in 1:length(lag)) {
      if(lag[j]==round(lag[j])){
        ind_all <- seq(length(predvar)) + length(predvar) * (j - 1)
        Xpredall <- Xpredall + Xpred_full[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
      }
    }


    ##########################
    logpredX <- as.numeric(Xpred_full %*% mode_cb)

    if( CI ){
      if(is.null(constraint_matrix)){
        sd_logpredX <- Predict_se(prec_xi,Xpred_full,which(index))
      }else{
        sd_logpredX <- Predict_se_con(prec_xi,Xpred_full,which(index), constraint_matrix)
      }

      quantiles <- mapply(function(mean_val, sd_val) {
        stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
      }, logpredX, sd_logpredX)

      Qlower_logpredX <- quantiles[1,]
      Qupper_logpredX <- quantiles[2,]

      Qlower_logpredX <- matrix(Qlower_logpredX, nrow = length(predvar), byrow = FALSE)
      Qupper_logpredX <- matrix(Qupper_logpredX, nrow = length(predvar), byrow = FALSE)

    }else{
      Qlower_logpredX = NULL
      Qupper_logpredX = NULL
    }


    pred_all <- exp(as.numeric(Xpredall %*% mode_cb))

    if(CI | CI.all){
      if(is.null(constraint_matrix)){
        sd_all <- Predict_se(prec_xi,Xpredall,which(index))

      }else{
        sd_all <- Predict_se_con(prec_xi,Xpredall,which(index), constraint_matrix)
      }
      quantiles_all <- mapply(function(mean_val, sd_val) {
        stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
      }, log(pred_all), sd_all)
      Qlower_all = exp(quantiles_all[1,])
      Qupper_all = exp(quantiles_all[2,])
    }else{
      Qlower_all = NULL
      Qupper_all = NULL
    }

    if(exc.prob){
      prob_all <- mapply(function(mean_val, sd_val) {
        stats::pnorm(log(threshold), mean = mean_val, sd = sd_val, lower.tail = FALSE)
      }, log(pred_all), sd_all)
      exc.prob_res = prob_all
    }else{
      exc.prob_res = NULL
    }

    output[[by.value]] <- list("logpredX" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                   "Qlower_logpredX" = Qlower_logpredX,
                   "Qupper_logpredX" = Qupper_logpredX,
                   "Qlower_all" = Qlower_all,
                   "Qupper_all" = Qupper_all,
                   "pred_all" = pred_all,
                   "exc.prob" = exc.prob_res,
                   "x" = at_x,
                   "y" = at_y,
                   "lag" = lag,
                   "by" = by[by.value])
    }

  }




  class(output) <- c("DLNMLPS_pred", kind_term)

  return(output)



}



# Construction of data matrix
data_pred_2d <- function(basis1, at_x, at_y, predlag, predvar, cen){
  ################################################################################
  if (is.matrix(at_x)){
    nrow <- dim(at_x)[1]
    if (!is.null(cen)) { #bind center to matrix
      at_x = rbind(at_x,rep(cen[1],dim(at_x)[2]))
      at_y = rbind(at_y,rep(cen[2],dim(at_x)[2]))
      nrow = nrow+1
    }
    at_x <- as.numeric(at_x)
    at_y <- as.numeric(at_y)
  } else{
    nrow <- length(at_x)
    if (!is.null(cen)) { #add center to vector
      at_x = append(at_x,cen[1])
      at_y = append(at_y,cen[2])
      nrow = nrow+1
    }
    at_x <- rep(at_x, length(predlag))
    at_y <- rep(at_y, length(predlag))
  }

  # Construct one basis in prediction values
    basisvar1 <- do.call("joint.onebasis", c(list(x = at_x, type="var"), attr(basis1,"argvar1")))
    basisvar2 <- do.call("joint.onebasis", c(list(x = at_y, type="var"), attr(basis1,"argvar2")))
    basislag <- do.call("joint.onebasis", c(list(x = predlag, lag = attributes(basis1)$lag[2], type="lag"), attr(basis1,"arglag")))

    # Construct design matrix

    if (attributes(basis1)$is_identifiable){
      XYpred <- matrix(0,nrow = length(at_x),ncol = (ncol(basisvar1)-1)*ncol(basislag)+(ncol(basisvar2)-1)*ncol(basislag)+(ncol(basisvar1)-1)*(ncol(basisvar2)-1)*ncol(basislag))
      for (l in seq(length = length(predlag))){
        for (v in seq(length = ncol(basisvar1))){
          for (w in seq(length = ncol(basisvar2))){
            if (!(v==w & v==1)){
              for (k in 1:ncol(basislag)){
                XYpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*ncol(basisvar2)*(v-1)+ncol(basislag)*(w-1)+k-ncol(basislag))] = (basisvar1[((l-1)*nrow+1):(nrow*l),v]*basisvar2[((l-1)*nrow+1):(nrow*l),w])*basislag[l,k]
              }
            }
          }
        }
      }

    } else {
      XYpred <- matrix(0,nrow = length(at_x),ncol = ncol(basisvar1)*ncol(basisvar2)*ncol(basislag))
      for (l in seq(length = length(predlag))){
        for (v in seq(length = ncol(basisvar1))){
          for (w in seq(length = ncol(basisvar2))){
            for (k in 1:ncol(basislag)){
              XYpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*ncol(basisvar2)*(v-1)+ncol(basislag)*(w-1)+k)] = (basisvar1[((l-1)*nrow+1):(nrow*l),v]*basisvar2[((l-1)*nrow+1):(nrow*l),w])*basislag[l,k]
            }
          }
        }
      }

    }




  # The predicted values are now expressed versus the center (which is the last row of each lag in XYpred)
  # Therefore, subtract this last row of each lag from every other row of that lag

  diff <- matrix(0, nrow = dim(XYpred)[1]-length(predlag), ncol=ncol(XYpred))
  for (i in seq(length(predlag))){
    ind <- seq(length(predvar)) + (length(predvar)) * (i - 1) # length predvar includes center
    ind_diff <- seq(length(predvar)-1) + (length(predvar)-1) * (i - 1)
    diff[ind_diff,] <- XYpred[ind[1:(length(ind)-1)],]-matrix(rep(XYpred[ind[length(ind)],],each=length(ind)-1),nrow=length(ind)-1)

  }

  return(diff)
}


# Construction of data matrix
data_pred_1d <- function(basis, at_x, predlag, cen){
  ################################################################################
  if (is.matrix(at_x)){
    nrow <- dim(at_x)[1]
    at_x <- as.numeric(at_x)
  } else{
    nrow <- length(at_x)
    at_x <- rep(at_x, length(predlag))
  }

  basisvar <- do.call(dlnm::onebasis, c(list(x = at_x), attr(basis,"argvar")))
  basislag <- do.call(dlnm::onebasis, c(list(x = predlag), attr(basis,"arglag")))

  # basis variables for center
  basiscen <- do.call(dlnm::onebasis, c(list(x = cen), attr(basis,"argvar")))
  # center data matrix
  Xpred_cen <- scale(basisvar,center = basiscen, scale = FALSE)

  # Prediction matrix
  Xpred <- matrix(0, nrow=length(at_x), ncol = ncol(basis))

  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(basisvar))){
      for (k in 1:ncol(basislag)){
        Xpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
      }
    }
  }

  return(Xpred)
}


