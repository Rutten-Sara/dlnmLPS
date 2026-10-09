#' Generate Predictions from a DLNMLPS Model
#'
#' This function generates predictions from a fitted
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
#'   probability calculations (at scale of linear predictor). Defaults to `0`.
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
#'   \item{matfit}{Matrix of predictions for the
#'   specified combinations of exposure and lag values.}
#'   \item{Qlower_matfit}{Matrix of lower credible interval bounds for the
#'   values in `matfit`.}
#'   \item{Qupper_matfit}{Matrix of upper credible interval bounds for the
#'   values in `matfit`.}
#'   \item{matRRfit}{Matrix of exponentiated associations from `matfit`.}
#'   \item{Qlower_matRRfit}{Matrix of exponentiated bounds from `Qlower_matfit`.}
#'   \item{Qupper_matRRfit}{Matrix of exponentiated bounds from `Qupper_matfit`.}
#'   \item{allfit}{Vector of estimated overall cumulative predictions.}
#'   \item{Qlower_allfit}{Vector of lower credible interval bounds for the overall
#'   cumulative predictions in `allfit`.}
#'   \item{Qupper_allfit}{Vector of upper credible interval bounds for the overall
#'   cumulative predictions in `allfit`.}
#'   \item{allRRfit}{Matrix of exponentiated associations from `allfit`.}
#'   \item{Qlower_allRRfit}{Matrix of exponentiated bounds from `Qlower_allfit`.}
#'   \item{Qupper_allRRfit}{Matrix of exponentiated bounds from `Qupper_allfit`.}
#'   \item{exc_prob_matfit}{Matrix of exceedance probabilities for matfit, if requested.}
#'   \item{exc_prob_allfit}{Vector of exceedance probabilities for allfit, if requested.}
#' }
#'
#' @examples
#' data(barcelona_temp)
#' data(shapefile_bcn)
#'
#' # Examples of applications
#'
#' dataset <- barcelona_temp |>
#'   dplyr::mutate(
#'     dow = weekdays(date),
#'     year = lubridate::year(date),
#'     month = lubridate::month(date)
#'   )
#'
#' group <- factor(paste0(dataset$District, "-", dataset$year))
#'
#'
#' ### Example 1: Penalized and common DLNM
#'
#' # Create the (penalized) crossbasis
#' basis.max <- dlnm::crossbasis(
#'     dataset$Tmax,
#'     argvar = list(df = 5, fun = "ps"),
#'     arglag = list(df = 5, fun = "ps"),
#'     lag = 14,
#'     group = group
#' )
#'
#' # Fit model
#' model_max_common <- DLNM_LPS(
#'     dtot ~ 1 + factor(dow) +
#'     sm(month, arg = list(type = "ns", df = 3)) +
#'     sm(year, arg = list(type = "rw1")) +
#'     sm(basis.max, arg = list(pen= 2)) +
#'     sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
#'     data = dataset,
#'     offset = log(pop),
#'     model.comparison = list(DIC = TRUE),
#'     family = "Poisson"
#' )
#'
#' # Make predictions
#' pred_max <- pred_crossbasis(
#'      basis = basis.max,
#'      model = model_max_common,
#'      at_x = 22:30,
#'      cen = 25,
#'      lag = c(0,14),
#'      CI = TRUE
#' )
#'
#' # Plot results (lag-specific and overall RR)
#' plot(pred_max, slice_lag = c(0,3,8), ptype = "slices", CI = TRUE, xlab = "temp", ylab = "RR")
#' plot(pred_max, ptype = "overall", CI = TRUE, xlab = "temp", ylab = "RR")
#'
#'
#' ### Example 2: Penalized and common joint DLNM
#'
#' # Create the (penalized) joint crossbasis
#' basis.max_min <- joint.crossbasis(
#'     cbind(dataset$Tmax, dataset$Tmin),
#'     argvar1 = list(df = 5, fun = "ps"),
#'     argvar2 = list(df = 5, fun = "ps"),
#'     arglag = list(df = 5, fun = "ps"),
#'     lag = 14,
#'     group = group
#' )
#'
#' # Fit model
#' model_joint_common <- DLNM_LPS(
#'     dtot ~ 1 + factor(dow) +
#'     sm(month, arg = list(type = "ns", df = 3)) +
#'     sm(year, arg = list(type = "rw1")) +
#'     sm(basis.max_min, arg = list(pen= 2)) +
#'     sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
#'     data = dataset,
#'     offset = log(pop),
#'     model.comparison = list(DIC = TRUE),
#'     family = "Poisson"
#' )
#'
#' # Make predictions
#' pred_max_min <- pred_crossbasis(
#'      basis = basis.max_min,
#'      model = model_joint_common,
#'      at_x = c(29,29,29,30,30,30,31,31,31),
#'      at_y = c(20,21,22,20,21,22,20,21,22),
#'      cen = c(25,15),
#'      lag = c(0,14),
#'      CI = TRUE
#' )
#'
#' # Plot results (contour plot overall RR)
#' contour(pred_max_min, ptype = "overall", cols = function(n) hcl.colors(n, "Inferno"))
#'
#'
#' ### Example 3: Unpenalized and spatially varying DLNM
#'
#' # Create the (unpenalized) crossbasis
#' knots_max <- quantile(dataset$Tmax, probs = c(0.10,0.50,0.90))
#' knots_lag <- dlnm::logknots(0:14, nk = 2)
#'
#' basis.max_unpen <- dlnm::crossbasis(
#'     dataset$Tmax,
#'     argvar=list(fun = "ns", knots = knots_max),
#'     arglag=list(fun="ns", knots = knots_lag),
#'     lag=14,
#'     group=group
#' )
#'
#' # Fit model
#' model_max_spat <- DLNM_LPS(
#'     dtot ~ 1 + factor(dow) +
#'     sm(month, arg = list(type = "ns", df = 3)) +
#'     sm(year, arg = list(type = "rw1")) +
#'     sm(basis.max_unpen, by = District,
#'             arg.by = list(type = "Leroux", map = shapefile_bcn)) +
#'     sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
#'     data = dataset,
#'     offset = log(pop),
#'     model.comparison = list(DIC = TRUE),
#'     family = "Poisson")
#'
#' summary(model_max_spat)
#'
#' # Make predictions
#' pred_spat_max <- pred_crossbasis(
#'     basis.max_unpen,
#'     model_max_spat,
#'     at_x = 22:30,
#'     cen = 25,
#'     lag = c(0,14),
#'     CI = TRUE,
#'     by = unique(dataset$District)
#' )
#'
#' # Plot results (region-specific overal RR)
#' plot(pred_spat_max, ptype = "overall", CI = FALSE, xlab = "temp", ylab = "RR", ylim = c(0.9,3.0),
#'     multiple.plots = FALSE)
#'
#'
#' ### Example 4: Unpenalized and deprivation-varying DLNM
#'
#' # Fit model
#' model_max_conf <- DLNM_LPS(
#'     dtot ~ 1 + factor(dow) +
#'     sm(month, arg = list(type = "ns", df = 3)) +
#'     sm(year, arg = list(type = "rw1")) +
#'     sm(depriv, arg = list(type = "ns", df = 3)) +
#'     sm(basis.max_unpen, by = depriv,
#'             arg.by = list(type = "ns", df = 3)) +
#'     sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
#'     data = dataset,
#'     offset = log(pop),
#'     model.comparison = list(DIC = TRUE),
#'     family = "Poisson")
#'
#' # Make predictions
#' pred_conf_max <- pred_crossbasis(
#'     basis.max_unpen,
#'     model_max_conf,
#'     at_x = 22:30,
#'     cen = 25,
#'     lag = c(0,14),
#'     CI = TRUE,
#'     by = quantile(dataset$depriv,probs= c(0.10, 0.90))
#' )
#'
#' # Extract overall RR for 28 and 30 degrees Celsius for both deprivation indices
#' extract_predictions(
#'     pred_conf_max,
#'     slice_x = c(28,30),
#'     ptype = "allRRfit"
#')
#'
#' # Extract lag-specific RR for 28 and 30 degrees Celsius, lag 0 and 8 and
#' # for 10% deprivation percentile
#' extract_predictions(
#'     pred_conf_max,
#'     slice_x = c(28,30),
#'     slice_lag = c(0,8),
#'     slice_by = quantile(dataset$depriv,probs= c(0.10)),
#'     ptype = "matRRfit"
#')
#'
#' # Plot results (contour plot exposure-lag-response)
#' contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))
#'
#'
#' @export

pred_crossbasis <- function(basis, model, at_x, at_y = NULL, cen, lag = NULL, by.lag = 1, CI.all = FALSE, CI = FALSE,
                   ci.level = 0.95, exc.prob = FALSE, threshold = 0, mode = NULL, by = NULL){


  if(missing(cen)){
    stop("Please specify a centering value")
  }


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
    CI = TRUE
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
  matRRfit <-exp(logpredX)



  if(CI){
    if(is.null(constraint_matrix)){
      sd_logpredX <- Predict_se(prec_xi,Xpred,ind_cb)
    }else{
      sd_logpredX <- Predict_se_con(prec_xi,Xpred,ind_cb, constraint_matrix)
    }

    quantiles <- mapply(function(mean_val, sd_val) {
      stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
    }, logpredX, sd_logpredX)
    Qlower_matfit <- quantiles[1,]
    Qupper_matfit <- quantiles[2,]

    Qlower_matfit <- matrix(Qlower_matfit, nrow = length(predvar), byrow = FALSE)
    Qupper_matfit <- matrix(Qupper_matfit, nrow = length(predvar), byrow = FALSE)

    Qlower_matRRfit <- exp(Qlower_matfit)
    Qupper_matRRfit <- exp(Qupper_matfit)


  }else{
    sd_logpredX <- Qlower_matfit <- Qupper_matfit <- Qlower_matRRfit <- Qupper_matRRfit <- NULL

  }


  # Overall RR
  Xpredall <- 0
  for (j in 1:length(lag)) {
    if(lag[j]==round(lag[j])){
      ind_all <- seq(length(predvar)) + length(predvar) * (j - 1)
      Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }
  }


  pred_all <- as.numeric(Xpredall %*% xi_mode_cb)
  allRRfit <- exp(pred_all)

  if(CI | CI.all){
    if(is.null(constraint_matrix)){
      sd_all <- Predict_se(prec_xi,Xpredall,ind_cb)
    }else{
      sd_all <- Predict_se_con(prec_xi,Xpredall,ind_cb, constraint_matrix)
    }

    quantiles_all <- mapply(function(mean_val, sd_val) {
      stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
    }, pred_all, sd_all)
    Qlower_allfit = quantiles_all[1,]
    Qupper_allfit = quantiles_all[2,]

    Qlower_allRRfit <- exp(Qlower_allfit)
    Qupper_allRRfit <- exp(Qupper_allfit)

  }else{
    sd_all <- Qlower_allfit <- Qupper_allfit <- Qlower_allRRfit <- Qupper_allRRfit <- NULL
  }


  if(exc.prob){
    prob_all <- mapply(function(mean_val, sd_val) {
      stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
    }, pred_all, sd_all)
    exc.prob_all = prob_all

    prob_fit <- mapply(function(mean_val, sd_val) {
      stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
    }, logpredX, sd_logpredX)
    exc.prob_fit = matrix(prob_fit, nrow = length(predvar), byrow = FALSE)

  }else{
    exc.prob_all = NULL
    exc.prob_fit = NULL
  }

  output <- list("matfit" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                 "Qlower_matfit" =  Qlower_matfit,
                 "Qupper_matfit" =  Qupper_matfit,
                 "Qlower_allfit" = Qlower_allfit,
                 "Qupper_allfit" = Qupper_allfit,
                 "allfit" = pred_all,
                 "matRRfit" = matrix(matRRfit, nrow = length(predvar), byrow = FALSE),
                 "Qlower_matRRfit" =  Qlower_matRRfit,
                 "Qupper_matRRfit" =  Qupper_matRRfit,
                 "Qlower_allRRfit" = Qlower_allRRfit,
                 "Qupper_allRRfit" = Qupper_allRRfit,
                 "allRRfit" = allRRfit,
                 "exc.prob_allfit" = exc.prob_all,
                 "exc.prob_matfit" = exc.prob_fit,
                 "x" = at_x,
                 "y" = at_y,
                 "lag" = lag,
                 "family" = model$family)

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
      matRRfit <- exp(logpredX)

      pred_all <- as.numeric(Xpredall %*% c(xi_mode_cb, xi_mode_cb_spat_i))
      allRRfit <- exp(pred_all)

      output[[i]] <- list("matfit" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                          "allfit" = pred_all,
                          "matRRfit" = matrix(matRRfit, nrow = length(predvar), byrow = FALSE),
                          "allRRfit" = allRRfit,
                          "by" = i,
                          "x" = at_x,
                          "y" = at_y,
                          "lag" = lag,
                          "family" = model$family)

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
        }, output[[i]]$allfit, sd_all[i,])
        output[[i]]$Qlower_allfit = quantiles_all[1,]
        output[[i]]$Qupper_allfit = quantiles_all[2,]
        output[[i]]$Qlower_allRRfit = exp(quantiles_all[1,])
        output[[i]]$Qupper_allRRfit = exp(quantiles_all[2,])

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
        }, output[[i]]$matfit, sd_logpredX[i,])
        output[[i]]$Qlower_matfit <- matrix(quantiles[1,], nrow = length(predvar), byrow = FALSE)
        output[[i]]$Qupper_matfit <- matrix(quantiles[2,], nrow = length(predvar), byrow = FALSE)
        output[[i]]$Qlower_matRRfit <- matrix(exp(quantiles[1,]), nrow = length(predvar), byrow = FALSE)
        output[[i]]$Qupper_matRRfit <- matrix(exp(quantiles[2,]), nrow = length(predvar), byrow = FALSE)
      }

      if(exc.prob){
        prob_all <- mapply(function(mean_val, sd_val) {
          stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
        }, output[[i]]$allfit, sd_all[i,])
        output[[i]]$exc.prob_allfit = prob_all

        prob_fit <- mapply(function(mean_val, sd_val) {
          stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
        }, output[[i]]$matfit, sd_logpredX[i,])
        output[[i]]$exc.prob_matfit = matrix(prob_fit, nrow = length(predvar), byrow = FALSE)
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
    matRRfit <- exp(logpredX)

    if( CI ){
      if(is.null(constraint_matrix)){
        sd_logpredX <- Predict_se(prec_xi,Xpred_full,which(index))
      }else{
        sd_logpredX <- Predict_se_con(prec_xi,Xpred_full,which(index), constraint_matrix)
      }

      quantiles <- mapply(function(mean_val, sd_val) {
        stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
      }, logpredX, sd_logpredX)

      Qlower_matfit <- quantiles[1,]
      Qupper_matfit <- quantiles[2,]

      Qlower_matRRfit <- matrix(exp(Qlower_matfit), nrow = length(predvar), byrow = FALSE)
      Qupper_matRRfit <- matrix(exp(Qupper_matfit), nrow = length(predvar), byrow = FALSE)
      Qlower_matfit <- matrix(Qlower_matfit, nrow = length(predvar), byrow = FALSE)
      Qupper_matfit <- matrix(Qupper_matfit, nrow = length(predvar), byrow = FALSE)


    }else{
      Qlower_matfit <- Qupper_matfit <- Qlower_matRRfit <- Qupper_matRRfit <- NULL
    }


    pred_all <- as.numeric(Xpredall %*% mode_cb)
    allRRfit <- exp(pred_all)

    if(CI | CI.all){
      if(is.null(constraint_matrix)){
        sd_all <- Predict_se(prec_xi,Xpredall,which(index))

      }else{
        sd_all <- Predict_se_con(prec_xi,Xpredall,which(index), constraint_matrix)
      }
      quantiles_all <- mapply(function(mean_val, sd_val) {
        stats::qnorm(p = c((1-ci.level)/2, (1+ci.level)/2), mean = mean_val, sd = sd_val)
      }, pred_all, sd_all)
      Qlower_allfit = quantiles_all[1,]
      Qupper_allfit = quantiles_all[2,]
      Qlower_allRRfit = exp(quantiles_all[1,])
      Qupper_allRRfit = exp(quantiles_all[2,])

    }else{
      Qlower_allfit <- Qupper_allfit <- Qlower_allRRfit <- Qupper_allRRfit <- NULL
    }

    if(exc.prob){
      prob_all <- mapply(function(mean_val, sd_val) {
        stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
      }, pred_all, sd_all)
      exc.prob_all = prob_all

      prob_fit <- mapply(function(mean_val, sd_val) {
        stats::pnorm(threshold, mean = mean_val, sd = sd_val, lower.tail = FALSE)
      }, logpredX, sd_logpredX)
      exc.prob_fit = matrix(prob_fit, nrow = length(predvar), byrow = FALSE)

    }else{
      exc.prob_all = NULL
      exc.prob_fit = NULL
    }

    output[[by.value]] <- list("matfit" = matrix(logpredX, nrow = length(predvar), byrow = FALSE),
                   "Qlower_matfit" = Qlower_matfit,
                   "Qupper_matfit" = Qupper_matfit,
                   "Qlower_allfit" = Qlower_allfit,
                   "Qupper_allfit" = Qupper_allfit,
                   "allfit" = pred_all,
                   "matRRfit" = matrix(matRRfit, nrow = length(predvar), byrow = FALSE),
                   "Qlower_matRRfit" =  Qlower_matRRfit,
                   "Qupper_matRRfit" =  Qupper_matRRfit,
                   "Qlower_allRRfit" = Qlower_allRRfit,
                   "Qupper_allRRfit" = Qupper_allRRfit,
                   "allRRfit" = allRRfit,
                   "exc.prob_allfit" = exc.prob_all,
                   "exc.prob_matfit" = exc.prob_fit,
                   "x" = at_x,
                   "y" = at_y,
                   "lag" = lag,
                   "by" = by[by.value],
                   "family" = model$family)
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


#' Construction of data matrix
#' @noRd
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



#' Extract predictions
#'
#' Help function to extract requested predictions from an object returned by [pred_crossbasis()].
#'
#' @param x An object returned by [pred_crossbasis()].
#' @param slice_x Exposure values of the first exposure variable at which to
#' extract predictions. Integer or numeric vector.
#' @param slice_y Exposure values of the second exposure variable at which to
#' extract predictions (bivariate exposures only). Integer or numeric vector.
#' @param slice_lag  Lag value(s) at which to extract predictions. Integer or numeric
#' vector.
#' @param slice_by  Optional argument to specify the levels of the by-variable
#' at which to extract predictions. One value or multiple values supplied through a vector.
#' @param ptype Character string specifying the type of prediction. One of:
#' \describe{
#'   \item{`"matfit"`}{Predictions for specified combinations of exposure and
#'   lag variables.}
#'   \item{`"matRRfit"`}{Predictions for exponentiated associations for
#'   specified combinations of exposure and lag variables.}
#'   \item{`"allfit"`}{Predictions for overall cumulative associations, default value.}
#'   \item{`"allRRfit"`}{Predictions for exponentiated associations for
#'   overall cumulative risks.}
#' }
#'
#' @returns
#' A dataframe containing the requested exposures (and lags), predicted value,
#' lower and upper bounds on the credible intervals (if requested in predictions),
#' exceedance probabilities (if requested in predictions), and levels of the
#' by-variable (if applicable).
#'
#' @examples
#' data(barcelona_temp)
#' data(shapefile_bcn)
#'
#' # Examples of applications
#'
#' dataset <- barcelona_temp |>
#'   dplyr::mutate(
#'     dow = weekdays(date),
#'     year = lubridate::year(date),
#'     month = lubridate::month(date)
#'   )
#'
#' group <- factor(paste0(dataset$District, "-", dataset$year))
#'
#'
#' ### Example 1: Penalized and common DLNM
#'
#' # Create the (unpenalized) crossbasis
#'
#' knots_max <- quantile(dataset$Tmax, probs = c(0.10,0.50,0.90))
#' knots_lag <- dlnm::logknots(0:14, nk = 2)
#'
#' basis.max_unpen <- dlnm::crossbasis(
#'     dataset$Tmax,
#'     argvar=list(fun = "ns", knots = knots_max),
#'     arglag=list(fun="ns", knots = knots_lag),
#'     lag=14,
#'     group=group
#' )
#'
#' # Fit model
#' model_max_conf <- DLNM_LPS(
#'     dtot ~ 1 + factor(dow) +
#'     sm(month, arg = list(type = "ns", df = 3)) +
#'     sm(year, arg = list(type = "rw1")) +
#'     sm(depriv, arg = list(type = "ns", df = 3)) +
#'     sm(basis.max_unpen, by = depriv,
#'             arg.by = list(type = "ns", df = 3)) +
#'     sm(District, arg = list(type = "Leroux", map = shapefile_bcn)),
#'     data = dataset,
#'     offset = log(pop),
#'     model.comparison = list(DIC = TRUE),
#'     family = "Poisson")
#'
#' # Make predictions
#' pred_conf_max <- pred_crossbasis(
#'     basis.max_unpen,
#'     model_max_conf,
#'     at_x = 22:30,
#'     cen = 25,
#'     lag = c(0,14),
#'     CI = TRUE,
#'     by = quantile(dataset$depriv,probs= c(0.10, 0.90))
#' )
#'
#' # Extract overall RR for 28 and 30 degrees Celsius for both deprivation indices
#' extract_predictions(
#'     pred_conf_max,
#'     slice_x = c(28,30),
#'     ptype = "allRRfit"
#')
#'
#' # Extract lag-specific RR for 28 and 30 degrees Celsius, lag 0 and 8 and
#' # for 10% deprivation percentile
#' extract_predictions(
#'     pred_conf_max,
#'     slice_x = c(28,30),
#'     slice_lag = c(0,8),
#'     slice_by = quantile(dataset$depriv,probs= c(0.10)),
#'     ptype = "matRRfit"
#')
#'
#' @export
#'
extract_predictions <- function(x, slice_x = NULL, slice_y = NULL, slice_lag = NULL, slice_by = NULL,
                                ptype = "allfit"){

  nested <- !("allfit" %in% names(x))


  if(nested){
    if(!is.null(slice_by)){
      idx <- which(sapply(x, function(x) which(x$by == slice_by)) != 0)
      x <- x[idx]
    }

    extracted <- list()
    for (xi in 1:length(x)){
      extracted[[as.character(x[[xi]]$by)]] <- extract_predictions_one(x[[xi]],
                                                                       slice_x, slice_y, slice_lag, ptype)
    }

    extracted <- do.call(rbind, Map(
      function(df, name) {
        df <- as.data.frame(df)
        df$by <- name
        df
      },
      extracted,
      names(extracted)
    ))
    rownames(extracted) <- NULL

  }else{
    extracted <- extract_predictions_one(x, slice_x, slice_y, slice_lag, ptype)
  }

  return(extracted)
}


#' Extraction for one single by value
#' @noRd
extract_predictions_one <- function(x, slice_x = NULL, slice_y = NULL, slice_lag = NULL,
                                ptype = "allfit"){


  if(!is.null(slice_y) & is.null(x$y)){
    stop("Can only slice on y if predictions are made for y")
  }
  if(ptype == "allfit"){
    estimate <- x$allfit
    lower <- x$Qlower_allfit
    upper <- x$Qupper_allfit
    exc <- x$exc.prob_allfit
  } else if (ptype == "allRRfit"){
    estimate <- x$allRRfit
    lower <- x$Qlower_allRRfit
    upper <- x$Qupper_allRRfit
    exc <- x$exc.prob_allfit
  } else if (ptype == "matfit"){
    estimate <- x$matfit
    lower <- x$Qlower_matfit
    upper <- x$Qupper_matfit
    exc <- x$exc.prob_matfit
  } else if (ptype == "matRRfit"){
    estimate <- x$matRRfit
    lower <- x$Qlower_matRRfit
    upper <- x$Qupper_matRRfit
    exc <- x$exc.prob_matfit
  }

  if(ptype == "allfit" | ptype == "allRRfit"){
    at_x <- if(is.null(slice_x)){
      seq_len(length(estimate))
    }else{
      which(x$x %in% slice_x)
    }

    at_y <- if(is.null(slice_y)){
      seq_len(length(estimate))
    }else{
      which(x$y %in% slice_y)
    }


    if(length(at_x) == 0){
      stop("Requested x value not found in predictions")
    }

    if(length(at_y) == 0){
      stop("Requested y value not found in predictions")
    }

    if(!is.null(slice_lag)){
      stop("Cannot slice on lag for overall predictions")
    }

    if(!is.null(slice_x) & !is.null(slice_y)){
      at_exp <- mapply(function(xi, yi) {
        which(x$x == xi & x$y == yi)
      }, slice_x, slice_y)
    }else{
      at_exp = intersect(at_x, at_y)
    }


    result <- do.call(cbind, Filter(
      Negate(is.null),
      list(
        x = x$x[at_exp],
        y = if (!is.null(x$y)) x$y[at_exp] else NULL,
        fit = estimate[at_exp],
        lower = if (!is.null(lower)) lower[at_exp] else NULL,
        upper = if (!is.null(upper)) upper[at_exp] else NULL,
        exc.prob = if (!is.null(exc)) exc[at_exp] else NULL
      )
    ))


  }else{
    at_x <- if(is.null(slice_x)){
      seq_len(nrow(estimate))
    }else{
      which(x$x %in% slice_x)
    }

    at_y <- if(is.null(slice_y)){
      seq_len(nrow(estimate))
    }else{
      which(x$y %in% slice_y)
    }


    if(length(at_x) == 0){
      stop("Requested x value not found in predictions")
    }

    if(length(at_y) == 0){
      stop("Requested y value not found in predictions")
    }

    at_lag <- if(is.null(slice_lag)){
      seq_len(ncol(estimate))
    }else{
      which(x$lag %in% slice_lag)
    }
    if(length(at_lag) == 0){
      stop("Requested lag value not found in predictions")
    }

    if(!is.null(slice_x) & !is.null(slice_y)){
      at_exp <- mapply(function(xi, yi) {
        which(x$x == xi & x$y == yi)
      }, slice_x, slice_y)

      result <- expand.grid(x = x$x[at_exp], lag = x$lag[at_lag])
      result$y <- rep(x$y[at_exp], length(at_lag))
      result <- result[,c("x","y","lag")]
    }else{
      at_exp = intersect(at_x, at_y)
      result <- expand.grid(x = x$x[at_exp], lag = x$lag[at_lag])
    }


    result <- do.call(cbind, Filter(
      Negate(is.null),
      list(
        result,
        fit = as.vector(estimate[at_exp, at_lag, drop = FALSE]),
        lower = if (!is.null(lower)) as.vector(lower[at_exp, at_lag, drop = FALSE]) else NULL,
        upper = if (!is.null(upper)) as.vector(upper[at_exp, at_lag, drop = FALSE]) else NULL,
        exc.prob = if (!is.null(exc)) as.vector(exc[at_exp, at_lag, drop = FALSE]) else NULL
      )
    ))

  }



}

