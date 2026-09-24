#' Plot predicted outcomes
#'
#' Plot predicted outcomes obtained from [pred_crossbasis()]. For Poisson or
#' negative binomial models with a log link, relative risks (RRs) are plotted;
#' for binomial models with a logit link, odds ratios (ORs) are plotted; and
#' for Gaussian models with an identity link, additive differences in the
#' outcome are plotted. Depending on `ptype`, the function produces overall
#' exposure-response curves, exposure-lag-response surfaces, or slices
#' through the exposure-lag-response surface.
#'
#' @param x An object returned by [pred_crossbasis()].
#' @param slice_lag Lag value(s) at which to obtain slices. Integer or numeric
#' vector.
#' @param slice_x Exposure values of the first exposure variable at which to
#' obtain slices. Integer or numeric vector.
#' @param slice_y Exposure values of the second exposure variable at which to
#' obtain slices (bivariate exposures only). Integer or numeric vector.
#' @param ptype Character string specifying the type of plot. One of:
#' \describe{
#'   \item{`"overall"`}{Overall cumulative exposure-response plot.}
#'   \item{`"cross"`}{Complete exposure-lag-response surface (only available for
#'   univariate exposures).}
#'   \item{`"slices"`}{Slices through the exposure-lag-response surface.}
#' }
#' @param multiple.plots Logical. Should the results of different by-levels be plotted in
#' different windows (default = TRUE)?
#' @param show.legend Logical. Should the legend appear on the plot (default = TRUE)?
#' @param CI Logical. Should confidence intervals be displayed (default = FALSE)? Confidence
#' intervals are only available for two-dimensional plots.
#' @param xlim,ylim,zlim Numeric vectors of length two giving axis limits.
#' @param cols Colours used for lines, surfaces or contours.
#' @param xlab,ylab,zlab Axis labels.
#' @param ... Additional graphical arguments passed to the underlying plotting
#' functions.
#'
#' @details
#' The available plot types depend on whether the prediction object contains
#' univariate or bivariate exposures.
#'
#' **`ptype = "overall"`**
#'
#' For univariate exposures, the overall cumulative exposure-response curve is
#' plotted. No slicing variables may be specified.
#'
#' For bivariate exposures, the overall cumulative exposure-response is shown as
#' a three-dimensional surface if no slicing variable is specified. Alternatively,
#' either `slice_x` or `slice_y` (but not both) may be specified to obtain a
#' two-dimensional cumulative exposure-response curve.
#'
#' **`ptype = "cross"`**
#'
#' Plots the exposure-lag-response surface for univariate exposures. This plot
#' type is not available for bivariate exposures.
#'
#' **`ptype = "slices"`**
#'
#' Plots slices through the exposure-lag-response surface.
#'
#' For univariate exposures, at least one of `slice_lag` or `slice_x` must be
#' specified.
#'
#' For bivariate exposures, at least one of `slice_lag`, `slice_x` or
#' `slice_y` must be specified. If only one slicing variable is specified, the
#' result is a three-dimensional surface. If two slicing variables are specified,
#' the result is a two-dimensional curve.
#'
#' Three-dimensional plots can display only a single slice at a time. In
#' contrast, two-dimensional plots may contain multiple slices, resulting in
#' multiple curves in the same figure.
#'
#' Confidence intervals are only available for two-dimensional plots.
#'
#' @returns
#' Produces a plot. The function is called for its side effects and returns
#' `NULL` invisibly.
#'
#' @seealso [contour.DLNMLPS_pred()],[pred_crossbasis()]
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
#' # Plot results (contour plot exposure-lag-response)
#' contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))
#'
#'
#' @export

plot.DLNMLPS_pred <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall",
                              multiple.plots = TRUE, show.legend = TRUE, CI = FALSE,
                              xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                              ylab = NULL, zlab = NULL, ...){

  mar = c(5.1, 4.1, 4.1, 2.1)

  nested <- !("allfit" %in% names(x))

  if(nested){
    if(multiple.plots){
    graphics::par( mar = mar, mfrow = if(length(x)==1) c(1,1) else
      if(length(x)==2) c(1,2) else
        if(length(x)==3) c(1,3) else
          c(2,2))

    for (xi in 1:length(x)){
      DLNMLPS_plot_ind(x[[xi]], slice_lag, slice_x , slice_y , ptype , FALSE, show.legend, CI ,
                       xlim , ylim , zlim ,cols , xlab ,
                       ylab , zlab , main = paste('By = ', x[[xi]]$by), ...)
    }


    graphics::par(mar = mar, mfrow = c(1,1))
    }else{
      graphics::par(mar = mar, mfrow = c(1,1))
      DLNMLPS_plot_ind(x, slice_lag , slice_x , slice_y, ptype , TRUE, show.legend, CI ,
                       xlim , ylim , zlim ,cols , xlab ,
                       ylab , zlab, ...)
    }

  }else{
    graphics::par(mar = mar, mfrow = c(1,1))
    DLNMLPS_plot_ind(x, slice_lag , slice_x , slice_y, ptype , FALSE, show.legend, CI ,
     xlim , ylim , zlim ,cols , xlab ,
     ylab , zlab, ...)
  }
}


#' Plotting one curve of nested graph
#' @noRd
DLNMLPS_plot_ind <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall",
                             multiple.plots = FALSE, show.legend = TRUE, CI = FALSE,
                             xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                             ylab = NULL, zlab = NULL, ...){


  if(ptype == "slices"){
    if(is.null(x$y) | sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) == 2 ){
      plot_slices2D(x, slice_lag, slice_x, slice_y, multiple.plots, show.legend, CI, xlim, ylim, cols, xlab, ylab, ...)
    } else if (!is.null(x$y) & sum(!sapply(list(slice_lag, slice_x, slice_y), is.null))==1){
      plot_slices3D(x, slice_lag, slice_x, slice_y, CI, xlim, ylim, zlim, cols, xlab, ylab, zlab,FALSE, ...)
    } else{
      stop("Specify one or two slicing variables")
    }
  } else if (ptype == "overall"){
    if(is.null(x$y)){
      plot_overall_uni(x, xlim, ylim, cols, xlab, ylab,multiple.plots, show.legend, CI, ...)
      if(sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) >0){
        warning("No slices for overall RR with univariate exposure")
      }
    }else{
      if(sum(sapply(list(slice_lag, slice_x, slice_y), is.null)) == 3){
        plot_overall_3D(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, CI, FALSE, ...)
      } else if(!is.null(slice_x) | !is.null(slice_y)){
        plot_overall_2D(x, slice_x, slice_y, xlim, ylim, cols, xlab, ylab, multiple.plots,show.legend,CI, ...)
      } else{
        stop("Overall RR cannot be calculated for sliced lag")
      }
    }

  } else if (ptype == "cross"){
    if(!is.null(x$y)){
      stop("Cross cannot be specified for bivariate exposures (use slices or overall instead)")
    }else{
      plot_cross_3D(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, CI, FALSE, ...)
      if(sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) >0){
        warning("No slices for cross exposure-lag-response surface with univariate exposure")
      }
    }
  }else{
    stop("ptype should be slices, cross or overall")
  }

}



#' Contour plots of predicted outcomes
#'
#' Produce contour plots of predicted outcomes obtained from
#' [pred_crossbasis()]. In contrast to [plot.DLNMLPS_pred()], this function only
#' produces contour plots of surfaces and therefore is only available for
#' plot configurations that result in a three-dimensional surface. For Poisson or
#' negative binomial models with a log link, relative risks (RRs) are plotted;
#' for binomial models with a logit link, odds ratios (ORs) are plotted; and
#' for Gaussian models with an identity link, additive differences in the
#' outcome are plotted.
#'
#'
#' @inheritParams plot.DLNMLPS_pred
#'
#' @details
#' This function provides contour representations of the three-dimensional
#' surfaces available in [plot.DLNMLPS_pred()].
#'
#' The available plots are:
#' \itemize{
#'   \item the exposure-lag-response surface (`ptype = "cross"`) for
#'   univariate exposures;
#'   \item the overall cumulative exposure-response surface
#'   (`ptype = "overall"`) for bivariate exposures when no exposure slice is
#'   specified;
#'   \item slice surfaces (`ptype = "slices"`) whenever the selected slicing
#'   variables define a three-dimensional surface.
#' }
#'
#' Contour plots are not available for two-dimensional plots (e.g.,
#' cumulative exposure-response curves or exposure/lag slices shown as
#' multiple lines). Confidence intervals are therefore not supported.
#'
#' @returns
#' Produces a contour plot and returns `NULL` invisibly.
#'
#' @seealso [plot.DLNMLPS_pred()], [pred_crossbasis()]
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
#' ### Example 1: Unpenalized and spatially varying DLNM
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
#' plot(pred_spat_max, ptype = "overall", CI = TRUE, xlab = "temp", ylab = "RR", ylim = c(0.9,3.0))
#'
#'
#' ### Example 2: Unpenalized and deprivation-varying DLNM
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
#' # Plot results (contour plot exposure-lag-response)
#' contour(pred_conf_max, ptype = "cross", cols = function(n) hcl.colors(n, "Inferno"))
#'
#'
#' ### Example 3: Penalized and common DLNM
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
#' ### Example 4: Penalized and common joint DLNM
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
#' @export

contour.DLNMLPS_pred <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall",
                          xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                          ylab = NULL, zlab = NULL, ...){

  mar = c(5.1, 4.1, 4.1, 2.1)

  nested <- !("allfit" %in% names(x))

  if(nested){
    graphics::par(mar = mar, mfrow = if(length(x)==1) c(1,1) else
      if(length(x)==2) c(1,2) else
        if(length(x)==3) c(1,3) else
          c(2,2))

    for (xi in 1:length(x)){
      DLNMLPS_contour_ind(x[[xi]], slice_lag, slice_x , slice_y , ptype ,
                       xlim , ylim , zlim ,cols , xlab,
                       ylab , zlab, main = paste('By = ', x[[xi]]$by), ...)
    }

    graphics::par(mar = mar, mfrow = c(1,1))
  }else{
    graphics::par(mar = mar, mfrow = c(1,1))

    DLNMLPS_contour_ind(x, slice_lag , slice_x, slice_y , ptype  ,
                     xlim, ylim , zlim  ,cols , xlab  ,
                     ylab , zlab , ...)
  }


}




#' Plotting one contour plot of nested graph
#' @noRd
DLNMLPS_contour_ind <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall",
                                 xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                                 ylab = NULL, zlab = NULL, ...){

  if(ptype == "slices"){
    if(is.null(x$y) ){
      stop("Contour plot of sliced RR only possible for bivariate exposures")
    } else if (sum(!sapply(list(slice_lag, slice_x, slice_y), is.null))==1){
      plot_slices3D(x, slice_lag, slice_x, slice_y, FALSE, xlim, ylim, zlim, cols, xlab, ylab, zlab, TRUE, ...)
    } else{
      stop("Specify one slicing variables")
    }
  } else if (ptype == "overall"){
    if(is.null(x$y)){
      stop("Contour plot of overall RR only possible for bivariate exposures")
    }else{
      if(sum(sapply(list(slice_lag, slice_x, slice_y), is.null)) == 3){
        plot_overall_3D(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, FALSE, TRUE, ...)
      } else{
        warning("No slicing variables used for contour plot of overall RR")
      }
    }

  } else if (ptype == "cross"){
    if(!is.null(x$y)){
      stop("Cross cannot be specified for bivariate exposures (use slices or overall instead)")
    }else{
      plot_cross_3D(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, FALSE, TRUE, ...)
      if(sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) >0){
        warning("No slices for cross exposure-lag-response surface with univariate exposure")
      }
    }
  }else{
    stop("ptype should be slices, overall or cross")
  }
}



#' Plotting overall cumulative RR curve for one-exposure DLNM
#' @noRd
plot_overall_uni <- function(x, xlim, ylim, cols, xlab, ylab, multiple.plots,show.legend, CI, ... ){

  if(is.null(xlab)){
    xlab = "exposure"
  }



  if (multiple.plots) {
    ncurve <- length(x)

    if (is.null(cols))
      cols <- grDevices::hcl.colors(ncurve, "Dark 3")

    if (is.null(xlim)) {
      xlim <- range( unlist(lapply(x, function(z) z$x)), na.rm = TRUE )
    }


    if (is.null(ylim)) {
      if (x[[1]]$family != "Gaussian") {
        vals <- unlist(lapply(x, function(z) {
          if (CI) {
            c(z$allRRfit, z$Qlower_allRRfit, z$Qupper_allRRfit)
          } else {
            z$allRRfit
          }
        }))
      } else {
          vals <- unlist(lapply(x, function(z) {
            if (CI) {
              c(z$allfit, z$Qlower_allfit, z$Qupper_allfit)
            } else {
              z$allfit
            }
          }))
      }

      ylim <- range(vals, na.rm = TRUE)

    }

    if (is.null(ylab)) {
      if (x[[1]]$family %in% c("Poisson", "NB")) {
        ylab <- "RR"
      } else if (x[[1]]$family == "Binomial") {
        ylab <- "OR"
      } else {
        ylab <- "difference in outcome"
      }
    }


    # ---- First plot ----
    plot_labels <- paste("by = ",  unlist(unname(lapply(x, function(x_by) x_by$by))))
    z <- x[[1]]
    if (z$family != "Gaussian") {
      plot(
        z$x,
        z$allRRfit,
        type = "l",
        col = cols[1],
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        ...
      )

      for (i in seq_len(ncurve)) {
        z <- x[[i]]
        draw_slice2D(
          z$x,
          z$allRRfit,
          if (CI) z$Qlower_allRRfit,
          if (CI) z$Qupper_allRRfit,
          cols[i],
          CI
        )
      }



    }else{
      plot(
        z$x,
        z$allfit,
        type = "l",
        col = cols[1],
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        ...
      )


        for (i in seq_len(ncurve)) {
          z <- x[[i]]
          draw_slice2D(
            z$x,
            z$allfit,
            if (CI) z$Qlower_allfit,
            if (CI) z$Qupper_allfit,
            cols[i],
            CI
          )
      }
    }



  if(show.legend){
    graphics::legend(
      "top",
      legend = plot_labels,
      col = cols,
      lty = 1,
      bty = "n"
    )
  }


  }else{


  if(is.null(ylab)){
    if(x$family == "Poisson" | x$family == "NB"){
      ylab = "RR"
    }else if (x$family == "Binomial"){
      ylab = "OR"
    }else{
      ylab = "difference in outcome"
    }
  }

  if (is.null(xlim))
    xlim <- range(unlist(x$x), na.rm = TRUE)


  if (is.null(ylim)) {

    if (CI) {
      if(x$family != "Gaussian"){
        ylim <- range(c(
          x$allRRfit,
          x$Qlower_allRRfit,
          x$Qupper_allRRfit
        ), na.rm = TRUE)
      }else{
      ylim <- range(c(
        x$allfit,
        x$Qlower_allfit,
        x$Qupper_allfit
      ), na.rm = TRUE)
    }} else {
      if(x$family != "Gaussian"){
      ylim <- range(unlist(x$allRRfit), na.rm = TRUE)
      }else{
        ylim <- range(unlist(x$allfit), na.rm = TRUE)
      }
    }
  }

  if(is.null(cols)){
    cols = "steelblue"
  }


  if(x$family != "Gaussian"){
  plot(
    x$x,
    x$allRRfit,
    type = "l",
    col = cols,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  draw_slice2D(
      x$x,
      x$allRRfit,
      if (CI) x$Qlower_allRRfit,
      if (CI) x$Qupper_allRRfit,
      cols,
      CI
    )
  }else{
    plot(
      x$x,
      x$allfit,
      type = "l",
      col = cols,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      ...
    )

    draw_slice2D(
      x$x,
      x$allfit,
      if (CI) x$Qlower_allfit,
      if (CI) x$Qupper_allfit,
      cols,
      CI
    )
}
}


}



#' Plotting overall cumulative RR curve for two-exposure DLNM (in 2D)
#' @noRd
plot_overall_2D <- function(x, slice_x, slice_y, xlim, ylim, cols, xlab, ylab, multiple.plots,
                            show.legend, CI, ... ){

  if(multiple.plots & (length(slice_x)>1 | length(slice_y)>1)){
    warning("Multiple plots for different by-values required when slice is a vector")
    multiple.plots = FALSE
  }

  if(multiple.plots){
    ncurve <- length(x)


    if (is.null(ylab)) {
      if (x[[1]]$family %in% c("Poisson", "NB")) {
        ylab <- "RR"
      } else if (x[[1]]$family == "Binomial") {
        ylab <- "OR"
      } else {
        ylab <- "difference in outcome"
      }
    }

  plot_labels <- paste("by = ",  unlist(unname(lapply(x, function(x_by) x_by$by))))



  plot_x <- vector("list", ncurve)
  plot_y <- vector("list", ncurve)

  if (CI) {
    plot_lower <- vector("list", ncurve)
    plot_upper <- vector("list", ncurve)
  }

  if(!is.null(slice_x)){

    if(is.null(xlab)){
      xlab = "y-exposure"
    }


    for (i in seq_len(ncurve)) {
      z <- x[[i]]
      idx_x <- which(z$x == slice_x)

      if (length(idx_x) == 0 )
        stop("Some slicing values are not present in predictions.")

      plot_x[[i]] <- z$y[idx_x]

      if(z$family != "Gaussian"){
      plot_y[[i]] <- z$allRRfit[idx_x]
      }else{
        plot_y[[i]] <- z$allfit[idx_x]
      }


      if (CI) {
        if(z$family != "Gaussian"){
        plot_lower[[i]] <- z$Qlower_allRRfit[idx_x]
        plot_upper[[i]] <- z$Qupper_allRRfit[idx_x]
        }else{
          plot_lower[[i]] <- z$Qlower_allfit[idx_x]
          plot_upper[[i]] <- z$Qupper_allfit[idx_x]
        }

      }
    }

  }else{

    if(is.null(xlab)){
      xlab = "x-exposure"
    }



    for (i in seq_len(ncurve)) {
      z <- x[[i]]
      idx_x <- which(z$y == slice_y[i])

      if (length(idx_x) == 0 )
        stop("Some slice values are not present in predictions.")

      plot_x[[i]] <- z$x[idx_x]

      if(z$family != "Gaussian"){
        plot_y[[i]] <- z$allRRfit[idx_x]
      }else{
        plot_y[[i]] <- z$allfit[idx_x]
      }


      if (CI) {
        if(z$family != "Gaussian"){
        plot_lower[[i]] <- z$Qlower_allRRfit[idx_x]
        plot_upper[[i]] <- z$Qupper_allRRfit[idx_x]
        }else{
          plot_lower[[i]] <- z$Qlower_allfit[idx_x]
          plot_upper[[i]] <- z$Qupper_allfit[idx_x]
        }

      }
    }


}


  }else{
    if(is.null(ylab)){
      if(x$family == "Poisson" | x$family == "NB"){
        ylab = "RR"
      }else if (x$family == "Binomial"){
        ylab = "OR"
      }else{
        ylab = "difference in outcome"
      }
    }


    if(!is.null(slice_x)){

      if(is.null(xlab)){
        xlab = "y-exposure"
      }


      ncurve <- length(slice_x)

      plot_x <- vector("list", ncurve)
      plot_y <- vector("list", ncurve)

      if (CI) {
        plot_lower <- vector("list", ncurve)
        plot_upper <- vector("list", ncurve)
      }

      for (i in seq_len(ncurve)) {

        idx_x <- which(x$x == slice_x[i])

        if (length(idx_x) == 0 )
          stop("Some slicing values are not present in predictions.")

        plot_x[[i]] <- x$y[idx_x]

        if(x$family != "Gaussian"){
          plot_y[[i]] <- x$allRRfit[idx_x]
        }else{
          plot_y[[i]] <- x$allfit[idx_x]
        }


        if (CI) {
          if(x$family != "Gaussian"){
            plot_lower[[i]] <- x$Qlower_allRRfit[idx_x]
            plot_upper[[i]] <- x$Qupper_allRRfit[idx_x]
          }else{
            plot_lower[[i]] <- x$Qlower_allfit[idx_x]
            plot_upper[[i]] <- x$Qupper_allfit[idx_x]
          }

        }
      }

      plot_labels <- paste0(
        "x=", slice_x
      )

    }else{

      if(is.null(xlab)){
        xlab = "x-exposure"
      }

      ncurve <- length(slice_y)

      plot_x <- vector("list", ncurve)
      plot_y <- vector("list", ncurve)

      if (CI) {
        plot_lower <- vector("list", ncurve)
        plot_upper <- vector("list", ncurve)
      }

      for (i in seq_len(ncurve)) {

        idx_x <- which(x$y == slice_y[i])

        if (length(idx_x) == 0 )
          stop("Some slice values are not present in predictions.")

        plot_x[[i]] <- x$x[idx_x]

        if(x$family != "Gaussian"){
          plot_y[[i]] <- x$allRRfit[idx_x]
        }else{
          plot_y[[i]] <- x$allfit[idx_x]
        }


        if (CI) {
          if(x$family != "Gaussian"){
            plot_lower[[i]] <- x$Qlower_allRRfit[idx_x]
            plot_upper[[i]] <- x$Qupper_allRRfit[idx_x]
          }else{
            plot_lower[[i]] <- x$Qlower_allfit[idx_x]
            plot_upper[[i]] <- x$Qupper_allfit[idx_x]
          }

        }
      }


      plot_labels <- paste0(
        "y=", slice_y
      )

    }

}


  # -------------------------------
  # Plot
  # -------------------------------


  ncurve <- length(plot_y)


  if (is.null(cols))
    cols <- grDevices::hcl.colors(ncurve, "Dark 3")


  if (is.null(xlim))
    xlim <- range(unlist(plot_x), na.rm = TRUE)


  if (is.null(ylim)) {

    if (CI) {
      ylim <- range(c(
        unlist(plot_y),
        unlist(plot_lower),
        unlist(plot_upper)
      ), na.rm = TRUE)
    } else {
      ylim <- range(unlist(plot_y), na.rm = TRUE)
    }
  }


  plot(
    plot_x[[1]],
    plot_y[[1]],
    type = "l",
    col = cols[1],
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )


  for (i in seq_len(ncurve)) {

    draw_slice2D(
      plot_x[[i]],
      plot_y[[i]],
      if (CI) plot_lower[[i]],
      if (CI) plot_upper[[i]],
      cols[i],
      CI
    )
  }


  if(show.legend){
  graphics::legend(
    "top",
    legend = plot_labels,
    col = cols,
    lty = 1,
    bty = "n"
  )
  }

}

#' Plotting overall cumulative RR curve for two-exposure DLNM (in 3D)
#' @noRd
plot_overall_3D <- function(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, CI, contour = FALSE, ... ){


  if(is.null(xlab)){
    xlab = "x-exposure"
  }

  if(is.null(ylab)){
    ylab = "y-exposure"
  }

  if(is.null(zlab)){
    if(x$family == "Poisson" | x$family == "NB"){
      zlab = "RR"
    }else if (x$family == "Binomial"){
      zlab = "OR"
    }else{
      zlab = "difference in outcome"
    }
  }


  nx <- length(unique(x$x))
  ny <- length(unique(x$y))

  #if (nx * ny != length(x$x)) {
  #  stop("3D plots require a complete x-y prediction grid.")
  #}



  xvals <- sort(unique(x$x))
  yvals <- sort(unique(x$y))

  z <- matrix(
    NA_real_,
    nrow = length(xvals),
    ncol = length(yvals),
    dimnames = list(xvals, yvals)
  )

  for (i in seq_along(x$x)) {
    ix <- match(x$x[i], xvals)
    iy <- match(x$y[i], yvals)

    if(x$family != "Gaussian"){
    z[ix, iy] <- x$allRRfit[i]
    }else{
      z[ix, iy] <- x$allfit[i]
    }
  }



  if (CI) {
    warning("CI not possible for 3D surface")
  }


  if (is.null(xlim))
    xlim <- range(unlist(xvals), na.rm = TRUE)


  if (is.null(ylim))
    ylim <- range(unlist(yvals), na.rm = TRUE)

  if (is.null(zlim))
    zlim <- range(unlist(z), na.rm = TRUE)


  if(!contour){
    if(is.null(cols)){
      cols = "steelblue"
    }

    graphics::persp(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      col=cols,
      ...
    )
  }else{
    if(is.null(cols)){
      cols = function(n) grDevices::hcl.colors(n, "Viridis")
    }

    graphics::filled.contour(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      color.palette = cols,
      ...
    )
  }


}



#' Plotting exposure-lag-response curve
#' @noRd
plot_cross_3D <- function(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, CI, contour = FALSE, ... ){

  if(is.null(xlab)){
    xlab = "lag"
  }

  if(is.null(ylab)){
    ylab = "exposure"
  }

  if(is.null(zlab)){
    if(x$family == "Poisson" | x$family == "NB"){
      zlab = "RR"
    }else if (x$family == "Binomial"){
      zlab = "OR"
    }else{
      zlab = "difference in outcome"
    }
  }

  nx <- length(unique(x$lag))
  ny <- length(unique(x$x))


  xvals <- x$lag
  yvals <- sort(unique(x$x))

  z <- matrix(
    NA_real_,
    nrow = length(xvals),
    ncol = length(yvals),
    dimnames = list(xvals, yvals)
  )

  for (i in seq_along(x$x)) {
    for (j in seq_along(x$lag)){
    ix <- match(x$lag[j], xvals)
    iy <- match(x$x[i], yvals)

    if(x$family != "Gaussian"){
    z[ix, iy] <- exp(x$matRRfit[iy,ix])
    }else{
      z[ix, iy] <- exp(x$matfit[iy,ix])
    }
    }
  }



  if (CI) {
    warning("CI not possible for 3D surface")
  }


  if (is.null(xlim))
    xlim <- range(unlist(xvals), na.rm = TRUE)


  if (is.null(ylim))
    ylim <- range(unlist(yvals), na.rm = TRUE)

  if (is.null(zlim))
    zlim <- range(unlist(z), na.rm = TRUE)


  # -------------------------------
  # Plot
  # -------------------------------
  if(!contour){
    if(is.null(cols)){
      cols = "steelblue"
    }

    graphics::persp(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      col=cols,
      ...
    )
  }else{
    if(is.null(cols)){
      cols = function(n) grDevices::hcl.colors(n, "Viridis")
    }

    graphics::filled.contour(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      color.palette = cols,
      ...
    )
  }




}


#' Plotting sliced RR-curve (in 3D)
#' @noRd
plot_slices3D <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL,
                          CI = TRUE, xlim = NULL, ylim = NULL, zlim = NULL,
                          cols = NULL, xlab = NULL, ylab = NULL, zlab = NULL, contour = FALSE, ...) {


  if(is.null(zlab)){
    if(x$family == "Poisson" | x$family == "NB"){
      zlab = "RR"
    }else if (x$family == "Binomial"){
      zlab = "OR"
    }else{
      zlab = "difference in outcome"
    }
  }

  if (is.matrix(x$x)) {
    stop("Slices can only be plotted for a vector of constant exposures")
  }


  # -------------------------------
  # Extract curves
  # -------------------------------


    if (!is.null(slice_x)) {

      if(is.null(xlab)){
        xlab = "lag"
      }


      if(is.null(ylab)){
        ylab = "y-exposure"
      }

      if(length(slice_x)>1){
        stop("Only one x slice can be specified with 3D plot of two-dimensional exposures")
      }

      idx <- which(x$x == slice_x)

      xvals <- x$lag
      yvals <- x$y[idx]

      ord <- order(yvals)

      if (length(idx)==0)
        stop("Slice x must be present in predictions.")

      if(x$family != "Gaussian"){
      z <- t(exp(x$matRRfit[idx, ]))
      }else{
        z <- t(exp(x$matfit[idx, ]))
      }

      yvals <- yvals[ord]
      z <- z[, ord , drop = FALSE]


    } else if (!is.null(slice_y)) {

      if(is.null(xlab)){
        xlab = "lag"
      }


      if(is.null(ylab)){
        ylab = "x-exposure"
      }

      if(length(slice_y)>1){
        stop("Only one y slice can be specified with 3D plot of two-dimensional exposures")
      }

      idx <- which(x$y == slice_y)

      xvals <- x$lag
      yvals <- x$x[idx]

      ord <- order(yvals)

      if (length(idx)==0)
        stop("Slice y must be present in predictions.")

      if(x$family != "Gaussian"){
        z <- t(exp(x$matRRfit[idx, ]))
      }else{
        z <- t(exp(x$matfit[idx, ]))
      }

      yvals <- yvals[ord]
      z <- z[, ord , drop = FALSE]

      }else if (!is.null(slice_lag)) {

        if(is.null(xlab)){
          xlab = "x-exposure"
        }

        if(is.null(ylab)){
          ylab = "y-exposure"
        }


      if(length(slice_lag)>1){
        stop("Only one lag slice can be specified with 3D plot of two-dimensional exposures")
      }

      idx_lag <- match(slice_lag, x$lag)


      if (is.na(idx_lag))
        stop("Slice lag must be present in predictions.")


      nx <- length(unique(x$x))
      ny <- length(unique(x$y))

      #if (nx * ny != length(x$x)) {
      #  stop("3D plots require a complete x-y prediction grid.")
      #}

      xvals <- sort(unique(x$x))
      yvals <- sort(unique(x$y))


      z <- matrix(
        NA_real_,
        nrow = length(xvals),
        ncol = length(yvals),
        dimnames = list(xvals, yvals)
      )

      for (i in seq_along(x$x)) {
        ix <- match(x$x[i], xvals)
        iy <- match(x$y[i], yvals)

        if(x$family != "Gaussian"){
          z[ix, iy]  <-exp(x$matRRfit[i, idx_lag])
        }else{
          z[ix, iy] <- exp(x$matfit[i, idx_lag])
        }

      }


    }


  #if (anyNA(z))
  #   stop("Incomplete x-y grid.")

  if (CI) {
    warning("CI not possible for 3D surface")
  }


  if (is.null(xlim))
    xlim <- range(unlist(xvals), na.rm = TRUE)


  if (is.null(ylim))
    ylim <- range(unlist(yvals), na.rm = TRUE)

  if (is.null(zlim))
    zlim <- range(unlist(z), na.rm = TRUE)



  # -------------------------------
  # Plot
  # -------------------------------
  if(!contour){
    if(is.null(cols)){
      cols = "steelblue"
    }

    graphics::persp(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      col=cols,
      ...
    )
  }else{
    if(is.null(cols)){
      cols = function(n) grDevices::hcl.colors(n, "Viridis")
    }

    graphics::filled.contour(
      xvals,
      yvals,
      z,
      xlab = xlab,
      ylab = ylab,
      color.palette = cols,
      ...
    )
  }
}






#' Plotting sliced RR-curve (in 2D)
#' @noRd
plot_slices2D <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL,
                        multiple.plots = TRUE, show.legend = TRUE,
                        CI = TRUE, xlim = NULL, ylim = NULL,
                        cols = NULL, xlab = NULL, ylab = NULL, ...) {

  if(multiple.plots & (length(slice_x)>1 | length(slice_y)>1 | length(slice_lag)>1)){
    warning("Multiple plots for different by-values required when slice is a vector")
    multiple.plots = FALSE
  }

  if (is.null(slice_lag) && is.null(slice_x) && is.null(slice_y)) {
    stop("Specify a slicing variable")
  }

  if(multiple.plots){
    ncurve <- length(x)


    if (is.null(ylab)) {
      if (x[[1]]$family %in% c("Poisson", "NB")) {
        ylab <- "RR"
      } else if (x[[1]]$family == "Binomial") {
        ylab <- "OR"
      } else {
        ylab <- "difference in outcome"
      }
    }

    plot_labels <- paste("by = ",  unlist(unname(lapply(x, function(x_by) x_by$by))))

    if (is.matrix(x[[1]]$x)) {
      stop("Slices can only be plotted for a vector of constant exposures")
    }


    plot_x <- vector("list", ncurve)
    plot_y <- vector("list", ncurve)

    if (CI) {
      plot_lower <- vector("list", ncurve)
      plot_upper <- vector("list", ncurve)
    }


    if (is.null(x[[1]]$y)) {

      if (!is.null(slice_y)) {
        stop("Can only slice on y if predictions are made for y")
      }
      if (!is.null(slice_lag)) {

        if(is.null(xlab)){
          xlab = "exposure"
        }

        idx <- match(slice_lag, x[[1]]$lag)

        if (anyNA(idx))
          stop("Slice lag must be present in predictions.")

        for (i in seq_len(ncurve)) {
        z = x[[i]]
        if(z$family != "Gaussian"){
          plot_y[[i]] <- exp(z$matRRfit[, idx, drop = TRUE])
        }else{
          plot_y[[i]] <- exp(z$matfit[, idx, drop = TRUE])
        }
        plot_x[[i]] <- z$x


        if (CI) {
          if(z$family != "Gaussian"){
            plot_lower[[i]] <- exp(z$Qlower_matRRfit[, idx, drop = TRUE])
            plot_upper[[i]] <- exp(z$Qupper_matRRfit[, idx, drop = TRUE])
          }else{
            plot_lower[[i]] <- exp(z$Qlower_matfit[, idx, drop = TRUE])
            plot_upper[[i]] <- exp(z$Qupper_matfit[, idx, drop = TRUE])
          }

        }


        }

      }else{
        if(is.null(xlab)){
          xlab = "lag"
        }


        idx <- unlist(lapply(slice_x, function(s) which(x[[1]]$x == s)))

        if (anyNA(idx))
          stop("Slice exposure must be present in predictions")

        for (i in seq_len(ncurve)) {
          z = x[[i]]
          if(z$family != "Gaussian"){
            plot_y[[i]] <- exp(z$matRRfit[idx, , drop = TRUE])
          }else{
            plot_y[[i]] <- exp(z$matfit[idx, , drop = TRUE])
          }
          plot_x[[i]] <- z$lag


          if (CI) {
            if(z$family != "Gaussian"){
              plot_lower[[i]] <- exp(z$Qlower_matRRfit[idx, , drop = TRUE])
              plot_upper[[i]] <- exp(z$Qupper_matRRfit[idx, , drop = TRUE])
            }else{
              plot_lower[[i]] <- exp(z$Qlower_matfit[idx, , drop = TRUE])
              plot_upper[[i]] <- exp(z$Qupper_matfit[idx, , drop = TRUE])
            }

          }

        }



      }

    }else{
      if (!is.null(slice_x) && !is.null(slice_y)) {


        if(is.null(xlab)){
          xlab = "lag"
        }

        idx <- mapply(function(xi, yi) {
          which(x[[1]]$x == xi & x[[1]]$y == yi)
        }, slice_x, slice_y)

        if (any(lengths(idx) == 0))
          stop("Some x-y slices are not present in predictions.")

        idx <- unlist(idx)


        for (i in seq_len(ncurve)) {
          z <- x[[i]]
        if(z$family != "Gaussian"){
          plot_y[[i]] <- exp(z$matRRfit[idx, , drop = TRUE])
        }else{
          plot_y[[i]] <- exp(z$matfit[idx, , drop = TRUE])
        }


         plot_x[[i]] <- z$lag

        if (CI) {
          if(z$family != "Gaussian"){
            plot_lower <- exp(z$Qlower_matRRfit[idx, , drop = TRUE])
            plot_upper <- exp(z$Qupper_matRRfit[idx, , drop = TRUE])
          }else{
            plot_lower <- exp(z$Qlower_matfit[idx, , drop = TRUE])
            plot_upper <- exp(z$Qupper_matfit[idx, , drop = TRUE])
          }

        }

        }

      } else if (!is.null(slice_x) && !is.null(slice_lag)) {


        if(is.null(xlab)){
          xlab = "y-exposure"
        }


        if (length(slice_x) != length(slice_lag))
          stop("slice_x and slice_lag must have the same length.")


        idx_x <- which(x[[1]]$x == slice_x)
        idx_lag <- match(slice_lag, x[[1]]$lag)


        if (length(idx_x) == 0 || is.na(idx_lag))
          stop("Some slice values are not present in predictions.")


        for (i in seq_len(ncurve)) {
          z <- x[[i]]


          # free dimension is y
          plot_x[[i]] <- z$y[idx_x]

          if(z$family != "Gaussian"){
            plot_y[[i]] <- exp(
              z$matRRfit[idx_x, idx_lag]
            )
          }else{
            plot_z[[i]] <- exp(
              z$matfit[idx_x, idx_lag]
            )
          }


          if (CI) {
            if(z$family != "Gaussian"){
              plot_lower[[i]] <- exp(
                z$Qlower_matRRfit[idx_x, idx_lag]
              )

              plot_upper[[i]] <- exp(
                z$Qupper_matRRfit[idx_x, idx_lag]
              )
            }else{
              plot_lower[[i]] <- exp(
                z$Qlower_matfit[idx_x, idx_lag]
              )

              plot_upper[[i]] <- exp(
                z$Qupper_matfit[idx_x, idx_lag]
              )
            }
          }
        }


      }else if (!is.null(slice_y) && !is.null(slice_lag)) {


        if(is.null(xlab)){
          xlab = "x-exposure"
        }


        if (length(slice_y) != length(slice_lag))
          stop("slice_y and slice_lag must have the same length.")


        idx_y <- which(x[[1]]$y == slice_y)
        idx_lag <- match(slice_lag, x[[1]]$lag)

        if (length(idx_y) == 0 || is.na(idx_lag))
          stop("Some slice values are not present in predictions.")

        for (i in seq_len(ncurve)) {
          z <- x[[i]]

          # free dimension is x
          plot_x[[i]] <- z$x[idx_y]

          if(z$family != "Gaussian"){
            plot_y[[i]] <- exp(
              z$matRRfit[idx_y, idx_lag]
            )
          }else{
            plot_y[[i]] <- exp(
              z$matfit[idx_y, idx_lag]
            )
          }


          if (CI) {
            if(z$family != "Gaussian"){
              plot_lower[[i]] <- exp(
                z$Qlower_matRRfit[idx_y, idx_lag]
              )

              plot_upper[[i]] <- exp(
                z$Qupper_matRRfit[idx_y, idx_lag]
              )
            }else{
              plot_lower[[i]] <- exp(
                z$Qlower_matfit[idx_y, idx_lag]
              )

              plot_upper[[i]] <- exp(
                z$Qupper_matfit[idx_y, idx_lag]
              )
            }
          }
        }

      }
    }



  }else{

  if(is.null(ylab)){
    if(x$family == "Poisson" | x$family == "NB"){
      ylab = "RR"
    }else if (x$family == "Binomial"){
      ylab = "OR"
    }else{
      ylab = "difference in outcome"
    }
  }



  if (is.matrix(x$x)) {
    stop("Slices can only be plotted for a vector of constant exposures")
  }


  # -------------------------------
  # Extract curves
  # -------------------------------

  if (is.null(x$y)) {

    if (!is.null(slice_y)) {
      stop("Can only slice on y if predictions are made for y")
    }

    if (!is.null(slice_lag)) {

      if(is.null(xlab)){
        xlab = "exposure"
      }


      idx <- match(slice_lag, x$lag)

      if (anyNA(idx))
        stop("Slice lag must be present in predictions.")

      if(x$family != "Gaussian"){
      plot_y <- exp(x$matRRfit[, idx, drop = FALSE])
      }else{
        plot_y <- exp(x$matfit[, idx, drop = FALSE])
      }

      plot_y <- lapply(seq_len(ncol(plot_y)),
                       function(i) plot_y[, i])
      plot_x <- rep(list(x$x), length(plot_y))

      if (CI) {
        if(x$family != "Gaussian"){
        plot_lower <- exp(x$Qlower_matRRfit[, idx, drop = FALSE])
        plot_upper <- exp(x$Qupper_matRRfit[, idx, drop = FALSE])
        }else{
          plot_lower <- exp(x$Qlower_matfit[, idx, drop = FALSE])
          plot_upper <- exp(x$Qupper_matfit[, idx, drop = FALSE])
        }
        plot_lower <- lapply(seq_len(ncol(plot_lower)),
                             function(i) plot_lower[, i])

        plot_upper <- lapply(seq_len(ncol(plot_upper)),
                             function(i) plot_upper[, i])
      }

      plot_labels <- paste0("lag = ", slice_lag)

    } else {

      if(is.null(xlab)){
        xlab = "lag"
      }


      idx <- unlist(lapply(slice_x, function(s) which(x$x == s)))

      if (anyNA(idx))
        stop("Slice exposure must be present in predictions")

      if(x$family != "Gaussian"){
      plot_y <- exp(x$matRRfit[idx, , drop = FALSE])
      }else{
        plot_y <- exp(x$matfit[idx, , drop = FALSE])
      }
      plot_y <- lapply(seq_len(nrow(plot_y)),
                       function(i) plot_y[i, ])
      plot_x <- rep(list(x$lag), length(plot_y))


      if (CI) {
        if(x$family != "Gaussian"){
        plot_lower <- exp(x$Qlower_matRRfit[idx, , drop = FALSE])
        plot_upper <- exp(x$Qupper_matRRfit[idx, , drop = FALSE])
        }else{
          plot_lower <- exp(x$Qlower_matfit[idx, , drop = FALSE])
          plot_upper <- exp(x$Qupper_matfit[idx, , drop = FALSE])
        }

        plot_lower <- lapply(seq_len(nrow(plot_lower)),
                             function(i) plot_lower[i, ])

        plot_upper <- lapply(seq_len(nrow(plot_upper)),
                             function(i) plot_upper[i, ])
      }

      plot_labels <- paste0("x = ", slice_x)
    }


  } else {


    # -------------------------------
    # Two-dimensional exposure
    # -------------------------------


    if (!is.null(slice_x) && !is.null(slice_y)) {


      if(is.null(xlab)){
        xlab = "lag"
      }

      idx <- mapply(function(xi, yi) {
        which(x$x == xi & x$y == yi)
      }, slice_x, slice_y)

      if (any(lengths(idx) == 0))
        stop("Some x-y slices are not present in predictions.")

      idx <- unlist(idx)

      if(x$family != "Gaussian"){
        plot_y <- exp(x$matRRfit[idx, , drop = FALSE])
      }else{
        plot_y <- exp(x$matfit[idx, , drop = FALSE])
      }

      plot_y <- lapply(seq_len(nrow(plot_y)),
                       function(i) plot_y[i, ])
      plot_x <- rep(list(x$lag), length(plot_y))

      if (CI) {
        if(x$family != "Gaussian"){
        plot_lower <- exp(x$Qlower_matRRfit[idx, , drop = FALSE])
        plot_upper <- exp(x$Qupper_matRRfit[idx, , drop = FALSE])
        }else{
          plot_lower <- exp(x$Qlower_matfit[idx, , drop = FALSE])
          plot_upper <- exp(x$Qupper_matfit[idx, , drop = FALSE])
        }

        plot_lower <- lapply(seq_len(nrow(plot_lower)),
                             function(i) plot_lower[i, ])

        plot_upper <- lapply(seq_len(nrow(plot_upper)),
                             function(i) plot_upper[i, ])
      }

      plot_labels <- paste0(
        "x=", slice_x,
        ", y=", slice_y
      )


    } else if (!is.null(slice_x) && !is.null(slice_lag)) {


      if(is.null(xlab)){
        xlab = "y-exposure"
      }


      if (length(slice_x) != length(slice_lag))
        stop("slice_x and slice_lag must have the same length.")


      ncurve <- length(slice_x)

      plot_x <- vector("list", ncurve)
      plot_y <- vector("list", ncurve)

      if (CI) {
        plot_lower <- vector("list", ncurve)
        plot_upper <- vector("list", ncurve)
      }


      for (i in seq_len(ncurve)) {

        idx_x <- which(x$x == slice_x[i])
        idx_lag <- match(slice_lag[i], x$lag)

        if (length(idx_x) == 0 || is.na(idx_lag))
          stop("Some slice values are not present in predictions.")


        # free dimension is y
        plot_x[[i]] <- x$y[idx_x]

        if(x$family != "Gaussian"){
        plot_y[[i]] <- exp(
          x$matRRfit[idx_x, idx_lag]
        )
        }else{
          plot_y[[i]] <- exp(
            x$matfit[idx_x, idx_lag]
          )
        }


        if (CI) {
          if(x$family != "Gaussian"){
          plot_lower[[i]] <- exp(
            x$Qlower_matRRfit[idx_x, idx_lag]
          )

          plot_upper[[i]] <- exp(
            x$Qupper_matRRfit[idx_x, idx_lag]
          )
          }else{
            plot_lower[[i]] <- exp(
              x$Qlower_matfit[idx_x, idx_lag]
            )

            plot_upper[[i]] <- exp(
              x$Qupper_matfit[idx_x, idx_lag]
            )
          }
        }
      }


      plot_labels <- paste0(
        "x=", slice_x,
        ", lag=", slice_lag
      )


    } else if (!is.null(slice_y) && !is.null(slice_lag)) {


      if(is.null(xlab)){
        xlab = "x-exposure"
      }


      if (length(slice_y) != length(slice_lag))
        stop("slice_y and slice_lag must have the same length.")


      ncurve <- length(slice_y)

      plot_x <- vector("list", ncurve)
      plot_y <- vector("list", ncurve)

      if (CI) {
        plot_lower <- vector("list", ncurve)
        plot_upper <- vector("list", ncurve)
      }


      for (i in seq_len(ncurve)) {

        idx_y <- which(x$y == slice_y[i])
        idx_lag <- match(slice_lag[i], x$lag)

        if (length(idx_y) == 0 || is.na(idx_lag))
          stop("Some slice values are not present in predictions.")


        # free dimension is x
        plot_x[[i]] <- x$x[idx_y]

        if(x$family != "Gaussian"){
        plot_y[[i]] <- exp(
          x$matRRfit[idx_y, idx_lag]
        )
        }else{
          plot_y[[i]] <- exp(
            x$matfit[idx_y, idx_lag]
          )
        }


        if (CI) {
          if(x$family != "Gaussian"){
          plot_lower[[i]] <- exp(
            x$Qlower_matRRfit[idx_y, idx_lag]
          )

          plot_upper[[i]] <- exp(
            x$Qupper_matRRfit[idx_y, idx_lag]
          )
          }else{
            plot_lower[[i]] <- exp(
              x$Qlower_matfit[idx_y, idx_lag]
            )

            plot_upper[[i]] <- exp(
              x$Qupper_matfit[idx_y, idx_lag]
            )
          }
        }
      }


      plot_labels <- paste0(
        "y=", slice_y,
        ", lag=", slice_lag
      )

    }
  }

}

  # -------------------------------
  # Plot
  # -------------------------------


  ncurve <- length(plot_y)


  if (is.null(cols))
    cols <- grDevices::hcl.colors(ncurve, "Dark 3")


  if (is.null(xlim))
    xlim <- range(unlist(plot_x), na.rm = TRUE)


  if (is.null(ylim)) {

    if (CI) {
      ylim <- range(c(
        unlist(plot_y),
        unlist(plot_lower),
        unlist(plot_upper)
      ), na.rm = TRUE)
    } else {
      ylim <- range(unlist(plot_y), na.rm = TRUE)
    }
  }


  plot(
    plot_x[[1]],
    plot_y[[1]],
    type = "l",
    col = cols[1],
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )


  for (i in seq_len(ncurve)) {

    draw_slice2D(
      plot_x[[i]],
      plot_y[[i]],
      if (CI) plot_lower[[i]],
      if (CI) plot_upper[[i]],
      cols[i],
      CI
    )
  }


  if(show.legend){
  graphics::legend(
    "top",
    legend = plot_labels,
    col = cols,
    lty = 1,
    bty = "n"
  )
  }
}


#' Plotting credible intervals
#' @noRd
fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- utils::modifyList(list(col = grDevices::grey(0.9), border = NA),
                            ci.arg)
  polygon.arg <- utils::modifyList(polygon.arg, list(x = c(x,
                                                    rev(x)), y = c(high, rev(low))))
  do.call(graphics::polygon, polygon.arg)
}


#' Draws different curves on a graph
#' @noRd
draw_slice2D <- function(x, y, lower = NULL, upper = NULL,
                       col, CI = TRUE) {

  if (CI) {
    graphics::lines(
      x,
      y,
      col = col,
      lwd = 1.5
    )

    fci(
      x = x,
      high = upper,
      low = lower,
      ci.arg = list(
        col = grDevices::adjustcolor(col, alpha.f = 0.2)
      ),
      plot.arg = list(
        type = "l",
        col = col,
        lwd = 1.5
      )
    )
  } else {
    graphics::lines(
      x,
      y,
      col = col,
      lwd = 1.5
    )
  }
}





