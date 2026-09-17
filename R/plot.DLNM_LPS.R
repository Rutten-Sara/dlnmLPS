#' Plot predicted relative risks
#'
#' Plot predicted relative risks obtained from [predRR()]. Depending on
#' `ptype`, the function produces overall exposure-response curves,
#' exposure-lag-response surfaces, or slices through the exposure-lag-response
#' surface.
#'
#' @param x An object returned by [predRR()].
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
#' @param CI Logical. Should confidence intervals be displayed (defaults = FALSE)? Confidence
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
#' @seealso [contour.DLNMLPS_pred()],[predRR()]
#' @export

plot.DLNMLPS_pred <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall", CI = FALSE,
                          xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                          ylab = NULL, zlab = NULL, ...){

  mar = c(5.1, 4.1, 4.1, 2.1)

  nested <- !("pred_all" %in% names(x))

  if(nested){
    graphics::par( mar = mar, mfrow = if(length(x)==1) c(1,1) else
      if(length(x)==2) c(1,2) else
        if(length(x)==3) c(1,3) else
          c(2,2))

    for (xi in 1:length(x)){
      DLNMLPS_plot_ind(x[[xi]], slice_lag, slice_x , slice_y , ptype , CI ,
                       xlim , ylim , zlim ,cols , xlab ,
                       ylab , zlab , main = paste('By = ', x[[xi]]$by), ...)
    }


    graphics::par(mar = mar, mfrow = c(1,1))

  }else{
    graphics::par(mar = mar, mfrow = c(1,1))
    DLNMLPS_plot_ind(x, slice_lag , slice_x , slice_y, ptype , CI ,
     xlim , ylim , zlim ,cols , xlab ,
     ylab , zlab, ...)
  }
}


#' Plotting one curve of nested graph
#' @noRd
DLNMLPS_plot_ind <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall", CI = FALSE,
                              xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                              ylab = NULL, zlab = NULL, ...){


  if(ptype == "slices"){
    if(is.null(x$y) | sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) == 2 ){
      plot_slices2D(x, slice_lag, slice_x, slice_y, CI, xlim, ylim, cols, xlab, ylab, ...)
    } else if (!is.null(x$y) & sum(!sapply(list(slice_lag, slice_x, slice_y), is.null))==1){
      plot_slices3D(x, slice_lag, slice_x, slice_y, CI, xlim, ylim, zlim, cols, xlab, ylab, zlab,FALSE, ...)
    } else{
      stop("Specify one or two slicing variables")
    }
  } else if (ptype == "overall"){
    if(is.null(x$y)){
      plot_overall_uni(x, xlim, ylim, cols, xlab, ylab, CI, ...)
      if(sum(!sapply(list(slice_lag, slice_x, slice_y), is.null)) >0){
        warning("No slices for overall RR with univariate exposure")
      }
    }else{
      if(sum(sapply(list(slice_lag, slice_x, slice_y), is.null)) == 3){
        plot_overall_3D(x, xlim, ylim, zlim, cols, xlab, ylab, zlab, CI, FALSE, ...)
      } else if(!is.null(slice_x) | !is.null(slice_y)){
        plot_overall_2D(x, slice_x, slice_y, xlim, ylim, cols, xlab, ylab, CI, ...)
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



#' Contour plots of predicted relative risks
#'
#' Produce contour plots of predicted relative risks obtained from
#' [predRR()]. In contrast to [plot.DLNMLPS_pred()], this function only
#' produces contour plots of surfaces and therefore is only available for
#' plot configurations that result in a three-dimensional surface.
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
#' @seealso [plot.DLNMLPS_pred()], [predRR()]
#' @export

contour.DLNMLPS_pred <- function(x, slice_lag = NULL, slice_x = NULL, slice_y = NULL, ptype = "overall",
                          xlim = NULL, ylim = NULL, zlim = NULL,cols = NULL, xlab = NULL,
                          ylab = NULL, zlab = NULL, ...){

  mar = c(5.1, 4.1, 4.1, 2.1)

  nested <- !("pred_all" %in% names(x))

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
plot_overall_uni <- function(x, xlim, ylim, cols, xlab, ylab, CI, ... ){

  if(is.null(xlab)){
    xlab = "exposure"
  }

  if(is.null(ylab)){
    ylab = "RR"
  }

  if (is.null(xlim))
    xlim <- range(unlist(x$x), na.rm = TRUE)


  if (is.null(ylim)) {

    if (CI) {
      ylim <- range(c(
        x$pred_all,
        x$Qlower_all,
        x$Qupper_all
      ), na.rm = TRUE)
    } else {
      ylim <- range(unlist(x$pred_all), na.rm = TRUE)
    }
  }

  if(is.null(cols)){
    cols = "steelblue"
  }


  plot(
    x$x,
    x$pred_all,
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
      x$pred_all,
      if (CI) x$Qlower_all,
      if (CI) x$Qupper_all,
      cols,
      CI
    )




}



#' Plotting overall cumulative RR curve for two-exposure DLNM (in 2D)
#' @noRd
plot_overall_2D <- function(x, slice_x, slice_y, xlim, ylim, cols, xlab, ylab, CI, ... ){


  if(is.null(ylab)){
    ylab = "RR"
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

      plot_y[[i]] <- x$pred_all[idx_x]


      if (CI) {

        plot_lower[[i]] <- x$Qlower_all[idx_x]
        plot_upper[[i]] <- x$Qupper_all[idx_x]

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

      plot_y[[i]] <- x$pred_all[idx_x]


      if (CI) {

        plot_lower[[i]] <- x$Qlower_all[idx_x]
        plot_upper[[i]] <- x$Qupper_all[idx_x]

      }
    }


    plot_labels <- paste0(
      "y=", slice_y
    )

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


  graphics::legend(
    "top",
    legend = plot_labels,
    col = cols,
    lty = 1,
    bty = "n"
  )




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
    zlab = "RR"
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

    z[ix, iy] <- x$pred_all[i]
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
    zlab = "RR"
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

    z[ix, iy] <- exp(x$logpredX[iy,ix])
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
    zlab = "RR"
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

      z <- t(exp(x$logpredX[idx, ]))

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

      z <- t(exp(x$logpredX[idx, ]))

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

        z[ix, iy] <- exp(x$logpredX[i, idx_lag])
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
                        CI = TRUE, xlim = NULL, ylim = NULL,
                        cols = NULL, xlab = NULL, ylab = NULL, ...) {


  if(is.null(ylab)){
    ylab = "RR"
  }


  if (is.null(slice_lag) && is.null(slice_x) && is.null(slice_y)) {
    stop("Specify a slicing variable")
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

      plot_y <- exp(x$logpredX[, idx, drop = FALSE])

      plot_y <- lapply(seq_len(ncol(plot_y)),
                       function(i) plot_y[, i])
      plot_x <- rep(list(x$x), length(plot_y))

      if (CI) {
        plot_lower <- exp(x$Qlower_logpredX[, idx, drop = FALSE])
        plot_upper <- exp(x$Qupper_logpredX[, idx, drop = FALSE])

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

      plot_y <- exp(x$logpredX[idx, , drop = FALSE])
      plot_y <- lapply(seq_len(nrow(plot_y)),
                       function(i) plot_y[i, ])
      plot_x <- rep(list(x$lag), length(plot_y))


      if (CI) {
        plot_lower <- exp(x$Qlower_logpredX[idx, , drop = FALSE])
        plot_upper <- exp(x$Qupper_logpredX[idx, , drop = FALSE])

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

      plot_y <- exp(x$logpredX[idx, , drop = FALSE])
      plot_y <- lapply(seq_len(nrow(plot_y)),
                       function(i) plot_y[i, ])
      plot_x <- rep(list(x$lag), length(plot_y))

      if (CI) {
        plot_lower <- exp(x$Qlower_logpredX[idx, , drop = FALSE])
        plot_upper <- exp(x$Qupper_logpredX[idx, , drop = FALSE])

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

        plot_y[[i]] <- exp(
          x$logpredX[idx_x, idx_lag]
        )


        if (CI) {

          plot_lower[[i]] <- exp(
            x$Qlower_logpredX[idx_x, idx_lag]
          )

          plot_upper[[i]] <- exp(
            x$Qupper_logpredX[idx_x, idx_lag]
          )
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

        plot_y[[i]] <- exp(
          x$logpredX[idx_y, idx_lag]
        )


        if (CI) {

          plot_lower[[i]] <- exp(
            x$Qlower_logpredX[idx_y, idx_lag]
          )

          plot_upper[[i]] <- exp(
            x$Qupper_logpredX[idx_y, idx_lag]
          )
        }
      }


      plot_labels <- paste0(
        "y=", slice_y,
        ", lag=", slice_lag
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


  graphics::legend(
    "top",
    legend = plot_labels,
    col = cols,
    lty = 1,
    bty = "n"
  )
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





