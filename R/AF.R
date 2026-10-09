#' Attributable risk from DLNMLPS model
#'
#' This function calculated the attributable risk (attributable fraction and
#' attributable number) from a fitted [DLNM_LPS] model.
#'
#' @param basis An object of class `"crossbasis"` or `"joint.crossbasis"`
#'   defining the exposure-lag basis used in the model.
#' @param cases A vector of cases or, only if `dir = "forw"`, a matrix of
#'   future cases.
#' @param model An object of class `"DLNMLPS"` returned by the [DLNM_LPS] function.
#' @param cen A value  (or a vector of length 2 for `"joint.crossbasis"`)
#'   specifying the centering (reference) exposure level used
#'   for calculating relative risk predictions.
#' @param x Vector of the first exposure or, only if `dir = "back"` a matrix
#'   of lagged exposures, for which the attributable fraction needs to be
#'   calculated.
#' @param y Vector of the second exposure or, only if `dir = "back"` a matrix
#'   of lagged exposures, for which the attributable fraction needs to be
#'   calculated. Only for `basis` of class`"joint.crossbasis"`.
#' @param by Vector of the `by` variable-levels (same length as cases) at
#'   which the DLNM relationship is evaluated. For a continuous `by` variable,
#'   this is a vector of numeric values defining the level at
#'   which the effect is evaluated. For a factor `by` variable,
#'   this is a vector of factor levels defining the corresponding stratum.
#' @param group Vector of region indicators to specify an optional grouping structure
#'   of the output. For spatially varying DLNMs, this argument automatically corresponds
#'   to the ID vector in `by` (i.e., region-specific attributable risks are calculated).
#'   For non-spatially varying DLNMs, a vector can also be specified (same length
#'   as cases) to group the output per region.
#' @param dir Either "back" (default) or "forw" for backward and forward
#'   perspective on attributable risk respectively.
#' @param tot If `True`, the total attributable risk over the entire period is calculated.
#' @param range_x The range of the first exposure for which the attributable risk is
#'   calculated. If `NULL`(default), the whole range is used.
#' @param range_y The range of the second exposure for which the attributable risk is
#'   calculated. If `NULL`(default), the whole range is used.
#'   Only for `basis` of class`"joint.crossbasis"`.
#' @param samples Set to `TRUE` if samples from the posterior distribution should be returned.
#'   Samples can only be calculated for `tot = TRUE`.
#' @param nsamples Number of samples from the posterior distribution to be returned.
#'   (only for `nsamples = TRUE`).
#' @param exc.threshold Numeric value defining the threshold used for exceedance
#'   probability calculations of the attributable fraction. If `NULL` (default), no
#'   exceedance probabilities are calculated.
#'
#' @returns A list containing the attributable fractions and numbers.
#'
#'   If `samples = TRUE`, the list is structured in two main elements.
#'   Note that this option is only possible with `tot = TRUE`. The two elements are:
#' \describe{
#'   \item{samples element}{List with samples from the posterior distribution of the total
#'   attributable fraction `af` and total attributable number `an`. The group `ID`
#'   (if applicable) and exceedance probabilities (if requested) are also reported.
#'   For a non-NULL grouping structures, this output is structured per group.}
#'   \item{fit element}{The fitted total `af`, `an`, and the group `ID`
#'   if applicable. For a non-NULL grouping structures, this output is structured per group. }
#' }
#'
#'   If `samples = FALSE`, the list contains the fitted total `af`, `an`, and the group `ID`
#'   if applicable, when `tot = TRUE`. When `tot = FALSE`, this list contains the day-specific
#'   `af`, `an` and `cases`, as well as the group `ID` if applicable. For non-NULL grouping
#'   structures, this output is structured per group.
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
#' # Calculate attributable fraction (forward perspective) for 2007
#'
#' data2007_district <- dataset |>
#'   dplyr::filter(year == 2007)
#'
#' # Calculate attributable risk
#' af2007_common <- attributable_risk(
#'   basis = basis.max,
#'   cases = data2007_district$dtot,
#'   model = model_max_common,
#'   cen = 25,
#'   x = data2007_district$Tmax,
#'   group = data2007_district$District, # Calculate af per district
#'   dir = "forw",
#'   tot = TRUE,
#'   samples = TRUE  # Samples from posterior distribution
#' )
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
#' # Calculate attributable fraction (backward perspective) for 2007
#'
#' af2007_joint <- attributable_risk(
#'   basis = basis.max_min,
#'   cases = data2007_district$dtot,
#'   model = model_joint_common,
#'   cen = c(25,15),
#'   x = data2007_district$Tmax,
#'   y = data2007_district$Tmin,
#'   group = data2007_district$District, # Calculate af per district
#'   dir = "back",
#'   tot = TRUE,
#'   samples = TRUE  # Samples from posterior distribution
#' )
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
#'
#' # Calculate attributable fraction (backward perspective) for 2007
#'
#' af2007_spat <- attributable_risk(
#'   basis = basis.max_unpen,
#'   cases = data2007_district$dtot,
#'   model = model_max_spat,
#'   cen = 25,
#'   x = data2007_district$Tmax,
#'   by = data2007_district$District, # Spatially varying DLNM by district
#'   dir = "back",
#'   tot = TRUE,
#'   samples = TRUE  # Samples from posterior distribution
#' )
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
#' # Calculate attributable fraction (forward perspective) for 2007
#'
#' af2007_conf <- attributable_risk(
#'   basis = basis.max_unpen,
#'   cases = data2007_district$dtot,
#'   model = model_max_conf,
#'   cen = 25,
#'   x = data2007_district$Tmax,
#'   by = data2007_district$depriv, # By-value for every district
#'   group = data2007_district$District, # Group output by district
#'   dir = "forw",
#'   tot = TRUE,
#'   samples = TRUE  # Samples from posterior distribution
#' )
#'
#'
#' @export
#'
attributable_risk <- function(basis, cases, model, cen, x, y = NULL, by = NULL, group = NULL,
                              dir ="back", tot = TRUE, range_x = NULL,
                              range_y = NULL, samples = FALSE, nsamples = 1000,
                              exc.threshold = NULL){

  dir <- match.arg(dir,c("back","forw"))

  if(missing(cen)){
    stop("Please specify a centering value")
  }

  if (!is.null(exc.threshold) & !samples){
    exc.threshold = NULL
    warning("Exceedance probabilities can only be calculated with samples = TRUE")
  }


  # Which type of basis?
  rhs <-  model$model[[3]]
  all_terms <- split_terms(rhs)

  index_term <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                      gsub("\\s+", "", all_terms),
                      fixed = TRUE)

  kind_term <- model$types[which(index_term)-1]

  if(kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
    group = by
  }

  if(kind_term %in% c("spatial DLNM", "joint spatial DLNM","modified DLNM", "joint modified DLNM") &&
     is.null(by)){
     stop("Please specify the level of the by-variable")
  }


  # Select range
  if("joint.crossbasis" %in% class(basis)){
    if(is.null(y)) stop("Specify the values of the y-exposure")
    if(!is.null(range_x) & is.null(range_y)){
      x[x<range_x[1]|x>range_x[2]] <- cen[1]
      y[x<range_x[1]|x>range_x[2]] <- cen[2]
    } else if(is.null(range_x) & !is.null(range_y)){
      x[y<range_y[1]|y>range_y[2]] <- cen[1]
      y[y<range_y[1]|y>range_y[2]] <- cen[2]
    } else if(!is.null(range_x) & !is.null(range_y)){
      x[x<range_x[1]|x>range_x[2]|y<range_y[1]|y>range_y[2]] <- cen[1]
      y[x<range_x[1]|x>range_x[2]|y<range_y[1]|y>range_y[2]] <- cen[2]
    }
  }else{
    if(!is.null(range_x)) x[x<range_x[1]|x>range_x[2]] <- cen
  }


  # Construct the lagged exposure
  #   - If dir="back" -> lagged exposures
  #   - If dir="forw" -> constant exposures
  lag <- attr(basis,"lag")
  if("joint.crossbasis" %in% class(basis)){
  if(NCOL(x)==1L & NCOL(y)==1L) {
    at_x <- if(dir=="back") tsModel::Lag(x,seq(lag[1],lag[2]),group=group) else
      matrix(rep(x,diff(lag)+1),length(x))
    at_y <- if(dir=="back") tsModel::Lag(y,seq(lag[1],lag[2]), group=group) else
      matrix(rep(y,diff(lag)+1),length(y))
  } else if(NCOL(x)!=1L & NCOL(y)!=1L) {
    if(dir=="forw") stop("'x' and 'y' must be a vector of dimension 1 when dir='forw'")
    if(ncol(at_x <- x)!=diff(lag)+1)
      stop("dimension of 'x' not compatible with 'basis'")
    if(ncol(at_y <- y)!=diff(lag)+1)
      stop("dimension of 'y' not compatible with 'basis'")
  } else{stop("dimension of x and y are not compatible")}
  }else{
    if(NCOL(x)==1L) {
      at_x <- if(dir=="back") tsModel::Lag(x,seq(lag[1],lag[2]),group=group) else
        matrix(rep(x,diff(lag)+1),length(x))
    } else{
      if(dir=="forw") stop("'x' must be a vector of dimension 1 when dir='forw'")
      if(ncol(at_x <- x)!=diff(lag)+1)
        stop("dimension of 'x' not compatible with 'basis'")
    }
  }


  # Number of cases
  # - if backwards: cases must be a vector
  # - if forward:
  #     - If cases matrix (today, tomorrow...) -> row average
  #     - If cases vector -> compute forward moving average yourself
  # Also compute denominator of total attributable number and fraction

  if(NROW(cases)!=NROW(at_x)) stop("'exposure' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    if(is.null(group)){
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
    }else{
      den <- tapply(1:nrow(cases), group,
                    function(idx) sum(rowMeans(cases[idx, , drop=FALSE], na.rm = TRUE), na.rm = TRUE))
      cases <- rowMeans(cases)
    }
  } else {
    if(!is.null(group)){
      den <- tapply(1:length(cases), group,
                    function(idx) sum(cases[idx], na.rm = TRUE))
    }else{
      den <- sum(cases,na.rm=TRUE)
    }
    if(dir=="forw")
      cases <- rowMeans(as.matrix(tsModel::Lag(cases,-seq(lag[1],lag[2]), group=group)))
  }

  # Extract the model
    index <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                   gsub("\\s+", "", colnames(model$Xpred)),
                   fixed = TRUE)

    coef <- model$xi_mode
    Prec <- model$Prec
    model.link <- ifelse(model$family %in% c("Poisson", "NB"), "log",
                         ifelse(model$family == "Binomial", "logit", "identity"))
    if(!model.link %in% c("log","logit"))
      stop("'model' must have a log or logit link function")

    constraint_matrix <- model$constraint_matrix


    # Basis transformation

    if("joint.crossbasis" %in% class(basis)){

      predvar <- c(seq(NROW(at_x)),"cen")
      predlag <- seq(lag[1],lag[2], by = 1)

      # Data matrix for prediction
      Xpred <- data_pred_2d(basis, at_x, at_y, predlag, predvar, cen)

      predvar <- predvar[-length(predvar)]

    }else{

      predvar <- seq(NROW(at_x))
      predlag <- seq(lag[1],lag[2], by = 1)

      # Data matrix for prediction
      Xpred <- data_pred_1d(basis, at_x, predlag, cen)

    }


    dim_basis <- ncol(basis)
    ind_cb <- which(index)[1:dim_basis]
    xi_mode_cb = coef[ind_cb]

    # Effect-modified basis: interaction with effect-modification
    if (kind_term %in% c("modified DLNM", "joint modified DLNM")){
      type_smooth <- model$smooth[[which(index_term)-1]]
      output = list()

      if("ns" %in% class(type_smooth)){
          z = as.matrix(splines::ns(by, knots = attributes(type_smooth)$knots,
                                    Boundary.knots = attributes(type_smooth)$Boundary.knots,
                                    intercept = attributes(type_smooth)$intercept))
       }else if("bs" %in% class(type_smooth)){
          z = as.matrix(splines::bs(by, knots = attributes(type_smooth)$knots,
                                    Boundary.knots = attributes(type_smooth)$Boundary.knots,
                                    intercept = attributes(type_smooth)$intercept,
                                    degree = attributes(type_smooth)$degree))
      } else if ("ps" %in% class(type_smooth)){
          z = as.matrix(dlnm::ps(by, knots = attributes(type_smooth)$knots,
                                 df = attributes(type_smooth)$df,
                                 degree = attributes(type_smooth)$degree,
                                 intercept = attributes(type_smooth)$intercept))
      } else if (!is.null(attributes(type_smooth)$Dimnames[[1]])){
        z = NULL
        for(by.value in 1:length(by)){
          index_by <- grepl(by[by.value], attributes(type_smooth)$Dimnames[[2]],
                     fixed = TRUE)
          z_part <-rep(0,length(attributes(type_smooth)$Dimnames[[2]]))
          z_part[index_by] =1
          z <- rbind(z,matrix(z_part, nrow = 1))
        }
      } else{
          z = matrix(by, ncol = 1)
        }



      at_inter <- Matrix::kronecker(Matrix::Matrix(1, nrow = length(predlag), ncol = 1), z)

      Mat_crossbasis_spat <- Matrix::Matrix(Xpred, sparse = TRUE) # cross basis

      interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
        at_inter * Mat_crossbasis_spat[, j]
      })

      interaction_matrix <- do.call(cbind, interaction_list)

      # Global design mtarix
      Xpred_full <- Matrix::Matrix(cbind(Xpred,interaction_matrix), sparse = TRUE)


      }else if (kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
        Xpred_full <- cbind(Xpred, Xpred)
      } else{
        Xpred_full <- Xpred
      }


    # OVerall matrix
    Xpredall <- 0
    for (j in seq(length = length(predlag))) {
      ind_all <- seq(length(predvar)) + length(predvar) * (j - 1) # first lag period for all observations
      Xpredall <- Xpredall + Xpred_full[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }



    # Compute AF and AN

    if(kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
      data_list <- data.frame(Xpredall,cases,group)
      dlist <- split(data_list, data_list$group)
      output_list = list()

      for(i in unique(names(dlist))){
        ind_spat_ID <- grepl(i, colnames(model$Xpred),
                             fixed = TRUE)

        ind_spat <- index & ind_spat_ID

        xi_mode_cb_spat_i = coef[ind_spat]
        ##########################

        # COMPUTE AF AND AN
        ind_dlist <-  !grepl("^(cases|group)", names(dlist[[i]]))
        af <- 1-exp(-drop(as.matrix(as.matrix(dlist[[i]][,ind_dlist])%*%c(xi_mode_cb, xi_mode_cb_spat_i))))
        ind_dlist_cases <-grepl("cases", names(dlist[[i]]))

        an <- af*dlist[[i]][,ind_dlist_cases]


        output_list[[i]] <- data.frame("af" = af,
                                       "an" = an,
                                       "ID" = i,
                                       "cases" = dlist[[i]][,ind_dlist_cases])

      }

    }else{
      af <- 1-exp(-drop(as.matrix(Xpredall%*%coef[index])))
      an <- af*cases
      if(!is.null(group)){
        output_list_frame = data.frame(af = af, an = an, ID = group, cases = cases)
        output_list <- split(output_list_frame, output_list_frame$ID)

      }else{
        output_list = list(data.frame(af = af, an = an, cases = cases))
      }
    }

  # Total AF/AN
    if(tot) {
      an_tot = af_tot = NULL
      for (k in 1:length(output_list)){
        isna <- is.na(output_list[[k]]$an)
        af_tot <- sum(output_list[[k]]$an[!isna])/sum(output_list[[k]]$cases[!isna])
        an_tot <- af_tot*den[k]

        if(!is.null(group)){
          output_list[[k]] = data.frame(af = af_tot, an = an_tot, ID = output_list[[k]]$ID[1])
        }else{
          output_list[[k]] = data.frame(af = af_tot, an = an_tot)
        }
      }

    }


    # Credible intervals
    if(!tot && samples) {
      samples <- FALSE
      warning("simulation samples and exceedance probabilities only returned for tot=T")
    }

    if(samples){
    if(is.null(constraint_matrix)){
      coefsim  = simulate_from_precision(Prec, nsim = nsamples)
    }else{
      coefsim  = simulate_from_precision_con(Prec,constraint_matrix, nsim = nsamples)
    }
    coefsim <- coefsim + coef

    output_list_sim = list()

    if(kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
      for(i in names(output_list)){
        ind_spat_ID <- grepl(i, colnames(model$Xpred),
                             fixed = TRUE)

        afsim <- apply(coefsim,2, function(coefi) {

          ind_cb_spat_af <-  index & ind_spat_ID
          coefi <- c(coefi[ind_cb],coefi[ind_cb_spat_af])


          ind_dlist <-  !grepl("^(cases|group)", names(dlist[[i]]))
          ind_dlist_cases <-grepl("cases", names(dlist[[i]]))

          ani <- (1-exp(-drop(as.matrix(as.matrix(dlist[[i]][,ind_dlist])%*%coefi))))*dlist[[i]][,ind_dlist_cases]

          ani_result = sum(ani[!is.na(ani)])/sum((dlist[[i]][,ind_dlist_cases])[!is.na(ani)])
          return(ani_result)
        })

        ansim = NULL
        for (r in 1:nsamples){
          ansim[r] <- as.numeric(afsim[r])*den[i]
        }

        output_list_sim[[i]] = list(af = afsim, an = ansim, ID = i)
      }
    }else{
      afsim <- apply(coefsim,2, function(coefi) {
        ani <- as.numeric((1-exp(-drop(Xpredall%*%coefi[index])))*cases)
        ani_result = NULL
        if(!is.null(group)){
          for (r in 1:length(names(output_list))){
            ind_ID = which(names(output_list)[r] == unique(group))
            ani_result[r] = sum(ani[group==unique(group)[ind_ID] & !is.na(ani)])/sum(cases[group==unique(group)[ind_ID] & !is.na(ani)])
          }}
        else{
          ani_result[1] = sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
        }
        return(ani_result)
      })

      if(!is.null(group)){
        ansim = matrix(NA, nrow = nrow(afsim), ncol = ncol(afsim))
        for(i in names(output_list)){
          ind_ID = which(i == names(output_list))
          ansim[ind_ID,] <- as.numeric(afsim[ind_ID,])*den[i]

          output_list_sim[[i]] = list(af = afsim[ind_ID,], an = ansim[ind_ID,], ID = i)
        }

      }else{
        ansim <- as.numeric(afsim)*den

        output_list_sim = list(af = afsim, an = ansim)
      }
    }

    }


    if (!is.null(exc.threshold) & samples){
      if(!is.null(group)){
      for (r in 1:length(output_list)){
        output_list_sim[[r]]$exc.prob = sum(output_list_sim[[r]]$af>exc.threshold, na.rm=T)/length(output_list_sim[[r]]$af)
      }
      }else{
        output_list_sim$exc.prob = sum(output_list_sim$af>exc.threshold, na.rm=T)/length(output_list_sim$af)
      }


    }

    if(length(output_list) == 1){
      output_list = output_list[[1]]
    }
    if(samples){
      return(list(samples = output_list_sim, fit = output_list))
    }else{
      return(output_list)
    }


}
