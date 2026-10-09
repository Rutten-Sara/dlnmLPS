attributable_risk <- function(basis, cases, model, cen, x, y = NULL, by = NULL, group = NULL,
                              type = "af", dir ="back", tot = TRUE, range_x = NULL,
                              range_y = NULL, sim = FALSE, nsim = 1000, exc.prob = NULL){

  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))

  if(missing(cen)){
    stop("Please specify a centering value")
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
    at_x <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2]),group=group) else
      matrix(rep(x,diff(lag)+1),length(x))
    at_y <- if(dir=="back") tsModel:::Lag(y,seq(lag[1],lag[2]), group=group) else
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
      at_x <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2]),group=group) else
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
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]), group=group)))
  }

  # Extract the model
    index <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                   gsub("\\s+", "", colnames(model$Xpred)),
                   fixed = TRUE)

    coef <- model$xi_mode
    Prec <- model$Prec[index,index]
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


    # Which type of basis?
    rhs <-  model$model[[3]]
    all_terms <- split_terms(rhs)

    index_term <- grepl(gsub("\\s+", "", deparse(substitute(basis))),
                        gsub("\\s+", "", all_terms),
                        fixed = TRUE)

    kind_term <- model$types[which(index_term)-1]

    dim_basis <- ncol(basis)
    ind_cb <- which(index)[1:dim_basis]
    xi_mode_cb = coef[ind_cb]

    # Effect-modified basis: interaction with effect-modification
    if (kind_term %in% c("modified DLNM", "joint modified DLNM")){
      if(is.null(by)){stop("Specify a by-variable")}
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
          z <-rep(0,length(attributes(type_smooth)$Dimnames[[2]]))
          z[index_by] =1
          z <- rbind(z,as.matrix(z, nrow = 1))
        }
      } else{
          z = as.matrix(by, ncol = 1)
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
    if(!tot && sim) {
      sim <- FALSE
      warning("simulation samples only returned for tot=T")
    }

    if(sim){
    if(is.null(constraint_matrix)){
      coefsim  = simulate_from_precision(Prec, nsim = 500)
    }else{
      coefsim  = simulate_from_precision_con(Prec,constraint_matrix, nsim = 500)
    }
    coefsim <- coefsim + coef[index]

    output_list_sim = list()

    if(kind_term %in% c("spatial DLNM", "joint spatial DLNM")){
      for(i in names(output_list)){
        ind_spat_ID <- grepl(i, colnames(model$Xpred),
                             fixed = TRUE)

        afsim <- apply(coefsim,2, function(coefi) {

          ind_cb_spat_af <-  which(index & ind_spat_ID) - max(ind_cb ) + dim_basis
          coefi <- coefi[c(1:dim_basis,ind_cb_spat_af)]


          ind_dlist <-  !grepl("^(cases|group)", names(dlist[[i]]))
          ind_dlist_cases <-grepl("cases", names(dlist[[i]]))

          ani <- (1-exp(-drop(as.matrix(as.matrix(dlist[[i]][,ind_dlist])%*%coefi))))*dlist[[i]][,ind_dlist_cases]

          ani_result = sum(ani[!is.na(ani)])/sum((dlist[[i]][,ind_dlist_cases])[!is.na(ani)])
          return(ani_result)
        })

        ansim = NULL
        for (r in 1:nsim){
          ansim[r] <- as.numeric(afsim[r])*den[i]
        }

        output_list_sim[[i]] = list(af = afsim, an = ansim, ID = i)
      }
    }else{
      afsim <- apply(coefsim,2, function(coefi) {
        ani <- as.numeric((1-exp(-drop(Xpredall%*%coefi)))*cases)
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


    if (!is.null(exc.prob) & sim){
      if(!is.null(group)){
      for (r in 1:length(output_list)){
        output_list_sim[[r]]$exc.prob = sum(output_list_sim[[r]]$af>exc.prob, na.rm=T)/length(output_list_sim[[r]]$af)
      }
      }else{
        output_list_sim$exc.prob = sum(output_list_sim$af>exc.prob, na.rm=T)/length(output_list_sim$af)
      }


    }

    if(sim){
      return(list(sim = output_list_sim, fit = output_list))
    }else{
      return(output_list)
    }


}
