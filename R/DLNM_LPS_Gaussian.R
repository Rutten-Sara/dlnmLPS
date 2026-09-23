#' Fitting the model using LPS for Gaussian distribution
#' @noRd
DLNM_LPS_Gaussian <- function(model,
                  data,
                  offset,
                  model.comparison,
                  refit,
                  base_env){


  types_term <- NULL

  # Split formula
  setup <- initialize_laplace_setup_constraint(model, data, TRUE)
  all_terms <- setup$all_terms

  rhs_fixed <- Reduce(
    function(x, y) call("+", x, y),
    if (setup$has_no_intercept) c(list(0), all_terms) else c(list(1), all_terms)
  )

  model <- stats::as.formula(
    call("~", model[[2]], rhs_fixed)
  )


  # constants
  constants <- list(
    zeta = 1e-05,
    a = 10^(-5),
    b = 10^(-5),
    nu = 3,
    a.rho = 0.5,
    b.rho = 0.5,
    data = data,
    env = base_env
  )


  rm(data)
  # Construct crossbasis
  y <- setup$y

  Wcb_init <- setup$Wcb
  if (is.null(Wcb_init) || length(Wcb_init) == 0) {
    Wcb_blocks <- list()
  } else if (ncol(Wcb_init) > 0) {
    Wcb_blocks <- list(Matrix::Matrix(Wcb_init, sparse = TRUE))
  } else {
    Wcb_blocks <- list()
  }

  rm(Wcb_init)

  setup$Wcb = Wcb_blocks

  all_fixed_terms = list()

  if (length(all_terms) == 0){
    if(setup$has_no_intercept){
      stop("Specify at least one term in your model")
    }
  }else{
    count_term = 1
    for(term in 1:length(all_terms)){
      setup$pen_dlnm <- c(setup$pen_dlnm,0)
      if(is_sm(all_terms[[term]])){
        count_term = count_term + 1
        sm_args_iter <- normalize_sm_call(all_terms[[term]])
        sm_args_iter <- check_sm_term_inputs(sm_args_iter, constants$data, constants$env)


        # Construct crossbasis
        if(sm_args_iter$arg$type=="rw1"){
          if(is.null(sm_args_iter$by)){
            setup <- build_component_rw(setup, sm_args_iter, term, all_terms, constants,1)
            types_term <- c(types_term, "RW1")
          }else if (sm_args_iter$arg.by$type %in% c("ind","ICAR","Leroux")){
            setup <- build_component_spatialrw(setup, sm_args_iter, term, all_terms, constants,1)
            types_term <- c(types_term, "spatial RW1")
          }
        }else if (sm_args_iter$arg$type=="rw2"){
          if(is.null(sm_args_iter$by)){
            setup <- build_component_rw(setup, sm_args_iter, term, all_terms, constants,2)
            types_term <- c(types_term, "RW2")
          }else if (sm_args_iter$arg.by$type %in% c("ind","ICAR","Leroux")){
            setup <- build_component_spatialrw(setup, sm_args_iter, term, all_terms, constants,2)
            types_term <- c(types_term, "spatial RW2")
          }
        }else if (sm_args_iter$arg$type %in% c("ind","ICAR","Leroux","Convolution")){
          if(is.null(sm_args_iter$by)){
            setup <- build_component_re(setup, sm_args_iter, term, all_terms, constants)
            types_term <- c(types_term, "re")
          }else{
            setup <-build_component_rwspatial(setup, sm_args_iter, term, all_terms, constants)
            types_term <- c(types_term, "temporal re")
          }
        } else if (sm_args_iter$arg$type %in% c("ps","ns","bs")){
          setup <- build_component_spline(setup, sm_args_iter, term, all_terms,constants)
          types_term <- c(types_term, "spline")
        } else if ("joint.crossbasis" %in% class(sm_args_iter$x)){
          if(is.null(sm_args_iter$arg.by$type)){
            setup <- build_component_2Ddlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "joint DLNM")
          }else if (sm_args_iter$arg.by$type %in% c("ind","ICAR","Leroux")){
            setup<-build_component_spatial2Ddlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "joint spatial DLNM")
          } else{
            setup<-build_component_modified2Ddlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "joint modified DLNM")
          }
        }else if(sm_args_iter$arg$type == "crossbasis" ){
          if(is.null(sm_args_iter$arg.by$type)){
            setup<-build_component_dlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "DLNM")
          } else if (sm_args_iter$arg.by$type %in% c("ind","ICAR","Leroux")){
            setup<-build_component_spatialdlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "spatial DLNM")
          } else{
            setup<-build_component_modifieddlnm(setup, sm_args_iter, term, all_terms,constants)
            types_term <- c(types_term, "modified DLNM")
          }
        }
      }else{
        all_fixed_terms <- c(all_fixed_terms, all_terms[[term]])
        types_term <- c(types_term, "fixed")
      }

    }
    if(length(all_fixed_terms) != 0){
      setup <- build_component_fixed(model,setup, all_fixed_terms, count_term, constants)

    }
  }


  # Bind design matrix
  Wcb_blocks <- setup$Wcb
  if (length(Wcb_blocks) == 0) {
    Wcb <- Matrix::Matrix(matrix(0, nrow = length(y), ncol = 0), sparse = TRUE)
  } else if (length(Wcb_blocks) == 1) {
    Wcb <- Wcb_blocks[[1]]
  } else {
    Wcb <- Reduce(Matrix::cbind2, Wcb_blocks)
  }

  rm(Wcb_blocks)

  pen_dlnm <- setup$pen_dlnm
  Pv.list <- setup$Pv.list
  logpv.fixed <- setup$logpv.fixed
  v_init <- setup$v_init
  sm_list <- setup$sm_list

  constraint_matrix <- setup$constraint
  has_intercept <- !setup$has_no_intercept
  rm(setup)


  #Bind intercept
  if(has_intercept){
    intercept_block <- Matrix::Matrix(matrix(1, nrow = length(y), ncol = 1), sparse = TRUE)
    colnames(intercept_block) <- "(Intercept)"
    Wcb <- Matrix::cbind2(intercept_block, Wcb)
    constraint_matrix <- c(list(Matrix::Matrix(0,1,1)), constraint_matrix)
    all_terms <- c(list(1), all_terms)

    logpv.fixed <- c(list(function(v) 0), logpv.fixed)
    Pv.list <- c(list(function(v){
      methods::as(Matrix::Diagonal(n = 1,
                                   x = constants$zeta), "generalMatrix")

    }), Pv.list)
  }



  constraint_matrix <- methods::as(do.call(Matrix::bdiag, constraint_matrix), "CsparseMatrix")

  keep <- Matrix::rowSums(constraint_matrix != 0) > 0
  constraint_matrix <- constraint_matrix[keep, , drop = FALSE]

  if(nrow(constraint_matrix) == 0){
    constraint_matrix <- NULL
  }





  # If no population offset --> put equal to 1 for every area
  if (is.null(offset)){
    offset <- rep(0, length(y))
  }



  ind <- (Matrix::rowSums(is.na(Wcb)) + is.na(y) + is.na(offset)) == 0
  Wcb <- Wcb[ind, , drop = FALSE]
  y <- y[ind]
  offset <- offset[ind]
  XtWcb <- Matrix::t(Wcb)
  rm(ind)

  #Precision matrix for parameter xi
  Qv <- function(v){
    blocks <- lapply(seq_along(Pv.list), function(terms_q) {
      Pv.list[[terms_q]](v)
    })
    if (length(blocks) == 1L) {
      return(methods::as(blocks[[1]], "CsparseMatrix"))
    }
    methods::as(do.call(Matrix::bdiag, blocks), "CsparseMatrix")
  }


  # Hyperparameters
    logpv.hyper <- function(v) {
      sum(sapply(seq_along(Pv.list),
                 function(i) logpv.fixed[[i]](v)))
    }

    # Posterior mode (Gaussian)
    XtX <- XX_func(Wcb)
    Xty <- as.numeric(XtWcb%*%(y-offset))
    Qv_init <- Qv(v_init)
    dimxi <- dim(Wcb)[2]

    if(is.null(constraint_matrix)){
      xi_init <- solve_gaussian_mode(
        Qv_init,
        XtX,
        Xty,
        exp(v_init[1])
      )
    }else{
      xi_init <- solve_gaussian_mode_constr(
        Qv_init,
        XtX,
        Xty,
        exp(v_init[1]),
        constraint_matrix
      )
    }

    error_y_init <- y-offset-as.numeric(Wcb%*%xi_init)
    n_y <- length(y)

    if(length(v_init)>0){

      log_pv <- function(v, xi, error_y){

        result <- tryCatch({
          tau <- exp(v[1])
          Qv_est <- Qv(v)
          a1 <- -0.5*determinant_Pv(tau*XtX+Qv_est)

          a4 <- 0.5 * sum((xi * Qv_est) %*% xi)

          # Gaussian likelihood
          a2 <- n_y*0.5*v[1]-0.5*tau*sum((error_y)^2)

          as.numeric(value <- a1+a2-a4+logpv.hyper(v)) #+a2-a3
        }, error = function(e) {
          -1e10  # penalize the optimizer
        })

        return(result)
      }

      # Mode a posteriori estimate of v
      v_mode <- stats::optim(par = v_init,
                             fn = log_pv,
                             xi = xi_init,
                             error_y = error_y_init,
                             method="Nelder-Mead",
                             control = list(fnscale = -1, reltol = 1e-8))$par

      # Estimate of xi (regression parameter)
      rm(Qv_init)

      Qv_mode <- Qv(v_mode)

      if(is.null(constraint_matrix)){
        xi_mode <- solve_gaussian_mode(
          Qv_mode,
          XtX,
          Xty,
          exp(v_mode[1])
        )
      }else{
        xi_mode <- solve_gaussian_mode_constr(
          Qv_mode,
          XtX,
          Xty,
          exp(v_mode[1]),
          constraint_matrix
        )
      }


      if(refit>0){
        for (refit_iter in 1:refit){
          error_y_mode <- y-offset-as.numeric(Wcb%*%xi_mode)


          v_mode <- stats::optim(par = v_mode,
                                 fn = log_pv,
                                 xi = xi_mode,
                                 error_y = error_y_mode,
                                 method="BFGS",
                                 control = list(fnscale = -1, reltol = 1e-8))$par


          Qv_mode <- Qv(v_mode)

          if(is.null(constraint_matrix)){
            xi_mode <- solve_gaussian_mode(
              Qv_mode,
              XtX,
              Xty,
              exp(v_mode[1])
            )
          }else{
            xi_mode <- solve_gaussian_mode_constr(
              Qv_mode,
              XtX,
              Xty,
              exp(v_mode[1]),
              constraint_matrix
            )
          }


        }



      }

    }else{
      xi_mode = xi_init
      v_mode = v_init
      Qv_mode <- Qv_init
      rm(Qv_init)
    }

    Prec <- exp(v_mode[1])*XtX+Qv_mode


    if(model.comparison$DIC == T){
      DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,Wcb, y, offset , v_mode, model.family = "Gaussian",
                                 constraint_matrix)
      DIC_result = DIC_pd$DIC
      pd_result = DIC_pd$pd
      rm(DIC_pd)
    }else{
      DIC_result = NULL
      pd_result = NULL
    }


    if(model.comparison$WAIC == T | model.comparison$CPO == T ){
      WAIC_CPO = calculate_CPO_WAIC(xi_mode, Prec, Qv_mode,Wcb, y, offset, v_mode, model.family = "Gaussian",
                                   model.comparison$CPO, model.comparison$WAIC, constraint_matrix)
      WAIC_result = WAIC_CPO$WAIC
      pd_WAIC_result = WAIC_CPO$pd
      CPO_result = WAIC_CPO$CPO_result
      rm(WAIC_CPO)
    }else{
      WAIC_result = NULL
      pd_WAIC_result = NULL
      CPO_result = NULL
    }



    names(xi_mode) = colnames(Wcb)
    output <- list("model" = model,
                   "xi_mode" = xi_mode,
                   "Prec" = Prec,
                   "Xpred" = Wcb,
                   "y" = y,
                   "offset"= offset,
                   "smooth" = sm_list,
                   "v_mode" = v_mode,
                   "DIC" = DIC_result,
                   "pd DIC" = pd_result,
                   "WAIC" = WAIC_result,
                   "pd WAIC" = pd_WAIC_result,
                   "CPO" = CPO_result,
                   "constraints" = constraint_matrix,
                   "types" = types_term,
                   "family" = "Gaussian")

    class(output) <- "DLNMLPS"
    return(output)

}
