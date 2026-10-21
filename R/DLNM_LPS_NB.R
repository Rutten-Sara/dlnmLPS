#' Fitting the model using LPS for negative binomial distribution
#' @noRd
DLNM_LPS_NB <- function(model,
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

  a.disp <- 1e-05
  b.disp <- 1e-05


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

  # For NB GLM with log-link
    Cv_xi <- function(xi, Xv){
      Cvxi <- exp(as.numeric(Xv %*% xi)+offset)
      return(Cvxi)
    }


    # Log conditional posterior of xi given v
    gamma.nb <- function( Cvxi,vphi) {
      res <- vphi * (log(Cvxi) - log(Cvxi + vphi))
      return(res)
    }
    bgamma.nb <- function( Cvxi, vphi) {
      res <- (-1)* (vphi^2) * (log(vphi)-log(vphi + Cvxi))
      return(res)
    }

    #log_pxi <- function(xi, Qv, Cvxi,v) {
    #   value <- (1/exp(v[1])) * sum((y * gamma.nb(xi, Cvxi,v)) - bgamma.nb(xi, Cvxi,v)) -
    #    .5 * Matrix::t(xi) %*% Qv %*% xi
    #  return(as.numeric(value))
    #}


    log_pxi <- function(xi, Qv, Cvxi,vphi) {
      invphi <- 1/vphi
      quad <- as.numeric(sum(xi * (Qv %*% xi)))
      value <- invphi * sum((y * gamma.nb(Cvxi,vphi)) - bgamma.nb(Cvxi,vphi)) - 0.5 * quad
      return(as.numeric(value))
    }


    # Gradient of log conditional posterior of xi given v
    var.nb <- function(Cvxi,invphi) {
      res <- Cvxi + invphi * (Cvxi ^ 2)
      return(res)
    }

    W.nb <- function(Cvxi, invphi, inv_varval) {
      res <- Matrix::Diagonal(x = ((Cvxi) ^ 2) * inv_varval)
      return(res)
    }


    Grad.logpxi <- function(xi,Qv, Cvxi,invphi, inv_varval){
      W.nbval <- W.nb( Cvxi, invphi, inv_varval)
      inv_Cxi <- 1/Cvxi
      D.nbval <- Matrix::Diagonal(x = inv_Cxi)
      value <- as.numeric(XtWcb%*%W.nbval%*%D.nbval%*%(y - Cvxi) -
                            Qv%*%xi)
      return(value)
    }

    # Hessian of log conditional posterior of xi given v
    M.nb <- function( Cvxi) Matrix::Diagonal(x = y - Cvxi)

    V.nb <- function( Cvxi, invphi, inv_varval){
      res <- Matrix::Diagonal(x = Cvxi * (inv_varval - (Cvxi*inv_varval^2) *
                                            (1 + 2*Cvxi*invphi)))
      return(res)
    }





    NR_xi <- function(xi0, Qv0, v0){

      epsilon <- 1e-03 # Stop criterion
      maxiter <- 300   # Maximum iterations
      iter <- 0        # Iteration counter
      #message("Newton Raphson")
      for (k in 1:maxiter) {
        Cvxi0 <- Cv_xi(xi0, Wcb)

        vphi <- exp(v0[1])
        inv0 <- 1/vphi
        inv_varval0 <- 1/var.nb(Cvxi0, inv0)

        grad_est =Grad.logpxi(xi0,Qv0, Cvxi0,inv0, inv_varval0)

        W.nb0 = methods::as(W.nb(Cvxi0,inv0, inv_varval0),"generalMatrix")
        MV.nb0 = methods::as(M.nb(Cvxi0)%*%V.nb(Cvxi0,inv0, inv_varval0),"generalMatrix")


        if(is.null(constraint_matrix)){
        dxi <- as.numeric((-1)* solve_sparse_cholesky_NB(Qv0, W.nb0, MV.nb0, Wcb, grad_est))
        }else{
          dxi <- as.numeric((-1)* solve_sparse_cholesky_NB_constr(Qv0, W.nb0, MV.nb0, Wcb, constraint_matrix, grad_est))
        }

        Xdxi <- as.numeric(Wcb %*% dxi)

        step <- 1
        iter.halving <- 1
        logpxi.current <- log_pxi(xi0, Qv0, Cvxi0, vphi)

        xi.new <- xi0 + dxi
        Cvxi.new <- Cv_xi(xi.new, Wcb)

        logpxi.new <- log_pxi(xi.new, Qv0, Cvxi.new, vphi)

        while (((is.na(logpxi.new) || logpxi.new <= logpxi.current) ||
                max(abs(step * dxi)) > 10) &&
               iter.halving <= 30) {
          step <- step * .5
          xi.new <- xi0 + (step * dxi)
          Cvxi.new <- Cv_xi(xi.new, Wcb)
          logpxi.new <- log_pxi(xi.new, Qv0, Cvxi.new, vphi)
          iter.halving <- iter.halving + 1
        }

        dist <- sqrt(sum((xi.new - xi0) ^ 2))
        iter <- iter + 1
        xi0 <- xi.new
        if(dist < epsilon) break
      }

      if(iter == maxiter){
        print("Warning: algorithm did not converge")
      }


      xistar <- xi0
      return(xistar)
    }

    # Initial values for log-penalty and log-overdispersion parameter

    dimxi <- dim(Wcb)[2]
    xi0_init = rep(0,dimxi)

    Qv_init <- Qv(v_init)

    xi_init <- NR_xi(xi0 = xi0_init, Qv0 = Qv_init, v0 = v_init)


    eta_temp <- as.numeric(Wcb %*% xi_init)+offset
    Cvxi_temp <- exp(eta_temp)


    # Log-posterior for the penalty- vector
    Mnb_v <- M.nb( Cvxi_temp)


    if(length(v_init)>0){
      log_pv <- function(v, xi){

        result <- tryCatch({
          phi = exp(v[1])
          inv_phi = 1/phi

          Qv_est <- Qv(v)


          varval <-  var.nb(Cvxi_temp, inv_phi)
          inv_varval <- 1/varval

          Vnb <- V.nb(Cvxi_temp, inv_phi, inv_varval)
          Wnb <- W.nb(Cvxi_temp, inv_phi, inv_varval)

          W_v_opt = Matrix::diag(-(Mnb_v%*%Vnb-Wnb))

          # I AM APPROXIMATING Wcb BY DIAGONAL TO SPEED UP CALCULATIONS
          XWX_est = XWX_func_approx(W_v_opt, methods::as(Matrix::Diagonal(x = Matrix::diag(Wcb)), "generalMatrix"))

          a1 <- tryCatch(
            -0.5 * determinant_Pv(XWX_est + Qv_est),
            error = function(e) {
              # fallback to exact
              XWX_exact <- XWX_func(W_v_opt, Wcb, XtWcb)
              -0.5 * determinant_Pv(XWX_exact + Qv_est)
            }
          )

          # Likelihood
          gammanb <- phi * (eta_temp-log(Cvxi_temp + phi))
          bgammanb <- bgamma.nb(Cvxi_temp, phi)

          a2 <- sum((1/phi)*((y * gammanb) - bgammanb) +
                      lgamma(y + phi) - lgamma(phi))


          a3 <- 0.5 * sum((xi * Qv_est) %*% xi)

          # Gamma prior overdispersion
          a4 <- a.disp*v[1] - b.disp*phi

          as.numeric(value <- a1+a2-a3+a4+logpv.hyper(v))
        }, error = function(e) {
          -1e10  # penalize the optimizer
        })

        return(result)
      }

      # Mode a posteriori estimate of v
      v_mode <- stats::optim(par = v_init,
                             fn = log_pv,
                             xi = xi_init,
                             method="Nelder-Mead",
                             control = list(fnscale = -1, reltol = 1e-8))$par
      rm(Qv_init)


      # Estimate of xi (regression parameter)
      Qv_mode <- Qv(v_mode)

      xi_mode <- NR_xi(xi0 = xi_init, Qv0 = Qv_mode, v0 = v_mode)

      if(refit>0){
        for (refit_iter in 1:refit){
          eta_temp <- as.numeric(Wcb %*% xi_mode)+offset
          Cvxi_temp <- exp(eta_temp)
          Mnb_v <- M.nb( Cvxi_temp)

          v_mode <- stats::optim(par = v_mode,
                                 fn = log_pv,
                                 xi = xi_mode,
                                 method="BFGS",
                                 control = list(fnscale = -1, reltol = 1e-8))$par


          Qv_mode <- Qv(v_mode)


          xi_mode <- NR_xi(xi0 = xi_mode, Qv0 = Qv_mode, v0 = v_mode)

        }



      }

    }else{
      xi_mode <- xi_init
      Qv_mode <- Qv_init
      v_mode <- v_init
      rm(Qv_init)
    }

    Cvxi_mode = Cv_xi(xi_mode, Wcb)
    vphi <- exp(v_mode[1])
    inv_phi <- 1/vphi
    inv_varval <- 1/var.nb(Cvxi_mode, inv_phi)

    W.nb_mode = methods::as(W.nb(Cvxi_mode,inv_phi, inv_varval),"generalMatrix")
    MV.nb_mode = methods::as(M.nb(Cvxi_mode)%*%V.nb(Cvxi_mode,inv_phi, inv_varval),"generalMatrix")


    Prec <- - Hess_logpxi_NB(Qv_mode,W.nb_mode, MV.nb_mode, Wcb)


    if(model.comparison$DIC == T){
      DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,Wcb, y, offset , v_mode,
                                 model.family = "NB", constraint_matrix)
      DIC_result = DIC_pd$DIC
      pd_result = DIC_pd$pd
      rm(DIC_pd)
    }else{
      DIC_result = NULL
      pd_result = NULL
    }


    if(model.comparison$WAIC == T | model.comparison$CPO == T ){
      WAIC_CPO = calculate_CPO_WAIC(xi_mode, Prec, Qv_mode,Wcb, y, offset, v_mode, model.family = "NB",
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
                   "family" = "NB")

    class(output) <- "dlnmLPS"
    return(output)

}
