# Helpers model comparison criteria   -----------------------------------------


# DIC calculation

calculate_DIC_tot <- function(xi_mode,
                              Prec,
                              Qv_mode,
                              Xv,
                              y,
                              offset,
                              v_mode = NULL,
                              model.family = "poisson",
                              constraint_matrix,
                              Ntrials = NULL){
  # DIC of mean

  if(model.family == "poisson"){
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+offset)
      return(sum((y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y))))
    }

    }else if (model.family == "binomial"){
      log_lik <- function(xi_mode,y, Xv, offset, v_mode){
        Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+offset)
        pxi_mean <- Cvxi_mean/(1+Cvxi_mean)
        return(sum(y * log(pxi_mean)+(Ntrials-y)*log(1-pxi_mean)))
      }

      }else if(model.family == "Gaussian"){
      log_lik <- function(xi_mode,y, Xv, offset, v_mode){
        tau <- exp(v_mode)[1]
        Cvxi_mean <-  as.numeric(Xv %*% xi_mode)+offset
        n <- length(y)
        return(-n*0.5*(log(2*pi)-log(tau))-0.5*tau*sum(((y-Cvxi_mean)^2)))
      }
    }else{
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+offset)
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))

      return( sum((1/exp(v_mode[1]))*((y * gammanb) - bgammanb) +
                    lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1)))
    }


  }

  Dev_of_mean <- -2*log_lik(xi_mode, y, Xv, offset, v_mode)

  # pd
  if(is.null(constraint_matrix)){
  pd_est <- dim(Prec)[1]-Hutchinson_pD(Prec, Qv_mode)
  }else{
    pd_est <- dim(Prec)[1]-Hutchinson_pD_con(Prec, Qv_mode, constraint_matrix)
  }

  DIC <- Dev_of_mean + 2*pd_est



  return(list(DIC = DIC, pd = pd_est))
}




# Combine WAIC and CPO in one function such that samples only have to be calculated once
# if both WAIC and CPO are requested.
calculate_CPO_WAIC <- function(xi_mode,
                               Prec,
                               Qv_mode,
                               Xv,
                               y,
                               offset,
                               v_mode = NULL,
                               model.family = "poisson",
                               CPO,
                               WAIC,
                               constraint_matrix,
                               Ntrials = NULL){
  # WAIC of mean

  if(model.family == "poisson"){
    loglik_samples_fun <- poisson_loglik_samples_cpp
  } else if (model.family == "binomial"){
    loglik_samples_fun <- logistic_loglik_samples_cpp
    }else if (model.family == "Gaussian"){
    loglik_samples_fun <- gaussian_loglik_samples_cpp
    }else{
    loglik_samples_fun <- nb_loglik_samples_cpp
  }


  if(is.null(constraint_matrix)){
  betas = simulate_from_precision(Prec, nsim = 500)
  }else{
    betas = simulate_from_precision_con(Prec,constraint_matrix, nsim = 500)
  }

  loglik_samples <- loglik_samples_fun(xi_mode, betas, Xv, y, offset, v_mode, Ntrials)

  if(WAIC){
    lppd <- sum(log(Matrix::rowMeans(exp(loglik_samples))))
      pWAIC2 <- sum(matrixStats::rowVars(loglik_samples))
    WAIC <- -2 * (lppd - pWAIC2)

    rm(lppd)

  }else{
    WAIC <- pWAIC2 <- NULL
  }

  if(CPO){
    m <- matrixStats::rowMaxs(-loglik_samples)
    harmonic_mean_inv <- Matrix::rowMeans(exp(-loglik_samples - m))
    logCPO <- -(log(harmonic_mean_inv) + m)
    rm(m) ; rm(harmonic_mean_inv)

    CPO_result <- exp(logCPO)
  }else{
    CPO_result <- NULL
  }


  rm(loglik_samples)

  return(list(CPO_result  = CPO_result, WAIC = WAIC, pd = pWAIC2))
}





