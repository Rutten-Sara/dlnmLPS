# Building crossbasis DLNM functions ----------------------------

#' Construction of DLNM components
#' @noRd
build_component_dlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu
  data <- constants$data
  env <- constants$env

  Wcb_part <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(Wcb_part) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Wcb_part), sep = ""
  )

  Wcb <- add_wcb_block(Wcb, Wcb_part)
  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part))

  penorder <- sm_args_iter$arg$pen

  pen_start <- cumsum(pen_dlnm)[term]


  if(!is.null(penorder)){
    vx <- attributes(sm_args_iter$x)$df[1]
    vl <- attributes(sm_args_iter$x)$df[2]


    # Penalty for exposure
    Dx <- Matrix::Diagonal(vx+1, x = NULL)
    for (k in 1:penorder) Dx <- Matrix::diff(Dx)
    Px <- Matrix::t(Dx) %*% Dx
    Px <- Px[-1,-1]
    Px <- Px + Matrix::Diagonal(vx, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    Pl_add <- Dl_add + diag(1e-12, vl)



    Pv.list[[term]] <- function(v) exp(v[pen_start + 1])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
      exp(v[pen_start + 2])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
      exp(v[pen_start + 3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)

    Pv_term <- Pv.list[[term]]
    logpv.fixed[[term]] <- function(v) 0.5 * nu * (v[pen_start + 1]+v[pen_start + 2]+v[pen_start + 3])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1]))+log(b + 0.5*nu*exp(v[pen_start + 2]))+
                      log(b + 0.5*nu*exp(v[pen_start + 3])))+
      0.5*determinant_Pv(Pv_term(v))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+3

    v_init <- c(v_init, rep(5, 3))
    names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                       paste("Penalty DLNM",deparse1(all_terms[[term]]), 1:3))

  }else{
    logpv.fixed[[term]] <- function(v) 0
    nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }
  }


  state$v_init<-v_init
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$Wcb<-Wcb
  state$constraint <- constraint_matrix

  return(state)
}




#' Construction of spatially varying DLNM components
#' @noRd
build_component_spatialdlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu
  a.rho <- constants$a.rho
  b.rho <- constants$b.rho
  data <- constants$data
  env <- constants$env

  # main effect
  Wcb_part <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(Wcb_part) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Wcb_part), sep = ""
  )

  Wcb <- add_wcb_block(Wcb, Wcb_part)
  penorder <- sm_args_iter$arg$pen

  #spatial effect
  crossbasis_spat = Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(crossbasis_spat) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(crossbasis_spat), sep = ""
  )


  Z.rand_spat <- Matrix::sparse.model.matrix(~ sm_args_iter$by + 0)
  colnames(Z.rand_spat) <- sub(
    "^sm_args_iter\\$by",
    deparse(normalize_sm_call(all_terms[[term]])$by),
    colnames(Z.rand_spat)
  )

  q.rand <- dim(Z.rand_spat)[2]

  Mat_crossbasis_spat <- Matrix::Matrix(crossbasis_spat) # cross basis
  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(Mat_crossbasis_spat), Matrix::t(Z.rand_spat)))

  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(crossbasis_spat), colnames(Z.rand_spat), paste, sep = ":")
  ))


  # Global design matrix
  Wcb <- add_wcb_block(Wcb, interaction_matrix)

  rm(crossbasis_spat, Mat_crossbasis_spat, Z.rand_spat, interaction_matrix)


  pen_start <- cumsum(pen_dlnm)[term]

  if(!is.null(penorder)){
    vx <- attributes(sm_args_iter$x)$df[1]
    vl <- attributes(sm_args_iter$x)$df[2]

    # Penalty for exposure
    Dx <- Matrix::Diagonal(vx+1, x = NULL)
    for (k in 1:penorder) Dx <- Matrix::diff(Dx)
    Px <- Matrix::t(Dx) %*% Dx
    Px <- Px[-1,-1]
    Px <- Px + Matrix::Diagonal(vx, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    Pl_add <- Dl_add + diag(1e-12, vl)


    Pv1_term <- function(v) exp(v[pen_start + 1])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
      exp(v[pen_start + 2])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
      exp(v[pen_start + 3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)

    logpv1 <- function(v) 0.5 * nu * (v[pen_start + 1]+v[pen_start + 2]+v[pen_start + 3])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1]))+log(b + 0.5*nu*exp(v[pen_start + 2]))+
                      log(b + 0.5*nu*exp(v[pen_start + 3])))+
      0.5*determinant_Pv(Pv1_term(v))

    pen_dlnm[term+1] = pen_dlnm[term+1]+3

    v_init <- c(v_init, rep(1, 3))
    names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                       paste("Penalty common DLNM",deparse1(all_terms[[term]]), 1:3))


    # Spatial interaction
    if (sm_args_iter$arg.by$type == "ind"){
      Sv_list <- Matrix::Diagonal(n = q.rand, x = 1)

      Pv2_term <- function(v) exp(v[pen_start + 4])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Sv_list) +
        exp(v[pen_start + 5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Sv_list) +
        exp(v[pen_start + 6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Sv_list)


      logpv2<- function(v) 0.5 * nu * (v[pen_start + 4]+v[pen_start + 5]+v[pen_start + 6])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 4]))+log(b + 0.5*nu*exp(v[pen_start + 5]))+
                        log(b + 0.5*nu*exp(v[pen_start + 6])))+
        0.5*determinant_Pv(Pv2_term(v))

      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      v_init <- c(v_init, rep(3, 3))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                         paste("Penalty spatial DLNM",deparse1(all_terms[[term]]), 1:3))

      pen_dlnm[term+1] = pen_dlnm[term+1]+3

      constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))

    } else if (sm_args_iter$arg.by$type == "ICAR"){
      neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
      Rn <- make_graph_precision(neig.map, q.rand)

      Pv2_term <- function(v) exp(v[pen_start+4])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Rn) +
        exp(v[pen_start+5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Rn) +
        exp(v[pen_start+6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Rn)

      logpv2 <- function(v) 0.5 * nu * (v[pen_start+4]+v[pen_start+5]+v[pen_start+6])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+4]))+log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+6])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1]+3
      v_init <- c(v_init, rep(3,3))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                         paste("Penalty spatial DLNM",deparse1(all_terms[[term]]), 1:3))

      ncol_Wcb <- ncol(Wcb_part)
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:ncol_Wcb, each = q.rand),
                                                         j =  (ncol_Wcb+1):(ncol_Wcb*(1+q.rand)),
                                                         x = 1,
                                                         dims = c(ncol_Wcb, ncol_Wcb*(1+q.rand)))


    }  else if (sm_args_iter$arg.by$type == "Leroux"){
      neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
      Rn <- make_graph_precision(neig.map, q.rand)

      Sv_list <- function(v) {
        make_leroux_precision(0, v[pen_start + 7], q.rand, Rn)
      }

      Pv2_term <- function( v) exp(v[pen_start+4])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Sv_list(v)) +
        exp(v[pen_start+5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Sv_list(v)) +
        exp(v[pen_start+6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Sv_list(v))


      logpv2<- function(v)  {
        value <- 0.5 * nu * (v[pen_start+4]+v[pen_start+5]+v[pen_start+6])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+4]))+log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+6])))+
          0.5*determinant_Pv(Pv2_term(v)) +
          a.rho*v[pen_start+7] - (a.rho + b.rho)*log(1 + exp(v[pen_start+7]))
        return(as.numeric(value))}


      logpv.fixed[[term]] <- function(v) logpv1( v) + logpv2(v)
      Pv.list[[term]] <- function( v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      v_init <- c(v_init, c(rep(3,3),1))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                         paste("Penalty DLNM",deparse1(all_terms[[term]]), 1:3),
                         paste("Rho DLNM",deparse1(all_terms[[term]])))

      pen_dlnm[term+1] <- pen_dlnm[term+1]+4

      constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))
    }}else{
      logpv1 <- function(v) 0
      nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
      Pv1_term <- function(v){
        methods::as(Matrix::Diagonal(n = nbasis,
                            x = zeta), "generalMatrix")
      }

      if (sm_args_iter$arg.by$type == "ind"){
        Sv_list <- Matrix::Diagonal(n = q.rand, x = 1)

        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Sv_list
        n_d <- ncol(sm_args_iter$x)*q.rand


        logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))

        v_init <- c(v_init, 5)
        names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])))

        logpv.fixed[[term]] <- function( v) logpv1( v) + logpv2(v)
        Pv.list[[term]] <- function( v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+1
        constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))

      } else if (sm_args_iter$arg.by$type == "ICAR"){
        neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
        Rn <- make_graph_precision(neig.map, q.rand)

        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Rn
        n_d <- ncol(sm_args_iter$x)*q.rand


        logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))

        v_init <- c(v_init, 5)
        names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])))

        logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
        Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+1

        ncol_Wcb <- ncol(Wcb_part)
        constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:ncol_Wcb, each = q.rand),
                                                           j =  (ncol_Wcb+1):(ncol_Wcb*(1+q.rand)),
                                                           x = 1,
                                                           dims = c(ncol_Wcb, ncol_Wcb*(1+q.rand)))

      }else if (sm_args_iter$arg.by$type == "Leroux"){
        neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
        Rn <- make_graph_precision(neig.map, q.rand)

        Sv_list <- function(v) {
          make_leroux_precision(0, v[pen_start + 2], q.rand, Rn)
        }

        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Sv_list(v)


        logpv2<- function(v)  {
          value <- 0.5 * nu * (v[pen_start+1])-
            (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))+
            0.5*determinant_Pv(Pv2_term(v)) +
            a.rho*v[pen_start+2] - (a.rho + b.rho)*log(1 + exp(v[pen_start+2]))
          return(as.numeric(value))}


        logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
        Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        v_init <- c(v_init, c(5,1))
        names(v_init) <- c(names(v_init)[1:(length(v_init)-2)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])),
                           paste("Rho DLNM",deparse1(all_terms[[term]])))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+2

        constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))
      }
    }

  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)
}



#' Construction of effect modified DLNM components
#' @noRd
build_component_modifieddlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init
  sm_list <- state$sm_list

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu

  data <- constants$data
  env <- constants$env

  # main effect
  crossbasis_em <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(crossbasis_em) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(crossbasis_em), sep = ""
  )

  Wcb <- add_wcb_block(Wcb, crossbasis_em)

  penorder <- sm_args_iter$arg$pen

  pen_start <- cumsum(pen_dlnm)[term]

  #modifying effect
  if(sm_args_iter$arg.by$type == "linear"){
    sm_list[[term]] = matrix(sm_args_iter$by, ncol = 1)
    vz <- 1

    colnames(sm_list[[term]]) <- deparse(normalize_sm_call(all_terms[[term]])$by)


  }else if (sm_args_iter$arg.by$type == "factor"){
    sm_list[[term]] = Matrix::sparse.model.matrix(~ as.factor(sm_args_iter$by))[,-1, drop = F]
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- sub(
      "^sm_args_iter\\$by",
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]])
    )

  }else if(sm_args_iter$arg.by$type == "ps"){
    sm_list[[term]] <- dlnm::ps(sm_args_iter$by,
                                diff = sm_args_iter$arg.by$pen,
                                df = sm_args_iter$arg.by$df,
                                intercept = F)

    D_z = Matrix::Matrix(attr(sm_list[[term]],"S"))+
      Matrix::Diagonal(n = ncol(sm_list[[term]]), x = 1e-12)

    vz <- ncol(sm_list[[term]])

    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }else if (sm_args_iter$arg.by$type == "ns"){
    sm_list[[term]] <- splines::ns(sm_args_iter$by,
                                   df = sm_args_iter$arg.by$df,
                                   knots = sm_args_iter$arg.by$knots,
                                   intercept =F,
                                   Boundary.knots = sm_args_iter$arg.by$Boundary.knots)
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }else{
    sm_list[[term]] <- splines::bs(sm_args_iter$by,
                                   df = sm_args_iter$arg.by$df,
                                   knots = sm_args_iter$arg.by$knots,
                                   degree = sm_args_iter$arg.by$degree,
                                   intercept = F,
                                   Boundary.knots = sm_args_iter$arg.by$Boundary.knots)
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }


  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(crossbasis_em), Matrix::Matrix(Matrix::t(sm_list[[term]]))))

  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(crossbasis_em), colnames(sm_list[[term]]), paste, sep = ":")
  ))

  # Global design matrix
  Wcb <- add_wcb_block(Wcb, interaction_matrix)
  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(interaction_matrix)+ncol(crossbasis_em))

  rm(interaction_matrix, crossbasis_em)



  if(!is.null(penorder)){

    vx <- attributes(sm_args_iter$x)$df[1]
    vl <- attributes(sm_args_iter$x)$df[2]

    # Penalty for exposure
    Dx <- Matrix::Diagonal(vx+1, x = NULL)
    for (k in 1:penorder) Dx <- Matrix::diff(Dx)
    Px <- Matrix::t(Dx) %*% Dx
    Px <- Px[-1,-1]
    Px <- Px + Matrix::Diagonal(vx, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    #Dl_add <- diag(rep(0:1,c(6,4)))
    Pl_add <- Dl_add + diag(1e-12, vl)


    Pv1_term <- function(v) exp(v[pen_start + 1])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
      exp(v[pen_start + 2])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
      exp(v[pen_start + 3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)

    logpv1 <- function(v) 0.5 * nu * (v[pen_start + 1]+v[pen_start + 2]+v[pen_start + 3])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1]))+log(b + 0.5*nu*exp(v[pen_start + 2]))+
                      log(b + 0.5*nu*exp(v[pen_start + 3])))+
      0.5*determinant_Pv(Pv1_term(v))

    pen_dlnm[term+1] = pen_dlnm[term+1]+3

    v_init <- c(v_init, rep(2, 3))

    names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                       paste("Penalty common DLNM",deparse1(all_terms[[term]]), 1:3))

    if(sm_args_iter$arg.by$type == "ps"){
      Pv2_term <- function(v) exp(v[pen_start+4])*(Px%x%Matrix::Diagonal(n = vl*vz, x = NULL)) +
        exp(v[pen_start+5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = vz, x = NULL)) +
        exp(v[pen_start+6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = vz, x = NULL))+
        exp(v[pen_start+7])*(Matrix::Diagonal(n = vx*vl, x = NULL)%x%D_z)

      logpv2 <- function(v) 0.5 * nu * (v[pen_start+5]+v[pen_start+6]+v[pen_start+7]+v[pen_start+4])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+5]))+
                        log(b + 0.5*nu*exp(v[pen_start+6]))+log(b + 0.5*nu*exp(v[pen_start+7]))+log(b + 0.5*nu*exp(v[pen_start+4])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 4
      v_init <- c(v_init, rep(2, 4))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]]), 1:4))


    }else{
      Pv2_term <- function(v) exp(v[pen_start+4])*(Px%x%Matrix::Diagonal(n = vl*vz, x = NULL)) +
        exp(v[pen_start+5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = vz, x = NULL)) +
        exp(v[pen_start+6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = vz, x = NULL))

      logpv2 <- function(v) 0.5 * nu * (v[pen_start+5]+v[pen_start+6]+v[pen_start+4])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+5]))+
                        log(b + 0.5*nu*exp(v[pen_start+6]))+log(b + 0.5*nu*exp(v[pen_start+4])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 3
      v_init <- c(v_init, rep(2, 3))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-3)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]]), 1:3))
    }




  }else{
    logpv1 <- function(v) 0
    nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
    Pv1_term <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }

    if(sm_args_iter$arg.by$type == "ps"){
      Pv2_term <- function(v) exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)%x%D_z)
      n_d <- ncol(sm_args_iter$x)*ncol(D_z)


      logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 1
      v_init <- c(v_init, 2)
      names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]])))

    }else{

      logpv2 <- function(v) 0
      Pv2_term <- function(v){
        methods::as(Matrix::Diagonal(n = ncol(sm_args_iter$x)*vz,
                            x = zeta), "generalMatrix")
      }

      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))


    }

  }

  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$sm_list<-sm_list
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)
}


#' Construction of joint DLNM components
#' @noRd
build_component_2Ddlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu
  data <- constants$data
  env <- constants$env

  Wcb_part <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(Wcb_part) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Wcb_part), sep = ""
  )

  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part))
  Wcb <- add_wcb_block(Wcb, Wcb_part)
  penorder <- sm_args_iter$arg$pen

  pen_start <-  cumsum(pen_dlnm)[term]

  if(!is.null(penorder)){
    vx1 <- attributes(sm_args_iter$x)$argvar1$df
    vx2 <- attributes(sm_args_iter$x)$argvar2$df
    vl <- attributes(sm_args_iter$x)$arglag$df


    # Penalty for exposure 1
    Dx1 <- Matrix::Diagonal(vx1, x = NULL)
    for (k in 1:penorder) Dx1 <- Matrix::diff(Dx1)
    Px1 <- Matrix::t(Dx1) %*% Dx1
    Px1 <- Px1 + Matrix::Diagonal(vx1, x = 1e-12)


    # Penalty for exposure 2
    Dx2 <- Matrix::Diagonal(vx2, x = NULL)
    for (k in 1:penorder) Dx2 <- Matrix::diff(Dx2)
    Px2 <- Matrix::t(Dx2) %*% Dx2
    Px2 <- Px2 + Matrix::Diagonal(vx2, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    Pl_add <- Dl_add + diag(1e-12, vl)


    if(attributes(sm_args_iter$x)$is_identifiable){

      Pv.list[[term]] <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)] +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)]  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)
    }else{

      Pv.list[[term]] <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)) +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)

    }

    Pv_term <- Pv.list[[term]]

    logpv.fixed[[term]] <- function(v) 0.5 * nu * (v[pen_start+1]+v[pen_start+2]+
                                                     v[pen_start+3]+v[pen_start+4])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1]))+log(b + 0.5*nu*exp(v[pen_start+2]))+
                      log(b + 0.5*nu*exp(v[pen_start+3]))+log(b + 0.5*nu*exp(v[pen_start+4])))+
      0.5*determinant_Pv(Pv_term(v))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+4

    v_init <- c(v_init, rep(5, 4))
    names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                       paste("Penalty DLNM",deparse1(all_terms[[term]]), 1:4))

  }else{
    logpv.fixed[[term]] <- function(v) 0
    nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }
  }

  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)

}




#' Construction of spatially varying joint DLNM components
#' @noRd
build_component_spatial2Ddlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu
  a.rho <- constants$a.rho
  b.rho <- constants$b.rho
  data <- constants$data
  env <- constants$env

  # main effect
  Wcb_part <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(Wcb_part) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Wcb_part), sep = ""
  )

  Wcb <- add_wcb_block(Wcb, Wcb_part)
  penorder <- sm_args_iter$arg$pen

  #spatial effect
  crossbasis_spat = Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(crossbasis_spat) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(crossbasis_spat), sep = ""
  )


  Z.rand_spat <- Matrix::sparse.model.matrix(~ sm_args_iter$by + 0)
  colnames(Z.rand_spat) <- sub(
    "^sm_args_iter\\$by",
    deparse(normalize_sm_call(all_terms[[term]])$by),
    colnames(Z.rand_spat)
  )

  q.rand <- dim(Z.rand_spat)[2]

  Mat_crossbasis_spat <- Matrix::Matrix(crossbasis_spat) # cross basis
  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(Mat_crossbasis_spat), Matrix::t(Z.rand_spat)))

  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(crossbasis_spat), colnames(Z.rand_spat), paste, sep = ":")
  ))


  # Global design matrix
  Wcb <- add_wcb_block(Wcb, interaction_matrix)
  rm(interaction_matrix, Z.rand_spat, crossbasis_spat, Mat_crossbasis_spat)

  pen_start <- cumsum(pen_dlnm)[term]

  if(!is.null(penorder)){

    vx1 <- attributes(sm_args_iter$x)$argvar1$df
    vx2 <- attributes(sm_args_iter$x)$argvar2$df
    vl <- attributes(sm_args_iter$x)$arglag$df

    # Penalty for exposure 1
    Dx1 <- Matrix::Diagonal(vx1, x = NULL)
    for (k in 1:penorder) Dx1 <- Matrix::diff(Dx1)
    Px1 <- Matrix::t(Dx1) %*% Dx1
    Px1 <- Px1 + Matrix::Diagonal(vx1, x = 1e-12)


    # Penalty for exposure 2
    Dx2 <- Matrix::Diagonal(vx2, x = NULL)
    for (k in 1:penorder) Dx2 <- Matrix::diff(Dx2)
    Px2 <- Matrix::t(Dx2) %*% Dx2
    Px2 <- Px2 + Matrix::Diagonal(vx2, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    Pl_add <- Dl_add + diag(1e-12, vl)



    if(attributes(sm_args_iter$x)$is_identifiable){


      Pv1_term <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)] +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)]  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)
    }else{

      Pv1_term <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)) +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)

    }


    logpv1 <- function(v) 0.5 * nu * (v[pen_start+1]+v[pen_start+2]+
                                        v[pen_start+3]+v[pen_start+4])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1]))+log(b + 0.5*nu*exp(v[pen_start+2]))+
                      log(b + 0.5*nu*exp(v[pen_start+3]))+log(b + 0.5*nu*exp(v[pen_start+4])))+
      0.5*determinant_Pv(Pv1_term(v))

    pen_dlnm[term+1] = pen_dlnm[term+1]+4

    v_init <- c(v_init, rep(1, 4))
    names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                       paste("Penalty DLNM",deparse1(all_terms[[term]]), 1:4))

    # Spatial interaction
    if (sm_args_iter$arg.by$type == "ind"){
      Sv_list <- Matrix::Diagonal(n = q.rand, x = 1)

      if(attributes(sm_args_iter$x)$is_identifiable){

        Pv2_term <- function(v) exp(v[pen_start+5])*((Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Sv_list +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Sv_list +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl)%x%Sv_list +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)%x%Sv_list
      }else{

        Pv2_term <- function(v) exp(v[pen_start+5])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)%x%Sv_list) +
          exp(v[pen_start+6])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL)%x%Sv_list)  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl)%x%Sv_list +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)%x%Sv_list

      }


      logpv2<- function(v) 0.5 * nu * (v[pen_start + 7]+v[pen_start + 8]+v[pen_start + 5]+v[pen_start + 6])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 7]))+log(b + 0.5*nu*exp(v[pen_start + 5]))+
                        log(b + 0.5*nu*exp(v[pen_start + 6]))+log(b + 0.5*nu*exp(v[pen_start + 8])))+
        0.5*determinant_Pv(Pv2_term(v))

      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      v_init <- c(v_init, rep(3, 4))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                         paste("Penalty spatial DLNM",deparse1(all_terms[[term]]), 1:4))

      pen_dlnm[term+1] = pen_dlnm[term+1]+4

      constraint_matrix[[term]] <-Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))


    } else if (sm_args_iter$arg.by$type == "ICAR"){
      neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
      Rn <- make_graph_precision(neig.map, q.rand)

      if(attributes(sm_args_iter$x)$is_identifiable){


        Pv2_term <- function(v) exp(v[pen_start+5])*((Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Rn +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Rn  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl)%x%Rn +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)%x%Rn
      }else{

        Pv2_term <- function(v) exp(v[pen_start+5])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)%x%Rn) +
          exp(v[pen_start+6])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL)%x%Rn)  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl)%x%Rn +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)%x%Rn

      }

      logpv2 <- function(v) 0.5 * nu * (v[pen_start+7]+v[pen_start+5]+v[pen_start+6]+v[pen_start+8])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+7]))+log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+6]))+log(b + 0.5*nu*exp(v[pen_start+8])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1]+4
      v_init <- c(v_init, rep(3,4))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                         paste("Penalty spatial DLNM",deparse1(all_terms[[term]]), 1:4))

      ncol_Wcb <- ncol(Wcb_part)
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:ncol_Wcb, each = q.rand),
                                                         j =  (ncol_Wcb+1):(ncol_Wcb*(1+q.rand)),
                                                         x = 1,
                                                         dims = c(ncol_Wcb, ncol_Wcb*(1+q.rand)))



    }  else if (sm_args_iter$arg.by$type == "Leroux"){
      neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
      Rn <- make_graph_precision(neig.map, q.rand)

      Sv_list <- function(v) {
        make_leroux_precision(0, v[pen_start + 9], q.rand, Rn)
      }


      if(attributes(sm_args_iter$x)$is_identifiable){


        Pv2_term <- function(v) exp(v[pen_start+5])*((Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Sv_list(v) +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Sv_list(v)  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl)%x%Sv_list(v) +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)%x%Sv_list(v)
      }else{

        Pv2_term <- function(v) exp(v[pen_start+5])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)%x%Sv_list(v)) +
          exp(v[pen_start+6])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL)%x%Sv_list(v))  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl)%x%Sv_list(v) +
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)%x%Sv_list(v)

      }


      logpv2<- function(v)  {
        value <- 0.5 * nu * (v[pen_start+5]+v[pen_start+6]+v[pen_start+7]+v[pen_start+8])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+8]))+log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+6]))+
                          log(b + 0.5*nu*exp(v[pen_start+7])))+
          0.5*determinant_Pv(Pv2_term(v)) +
          a.rho*v[pen_start+9] - (a.rho + b.rho)*log(1 + exp(v[pen_start+9]))
        return(as.numeric(value))}


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      v_init <- c(v_init, c(rep(3,4),1))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-5)],
                         paste("Penalty spatial DLNM",deparse1(all_terms[[term]]), 1:4),
                         paste("Rho DLNM",deparse1(all_terms[[term]])))

      pen_dlnm[term+1] <- pen_dlnm[term+1]+5

      constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))

    }}else{
      logpv1 <- function( v) 0
      nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
      Pv1_term <- function(v){
        methods::as(Matrix::Diagonal(n = nbasis,
                            x = zeta), "generalMatrix")
      }

      if (sm_args_iter$arg.by$type == "ind"){
        Sv_list <- Matrix::Diagonal(n = q.rand, x = 1)

        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Sv_list
        n_d <- q.rand*ncol(sm_args_iter$x)


        logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))

        v_init <- c(v_init, 5)
        names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])))

        logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
        Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+1

        constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))

      } else if (sm_args_iter$arg.by$type == "ICAR"){
        neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
        Rn <- make_graph_precision(neig.map, q.rand)

        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Rn
        n_d <- q.rand*ncol(sm_args_iter$x)



        logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))

        v_init <- c(v_init, 5)
        names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])))

        logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
        Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+1

        ncol_Wcb <- ncol(Wcb_part)
        constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:ncol_Wcb, each = q.rand),
                                                           j =  (ncol_Wcb+1):(ncol_Wcb*(1+q.rand)),
                                                           x = 1,
                                                           dims = c(ncol_Wcb, ncol_Wcb*(1+q.rand)))

      }else if (sm_args_iter$arg.by$type == "Leroux"){
        neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
        Rn <- make_graph_precision(neig.map, q.rand)


        Sv_list <- function(v) {
          make_leroux_precision(0, v[pen_start + 2], q.rand, Rn)
        }

          Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%Sv_list(v)


        logpv2<- function(v)  {
          value <- 0.5 * nu * (v[pen_start+1])-
            (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))+
            0.5*determinant_Pv(Pv2_term(v)) +
            a.rho*v[pen_start+2] - (a.rho + b.rho)*log(1 + exp(v[pen_start+2]))
          return(as.numeric(value))}


        logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
        Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

        v_init <- c(v_init, c(5,1))
        names(v_init) <- c(names(v_init)[1:(length(v_init)-2)],
                           paste("Precision spatial DLNM",deparse1(all_terms[[term]])),
                           paste("Rho DLNM",deparse1(all_terms[[term]])))

        pen_dlnm[term+1] <- pen_dlnm[term+1]+2

        constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Wcb_part)*(1+q.rand))

      }}

  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)
}




#' Construction of effect modified joint DLNM components
#' @noRd
build_component_modified2Ddlnm<-function(state, sm_args_iter, term, all_terms, constants){

  Wcb <- state$Wcb
  pen_dlnm <- state$pen_dlnm
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed
  v_init <- state$v_init
  sm_list <- state$sm_list

  constraint_matrix  <- state$constraint

  zeta<- constants$zeta
  a <- constants$a
  b <- constants$b
  nu <- constants$nu

  data <- constants$data
  env <- constants$env

  # main effect
  crossbasis_em <- Matrix::Matrix(sm_args_iter$x, sparse = TRUE)

  colnames(crossbasis_em) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(crossbasis_em), sep = ""
  )

  Wcb <- add_wcb_block(Wcb, crossbasis_em)
  penorder <- sm_args_iter$arg$pen

  pen_start <- cumsum(pen_dlnm)[term]

  #modifying effect
  if(sm_args_iter$arg.by$type == "linear"){
    sm_list[[term]] = matrix(sm_args_iter$by, ncol = 1)
    vz <- 1


    colnames(sm_list[[term]]) <- deparse(normalize_sm_call(all_terms[[term]])$by)


  }else if (sm_args_iter$arg.by$type == "factor"){
    sm_list[[term]] = Matrix::sparse.model.matrix(~ as.factor(sm_args_iter$by))[,-1, drop = F]
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- sub(
      "^sm_args_iter\\$by",
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]])
    )

  }else if(sm_args_iter$arg.by$type == "ps"){
    sm_list[[term]] <- dlnm::ps(sm_args_iter$by,
                                diff = sm_args_iter$arg.by$pen,
                                df = sm_args_iter$arg.by$df,
                                intercept = F)

    D_z = Matrix::Matrix(attr(sm_list[[term]],"S"))+
      Matrix::Diagonal(n = ncol(sm_list[[term]]), x = 1e-12)

    vz <- ncol(sm_list[[term]])

    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }else if (sm_args_iter$arg.by$type == "ns"){
    sm_list[[term]] <- splines::ns(sm_args_iter$by,
                                   df = sm_args_iter$arg.by$df,
                                   knots = sm_args_iter$arg.by$knots,
                                   intercept =F,
                                   Boundary.knots = sm_args_iter$arg.by$Boundary.knots)
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }else{
    sm_list[[term]] <- splines::bs(sm_args_iter$by,
                                   df = sm_args_iter$arg.by$df,
                                   knots = sm_args_iter$arg.by$knots,
                                   degree = sm_args_iter$arg.by$degree,
                                   intercept = F,
                                   Boundary.knots = sm_args_iter$arg.by$Boundary.knots)
    vz <- ncol(sm_list[[term]])


    colnames(sm_list[[term]]) <- paste(
      deparse(normalize_sm_call(all_terms[[term]])$by),
      colnames(sm_list[[term]]), sep = ""
    )

  }


  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(crossbasis_em), Matrix::Matrix(Matrix::t(sm_list[[term]]))))


  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(crossbasis_em), colnames(sm_list[[term]]), paste, sep = ":")
  ))


  # Global design matrix
  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(crossbasis_em)+ncol(interaction_matrix))

  Wcb <- add_wcb_block(Wcb, interaction_matrix)
  rm(interaction_matrix, crossbasis_em)



  if(!is.null(penorder)){

    vx1 <- attributes(sm_args_iter$x)$argvar1$df
    vx2 <- attributes(sm_args_iter$x)$argvar2$df
    vl <- attributes(sm_args_iter$x)$arglag$df

    # Penalty for exposure 1
    Dx1 <- Matrix::Diagonal(vx1, x = NULL)
    for (k in 1:penorder) Dx1 <- Matrix::diff(Dx1)
    Px1 <- Matrix::t(Dx1) %*% Dx1
    Px1 <- Px1 + Matrix::Diagonal(vx1, x = 1e-12)


    # Penalty for exposure 2
    Dx2 <- Matrix::Diagonal(vx2, x = NULL)
    for (k in 1:penorder) Dx2 <- Matrix::diff(Dx2)
    Px2 <- Matrix::t(Dx2) %*% Dx2
    Px2 <- Px2 + Matrix::Diagonal(vx2, x = 1e-12)

    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)

    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    Pl_add <- Dl_add + diag(1e-12, vl)


    if(attributes(sm_args_iter$x)$is_identifiable){


      Pv1_term <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)] +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)]  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)
    }else{

      Pv1_term <- function(v) exp(v[pen_start+1])*(Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL)) +
        exp(v[pen_start+2])*(Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))  +
        exp(v[pen_start+3])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl) +
        exp(v[pen_start+4])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)

    }


    logpv1 <- function(v) 0.5 * nu * (v[pen_start + 1]+v[pen_start + 2]+v[pen_start + 3]+v[pen_start + 4])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1]))+log(b + 0.5*nu*exp(v[pen_start + 2]))+
                      log(b + 0.5*nu*exp(v[pen_start + 3]))+log(b + 0.5*nu*exp(v[pen_start + 4])))+
      0.5*determinant_Pv(Pv1_term(v))

    pen_dlnm[term+1] = pen_dlnm[term+1]+4

    v_init <- c(v_init, rep(2, 4))
    names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                       paste("Penalty common DLNM",deparse1(all_terms[[term]]), 1:4))



    if(sm_args_iter$arg.by$type == "ps"){

      if(attributes(sm_args_iter$x)$is_identifiable){


        Pv2_term <- function(v) exp(v[pen_start+5])*((Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Matrix::Diagonal(n = vz, x = NULL) +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Matrix::Diagonal(n = vz, x = NULL)  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+9])*(Matrix::Diagonal(n = (vx1*vx2-1)*vl, x = NULL)%x%D_z)
      }else{

        Pv2_term <- function(v) exp(v[pen_start+5])*(Px1%x%(Matrix::Diagonal(n = vx2*vl*vz, x = NULL)))+
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl*vz, x = NULL))) +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+9])*(Matrix::Diagonal(n = vx1*vx2*vl, x = NULL)%x%D_z)

      }


      logpv2 <- function(v) 0.5 * nu * (v[pen_start+5]+v[pen_start+6]+v[pen_start+7]+v[pen_start+8]+
                                          +v[pen_start+9])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+6]))+
                        log(b + 0.5*nu*exp(v[pen_start+9]))+log(b + 0.5*nu*exp(v[pen_start+7]))+log(b + 0.5*nu*exp(v[pen_start+8])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 5
      v_init <- c(v_init, rep(2, 5))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-5)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]]), 1:5))

    }else{
      if(attributes(sm_args_iter$x)$is_identifiable){


        Pv2_term <- function(v) exp(v[pen_start+5])*((Px1%x%Matrix::Diagonal(n = vx2*vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Matrix::Diagonal(n = vz, x = NULL) +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl, x = NULL))[-c(1:vl),-c(1:vl)])%x%Matrix::Diagonal(n = vz, x = NULL)  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n =vx1*vx2-1, x = NULL)%x%Pl)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2-1, x = NULL)%x%Pl_add)%x%Matrix::Diagonal(n = vz, x = NULL)
      }else{

        Pv2_term <- function(v) exp(v[pen_start+5])*(Px1%x%(Matrix::Diagonal(n = vx2*vl*vz, x = NULL))) +
          exp(v[pen_start+6])*((Matrix::Diagonal(n = vx1, x = NULL)%x%Px2%x%Matrix::Diagonal(n = vl*vz, x = NULL)))  +
          exp(v[pen_start+7])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl)%x%Matrix::Diagonal(n = vz, x = NULL)+
          exp(v[pen_start+8])*(Matrix::Diagonal(n = vx1*vx2, x = NULL)%x%Pl_add)%x%Matrix::Diagonal(n = vz, x = NULL)

      }


      logpv2 <- function(v) 0.5 * nu * (v[pen_start+5]+v[pen_start+6]+v[pen_start+7]+v[pen_start+8])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+5]))+log(b + 0.5*nu*exp(v[pen_start+8]))+
                        log(b + 0.5*nu*exp(v[pen_start+6]))+log(b + 0.5*nu*exp(v[pen_start+7])))+
        0.5*determinant_Pv(Pv2_term(v))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 4
      v_init <- c(v_init, rep(2, 4))
      names(v_init) <- c(names(v_init)[1:(length(v_init)-4)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]]), 1:4))
    }




  }else{
    logpv1 <- function(v) 0
    nbasis <- ncol(eval_expression(normalize_sm_call(all_terms[[term]])$x,data,env))
    Pv1_term <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }

    if(sm_args_iter$arg.by$type == "ps"){
        Pv2_term <- function(v) methods::as(exp(v[pen_start+1])*(Matrix::Diagonal(n = ncol(sm_args_iter$x), x = NULL)), "generalMatrix")%x%D_z
        n_d <- ncol(sm_args_iter$x)*ncol(D_z)


      logpv2 <- function(v) 0.5 * (nu+n_d) * (v[pen_start+1])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))


      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))

      pen_dlnm[term+1] = pen_dlnm[term+1] + 1
      v_init <- c(v_init, 2)
      names(v_init) <- c(names(v_init)[1:(length(v_init)-1)],
                         paste("Penalty modified DLNM",deparse1(all_terms[[term]])))

    }else{

      logpv2 <- function(v) 0
      Pv2_term <- function(v){
        methods::as(Matrix::Diagonal(n = ncol(sm_args_iter$x)*vz,
                            x = zeta), "generalMatrix")
      }

      logpv.fixed[[term]] <- function(v) logpv1(v) + logpv2(v)
      Pv.list[[term]] <- function(v) Matrix::bdiag(Pv1_term(v), Pv2_term(v))


    }

  }




  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$sm_list<-sm_list
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)
}

