# Building crossbasis special functions ----------------------------

###########################################
# COMPONENT RE
###########################################

build_component_re<-function(state, sm_args_iter, term, all_terms, constants){

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

  Z.rand <- Matrix::sparse.model.matrix(~ sm_args_iter$x + 0)

  colnames(Z.rand) <- sub(
    "^sm_args_iter\\$x",
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Z.rand)
  )

  q.rand <- dim(Z.rand)[2]

  pen_start <- cumsum(pen_dlnm)[term]

  if(sm_args_iter$arg$type == "ind"){

    Pv.list[[term]] <- function(v) Matrix::Diagonal(n = q.rand, x = exp(v[pen_start+1]))
    logpv.fixed[[term]] <- function(v) 0.5 * (nu + q.rand) * v[pen_start+1] -
      (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1]))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+1
    v_init <- c(v_init, 0)

    constraint_matrix[[term]] <- Matrix::Matrix(0,1,q.rand)

  } else if (sm_args_iter$arg$type == "ICAR"){
    neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
    Rn <- make_graph_precision(neig.map, q.rand)

    Pv.list[[term]] <- function(v) exp(v[pen_start+1])*Rn

    logpv.fixed[[term]] <- function(v) 0.5 * (nu + q.rand) * v[pen_start+1] -
      (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1]))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+1
    v_init <- c(v_init, 0)

    constraint_matrix[[term]] <-  Matrix::Matrix(1, 1, q.rand)

  } else if (sm_args_iter$arg$type == "Convolution") {
    neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
    Rn <- make_graph_precision(neig.map, q.rand)

    Pv.list[[term]] <- function(v) Matrix::bdiag(exp(v[pen_start+1])*Rn,
                                                 Matrix::Diagonal(n = q.rand, x = exp(v[pen_start+2])))

    logpv.fixed[[term]] <- function(v) sum(0.5 * (nu + q.rand) * v[(pen_start+1):(pen_start+2)]) -
      sum((0.5*nu + a)*log(b + 0.5*nu*exp( v[(pen_start+1):(pen_start+2)])))
    Z.rand <- cbind(Z.rand, Z.rand)

    constraint_matrix[[term]] <-  Matrix::Matrix(c(rep(1, ncol(Z.rand)/2),rep(0, ncol(Z.rand)/2)),
                                                 1, ncol(Z.rand))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+2
    v_init <- c(v_init, 0,1)


  } else if (sm_args_iter$arg$type == "Leroux"){
    neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
    Rn <- make_graph_precision(neig.map, q.rand)


    Pv.list[[term]] <- function(v) {
      make_leroux_precision(v[pen_start + 1], v[pen_start + 2], q.rand, Rn)
    }

    Pv_term <- Pv.list[[term]]

    logpv.fixed[[term]]<- function(v)  {
      value <- 0.5 * nu * v[pen_start+1] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1])) +
        0.5*determinant_Pv(Pv_term(v)) +
        a.rho*v[pen_start+2] - (a.rho + b.rho)*log(1 + exp(v[pen_start+2]))
      return(as.numeric(value))}

    pen_dlnm[term+1] <- pen_dlnm[term+1]+2
    v_init <- c(v_init, 0,1)

    constraint_matrix[[term]] <- Matrix::Matrix(0,1,q.rand)


  }

  # Global design matrix
  Wcb <- add_wcb_block(state$Wcb, Z.rand)

  rm(Z.rand)


  state$Wcb<-Wcb
  state$pen_dlnm <- pen_dlnm
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$v_init<-v_init
  state$constraint <- constraint_matrix

  return(state)
}



###########################################
# COMPONENT RW
###########################################

build_component_rw<-function(state, sm_args_iter, term, all_terms, constants, penorder = 2){

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


  Z.rand <- Matrix::sparse.model.matrix(~ sm_args_iter$x + 0)

  colnames(Z.rand) <- sub(
    "^sm_args_iter\\$x",
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Z.rand)
  )

  q.rand <- dim(Z.rand)[2]

  pen_start <- cumsum(pen_dlnm)[term]



  Diff_mat <- Matrix::Diagonal(q.rand, x = NULL)
  for (k in 1:penorder) Diff_mat <- Matrix::diff(Diff_mat)

  if(penorder ==1){
    if(sm_args_iter$arg$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,2),
        j = c(1, q.rand),
        x = c(1, -1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }
    constraint_matrix[[term]] <-  Matrix::Matrix(1, 1, ncol(Z.rand))

  }else if (penorder == 2){
    if(sm_args_iter$arg$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,3),
        j = c(q.rand - 1, q.rand, 1),
        x = c(1, -2, 1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }
    constraint_matrix[[term]] <- Matrix::Matrix(c(rep(1, ncol(Z.rand)),
                                                  1:ncol(Z.rand)),
                                                2, ncol(Z.rand), byrow = TRUE)

  }
  Rn <- Matrix::t(Diff_mat) %*% Diff_mat


  Pv.list[[term]] <- function(v) exp(v[pen_start+1])*Rn

  logpv.fixed[[term]] <- function(v) 0.5 * (nu + q.rand) * v[pen_start+1] -
    (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1]))

  pen_dlnm[term+1] <- pen_dlnm[term+1]+1
  v_init <- c(v_init, 0)


  # Global design matrix
  Wcb <- add_wcb_block(state$Wcb, Z.rand)

  rm(Z.rand)


  state$Wcb<-Wcb
  state$pen_dlnm <- pen_dlnm
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$v_init<-v_init
  state$constraint <- constraint_matrix

  return(state)
}


###########################################
# COMPONENT SPLINE
###########################################


build_component_spline<-function(state, sm_args_iter, term, all_terms, constants){

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


  if(sm_args_iter$arg$type == "ps"){
    sm_list[[term]] <- dlnm::ps(sm_args_iter$x,
                                diff = sm_args_iter$arg$pen,
                                df = sm_args_iter$arg$df,
                                intercept = sm_args_iter$arg$intercept)

    D_list = Matrix::Matrix(attr(sm_list[[term]],"S"))+
      Matrix::Diagonal(n = sm_args_iter$arg$df, x = 1e-12)

    n_d <- ncol(D_list)

    pen_start <- cumsum(pen_dlnm)[term]
    Pv.list[[term]] <- function(v) methods::as(methods::as(exp(v[pen_start+1])*D_list,
                                         "generalMatrix"),"CsparseMatrix")


    logpv.fixed[[term]] <- function(v) 0.5 * (nu + n_d) * v[pen_start+1]-
      (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1]))


    pen_dlnm[term+1] <- pen_dlnm[term+1] + 1
    v_init <- c(v_init,5)

  }else if (sm_args_iter$arg$type == "ns"){
    sm_list[[term]] <- splines::ns(sm_args_iter$x,
                                   df = sm_args_iter$arg$df,
                                   knots = sm_args_iter$arg$knots,
                                   intercept =sm_args_iter$arg$intercept,
                                   Boundary.knots = sm_args_iter$arg$Boundary.knots)

    logpv.fixed[[term]] <- function(v) 0

    nbasis <- ncol(sm_list[[term]])
    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }

  }else{
    sm_list[[term]] <- splines::bs(sm_args_iter$x,
                                   df = sm_args_iter$arg$df,
                                   knots = sm_args_iter$arg$knots,
                                   degree = sm_args_iter$arg$degree,
                                   intercept = sm_args_iter$arg$intercept,
                                   Boundary.knots = sm_args_iter$arg$Boundary.knots)

    logpv.fixed[[term]] <- function(v) 0

    nbasis <- ncol(sm_list[[term]])
    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(n = nbasis,
                          x = zeta), "generalMatrix")
    }
  }

  colnames(sm_list[[term]]) <- paste(
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(sm_list[[term]]), sep = ""
  )

  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(sm_list[[term]]))
  Wcb <- add_wcb_block(Wcb, sm_list[[term]])


  state$sm_list<-sm_list
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$pen_dlnm<-pen_dlnm
  state$v_init<-v_init
  state$Wcb<-Wcb

  state$constraint <- constraint_matrix

  return(state)
}



###########################################
# COMPONENT RW - spat
###########################################

build_component_spatialrw<-function(state, sm_args_iter, term, all_terms, constants, penorder = 2){

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

  # RW effect
  lev <- sort(unique(sm_args_iter$x))
  Z.rand <- Matrix::sparse.model.matrix(~ sm_args_iter$x + 0)

  colnames(Z.rand) <- sub(
    "^sm_args_iter\\$x",
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Z.rand)
  )

  q.rand <- dim(Z.rand)[2]

  #spatial effect
  Z.rand_spat <- Matrix::sparse.model.matrix(~ sm_args_iter$by + 0)
  colnames(Z.rand_spat) <- sub(
    "^sm_args_iter\\$by",
    deparse(normalize_sm_call(all_terms[[term]])$by),
    colnames(Z.rand_spat)
  )

  q.rand_spat <- dim(Z.rand_spat)[2]


  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(Z.rand_spat), Matrix::t(Z.rand)))

  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(Z.rand_spat), colnames(Z.rand), paste, sep = ":")
  ))


  # Global design matrix
  Wcb <- add_wcb_block(Wcb, interaction_matrix)

  rm(Z.rand_spat, interaction_matrix, Z.rand)



  pen_start <- cumsum(pen_dlnm)[term]



  Diff_mat <- Matrix::Diagonal(q.rand, x = NULL)
  for (k in 1:penorder) Diff_mat <- Matrix::diff(Diff_mat)

  if(penorder ==1){
    if(sm_args_iter$arg$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,2),
        j = c(1, q.rand),
        x = c(1, -1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }


  }else if (penorder == 2){
    if(sm_args_iter$arg$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,3),
        j = c(q.rand - 1, q.rand, 1),
        x = c(1, -2, 1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }

  }



  Rn <- Matrix::t(Diff_mat) %*% Diff_mat


  if (sm_args_iter$arg.by$type == "ind"){
    Sv_list <- Matrix::Diagonal(n = q.rand_spat, x = 1)

    Pv.list[[term]] <- function(v) exp(v[pen_start + 1])*(Sv_list%x%Rn)


    logpv.fixed[[term]]<- function(v) 0.5 * (nu +q.rand*q.rand_spat) * (v[pen_start + 1])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1])))

    v_init <- c(v_init, 0)

    pen_dlnm[term+1] = pen_dlnm[term+1]+1

    if(penorder == 1){
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                         j =  1:(q.rand*q.rand_spat),
                                                         x = 1,
                                                         dims = c(q.rand_spat, q.rand*q.rand_spat))
    }else{
      constraint_matrix[[term]] <-  rbind(Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = 1,
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)),
                                          Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = rep(1:q.rand, q.rand_spat),
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)))

    }

  }else if (sm_args_iter$arg.by$type == "ICAR"){
    neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
    Rn_spat <- make_graph_precision(neig.map, q.rand_spat)

    Pv.list[[term]] <- function(v) exp(v[pen_start + 1])*(Rn_spat%x%Rn)

    logpv.fixed[[term]]<- function(v) 0.5 * (nu +q.rand*q.rand_spat)* (v[pen_start + 1])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1])))


    pen_dlnm[term+1] = pen_dlnm[term+1]+1
    v_init <- c(v_init, 0)

    if(penorder == 1){
      constraint_matrix[[term]] <-
        Matrix::sparseMatrix(
          i = rep(1, q.rand * q.rand_spat),
          j = 1:(q.rand * q.rand_spat),
          x = 1,
          dims = c(1, q.rand * q.rand_spat)
        )
    }else{
      constraint_matrix[[term]] <-
        Matrix::sparseMatrix(
          i = rep(1:2, each = q.rand * q.rand_spat),
          j = rep(1:(q.rand * q.rand_spat), 2),
          x = c(
            rep(1, q.rand * q.rand_spat),
            rep(seq_len(q.rand), q.rand_spat)
          ),
          dims = c(2, q.rand * q.rand_spat)
        )

    }

  }else if (sm_args_iter$arg.by$type == "Leroux"){
    neig.map <- build_neighbor_map(sm_args_iter$arg.by$map, connect = sm_args_iter$arg.by$connect)
    Rn_spat <- make_graph_precision(neig.map, q.rand_spat)

    Sv_list <- function(v) {
      make_leroux_precision(0, v[pen_start + 2], q.rand_spat, Rn_spat)
    }

    Pv.list[[term]] <- function( v) exp(v[pen_start+1])*(Sv_list(v)%x%Rn)

    Pv_term <- Pv.list[[term]]

    logpv.fixed[[term]] <- function(v)  {
      value <- 0.5 * nu * (v[pen_start+1])-
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start+1])))+
        0.5*determinant_Pv(Pv_term(v)) +
        a.rho*v[pen_start+2] - (a.rho + b.rho)*log(1 + exp(v[pen_start+2]))
      return(as.numeric(value))}


    v_init <- c(v_init, c(0,1))
    pen_dlnm[term+1] <- pen_dlnm[term+1]+2

    if(penorder == 1){
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                         j =  1:(q.rand*q.rand_spat),
                                                         x = 1,
                                                         dims = c(q.rand_spat, q.rand*q.rand_spat))
    }else{
      constraint_matrix[[term]] <-  rbind(Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = 1,
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)),
                                          Matrix::sparseMatrix(i = rep(1:q.rand_spat, each = q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = rep(1:q.rand, q.rand_spat),
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)))

    }

  }


  state$Wcb<-Wcb
  state$pen_dlnm <- pen_dlnm
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$v_init<-v_init

  state$constraint <- constraint_matrix

  return(state)
}





###########################################
# COMPONENT spat - rw
###########################################

build_component_rwspatial<-function(state, sm_args_iter, term, all_terms, constants){

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

  # RW effect

  Z.rand <- Matrix::sparse.model.matrix(~ sm_args_iter$by + 0)

  colnames(Z.rand) <- sub(
    "^sm_args_iter\\$by",
    deparse(normalize_sm_call(all_terms[[term]])$by),
    colnames(Z.rand)
  )

  q.rand <- dim(Z.rand)[2]

  #spatial effect
  Z.rand_spat <- Matrix::sparse.model.matrix(~ sm_args_iter$x + 0)
  colnames(Z.rand_spat) <- sub(
    "^sm_args_iter\\$x",
    deparse(normalize_sm_call(all_terms[[term]])$x),
    colnames(Z.rand_spat)
  )

  q.rand_spat <- dim(Z.rand_spat)[2]


  interaction_matrix <- Matrix::t(Matrix::KhatriRao(Matrix::t(Z.rand), Matrix::t(Z.rand_spat)))

  colnames(interaction_matrix) <- as.vector(t(
    outer(colnames(Z.rand), colnames(Z.rand_spat), paste, sep = ":")
  ))


  # Global design matrix
  Wcb <- add_wcb_block(Wcb, interaction_matrix)

  rm(Z.rand_spat, interaction_matrix, Z.rand)



  pen_start <- cumsum(pen_dlnm)[term]


  if (sm_args_iter$arg.by$type == "ind"){
    Rn <- Matrix::Diagonal(n = q.rand, x = 1)
  } else if (sm_args_iter$arg.by$type == "rw1"){
    Diff_mat <- Matrix::Diagonal(q.rand, x = NULL)
    Diff_mat <- Matrix::diff(Diff_mat)
    if(sm_args_iter$arg.by$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,2),
        j = c(1, q.rand),
        x = c(1, -1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }
    Rn <- Matrix::t(Diff_mat) %*% Diff_mat
  } else if (sm_args_iter$arg.by$type == "rw2"){
    Diff_mat <- Matrix::Diagonal(q.rand, x = NULL)
    for (k in 1:2) Diff_mat <- Matrix::diff(Diff_mat)
    if(sm_args_iter$arg.by$cyclic){
      cyc_row <- Matrix::sparseMatrix(
        i = rep(1,3),
        j = c(q.rand - 1, q.rand, 1),
        x = c(1, -2, 1),
        dims = c(1, q.rand)
      )
      Diff_mat <- rbind(Diff_mat, cyc_row)
    }


    Rn <- Matrix::t(Diff_mat) %*% Diff_mat
  }


  if(sm_args_iter$arg$type == "ind"){
    Rn_spat <- Matrix::Diagonal(n = q.rand_spat, x = 1)

    Pv.list[[term]] <- function(v) exp(v[pen_start + 1])*(Rn%x%Rn_spat)


    logpv.fixed[[term]]<- function(v) 0.5 * (nu +q.rand*q.rand_spat) * (v[pen_start + 1])-
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[pen_start + 1])))

    v_init <- c(v_init, 0)

    pen_dlnm[term+1] = pen_dlnm[term+1]+1


    if (sm_args_iter$arg.by$type == "ind"){
      constraint_matrix[[term]] <-  Matrix::Matrix(0,1,q.rand*q.rand_spat)
    }else if(sm_args_iter$arg.by$type == "rw1"){
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                         j =  1:(q.rand*q.rand_spat),
                                                         x = 1,
                                                         dims = c(q.rand_spat, q.rand*q.rand_spat))
    }else{
      constraint_matrix[[term]] <-  rbind(Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = 1,
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)),
                                          Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = rep(1:q.rand, q.rand_spat),
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)))

    }



  } else if (sm_args_iter$arg$type == "ICAR"){
    neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
    Rn_spat <- make_graph_precision(neig.map, q.rand_spat)

    Pv.list[[term]] <- function(v) exp(v[pen_start+1])*(Rn%x%Rn_spat)

    logpv.fixed[[term]] <- function(v) 0.5 * (nu + q.rand*q.rand_spat) * v[pen_start+1] -
      (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1]))

    pen_dlnm[term+1] <- pen_dlnm[term+1]+1
    v_init <- c(v_init, 0)

    if (sm_args_iter$arg.by$type == "ind"){
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:q.rand, each = q.rand_spat),
                                                         j =  1:(q.rand*q.rand_spat),
                                                         x = 1,
                                                         dims = c(q.rand, q.rand*q.rand_spat))
    }else if(sm_args_iter$arg.by$type == "rw1"){
      constraint_matrix[[term]] <-
        Matrix::sparseMatrix(
          i = rep(1, q.rand * q.rand_spat),
          j = 1:(q.rand * q.rand_spat),
          x = 1,
          dims = c(1, q.rand * q.rand_spat)
        )

    }else{
      constraint_matrix[[term]] <-
        Matrix::sparseMatrix(
          i = rep(1:2, each = q.rand * q.rand_spat),
          j = rep(1:(q.rand * q.rand_spat), 2),
          x = c(
            rep(1, q.rand * q.rand_spat),
            rep(1:q.rand, each = q.rand_spat)
          ),
          dims = c(2, q.rand * q.rand_spat)
        )

    }

  } else if (sm_args_iter$arg$type == "Convolution"){
    stop("Convolution cannot be specified in combination with a by argument")


  }else if (sm_args_iter$arg$type == "Leroux"){
    neig.map <- build_neighbor_map(sm_args_iter$arg$map, connect = sm_args_iter$arg$connect)
    Rn_spat <- make_graph_precision(neig.map, q.rand_spat)

    Sv_list <- function(v) {
      make_leroux_precision(0, v[pen_start + 2], q.rand_spat, Rn_spat)
    }


    Pv.list[[term]] <- function( v) exp(v[pen_start+1])*(Rn %x%Sv_list(v))

    Pv_term <- Pv.list[[term]]

    logpv.fixed[[term]]<- function(v)  {
      value <- 0.5 * nu * v[pen_start+1] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_start+1])) +
        0.5*determinant_Pv(Pv_term(v)) +
        a.rho*v[pen_start+2] - (a.rho + b.rho)*log(1 + exp(v[pen_start+2]))
      return(as.numeric(value))}

    pen_dlnm[term+1] <- pen_dlnm[term+1]+2
    v_init <- c(v_init, 0,1)

    if (sm_args_iter$arg.by$type == "ind"){
      constraint_matrix[[term]] <-  Matrix::Matrix(0,1,q.rand*q.rand_spat)
    }else if(sm_args_iter$arg.by$type == "rw1"){
      constraint_matrix[[term]] <-  Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                         j =  1:(q.rand*q.rand_spat),
                                                         x = 1,
                                                         dims = c(q.rand_spat, q.rand*q.rand_spat))
    }else{
      constraint_matrix[[term]] <-  rbind(Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = 1,
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)),
                                          Matrix::sparseMatrix(i = rep(1:q.rand_spat, q.rand),
                                                               j =  1:(q.rand*q.rand_spat),
                                                               x = rep(1:q.rand, q.rand_spat),
                                                               dims = c(q.rand_spat, q.rand*q.rand_spat)))

    }


  }


  state$Wcb<-Wcb
  state$pen_dlnm <- pen_dlnm
  state$Pv.list<-Pv.list
  state$logpv.fixed<-logpv.fixed
  state$v_init<-v_init

  state$constraint <- constraint_matrix

  return(state)
}





###########################################
# COMPONENT FIXED
###########################################


build_component_fixed<-function(model, state,  term, all_terms, constants){

  Wcb <- state$Wcb
  Pv.list <- state$Pv.list
  logpv.fixed <- state$logpv.fixed

  constraint_matrix  <- state$constraint

  zeta <- constants$zeta
  data <- constants$data

  other_formula <- stats::as.formula(call("~", model[[2]], all_terms[[term]]))
  if(!is.null(data)){
    mod <- stats::model.frame(other_formula, data = data, na.action = stats::na.pass)
    na_ind = (Matrix::rowSums(is.na(mod))>0)

    Z.linear_part <- Matrix::sparse.model.matrix(mod, data = data, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )[, -1, drop = FALSE]

    Z.NA <- Matrix::sparseMatrix(
      i = rep(which(na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(na_ind))),
      x = rep(NA_real_, sum(na_ind)*ncol(Z.linear_part)),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )[, -1, drop = FALSE]


    Z.linear <- Z.linear + Z.NA

    p_term <- ncol(Z.linear)
    colnames(Z.linear) <- colnames(Z.linear_part)[-1]
    Wcb <- add_wcb_block(Wcb, Z.linear)

    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(
        n = p_term,
        x = zeta), "generalMatrix")
    }


  } else{ # dit deel begrijp ik nog niet
    mod <- stats::model.frame(other_formula, na.action = stats::na.pass)

    na_ind = (Matrix::rowSums(is.na(mod))>0)

    Z.linear_part <- Matrix::sparse.model.matrix(mod, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )[, -1, drop = FALSE]

    Z.NA <- Matrix::sparseMatrix(
      i = rep(which(na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(na_ind))),
      x = rep(NA_real_, sum(na_ind)*ncol(Z.linear_part)),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )[, -1, drop = FALSE]


    Z.linear <- Z.linear + Z.NA


    p_term <- ncol(Z.linear)
    colnames(Z.linear) <- colnames(Z.linear_part)[-1]
    Wcb <- add_wcb_block(Wcb, Z.linear)

    Pv.list[[term]] <- function(v){
      methods::as(Matrix::Diagonal(n = p_term,
                          x = zeta), "generalMatrix")
    }

  }

  logpv.fixed[[term]] <- function(v) 0

  constraint_matrix[[term]] <- Matrix::Matrix(0,1,ncol(Z.linear))

  state$Wcb<-Wcb
  state$Pv.list<-Pv.list
  state$logpv.fixed <-logpv.fixed

  state$constraint <- constraint_matrix

  return(state)
}


