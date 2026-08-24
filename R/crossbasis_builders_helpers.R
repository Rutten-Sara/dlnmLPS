# Helpers building crossbasis   ------------------------------------------------

normalize_sm_call <- function(call) {
  # Match argument and remove function name
  mc <- match.call(definition = sm, call = call, expand.dots = FALSE)
  args <- as.list(mc)[-1]

  # Get defauls of sm
  defaults <- formals(sm)

  # Merge arguments with defaults
  out <- utils::modifyList(defaults, args)

  out
}




# Add block to design matrix
add_wcb_block <- function(Wcb, block) {
  if (is.null(block)) {
    return(Wcb)
  }
  Wcb[[length(Wcb) + 1]] <- Matrix::Matrix(block, sparse = TRUE)
  Wcb
}


# Construct neighbourhood matrix
build_neighbor_map <- function(map, connect = TRUE) {
  suppressWarnings(neig.map <- spdep::poly2nb(map))

  if (!isTRUE(connect)) {
    return(neig.map)
  }

  no_neigh <- which(spdep::card(neig.map) == 0)
  if (!length(no_neigh)) {
    return(neig.map)
  }

  suppressWarnings(centroids <- sf::st_centroid(map))
  coords <- sf::st_coordinates(centroids)

  for (i in no_neigh) {
    dists <- sqrt((coords[i, 1] - coords[, 1])^2 + (coords[i, 2] - coords[, 2])^2)
    dists[i] <- Inf
    nearest <- which.min(dists)
    neig.map[[i]] <- nearest
    neig.map[[nearest]] <- unique(c(neig.map[[nearest]], i))
  }

  neig.map
}

make_graph_precision <- function(neig.map, n) {
  deg <- lengths(neig.map)
  idx <- rep.int(seq_len(n), deg)
  nbrs <- unlist(neig.map, use.names = FALSE)

  Matrix::sparseMatrix(
    i = c(seq_len(n), idx),
    j = c(seq_len(n), nbrs),
    x = c(deg, rep.int(-1, length(nbrs))),
    dims = c(n, n)
  )
}


leroux_rho <- function(v) {
  exp_v <- exp(v)
  exp_v / (1 + exp_v)
}

make_leroux_precision <- function(log_tau, logit_rho, n, Rn_term) {
  rho <- leroux_rho(logit_rho)
  exp(log_tau) * (Matrix::Diagonal(n = n, x = 1 - rho) + rho * Rn_term)
}


