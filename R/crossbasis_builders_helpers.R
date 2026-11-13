# Helpers building crossbasis   ------------------------------------------------

#' Construction of standardized list of smooth function elements
#' @noRd
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




#' Construction of block of design matrix
#' @noRd
add_wcb_block <- function(Wcb, block, name_block) {
  if (is.null(block)) {
    return(Wcb)
  }
  Wcb[[length(Wcb) + 1]] <- Matrix::Matrix(block, sparse = TRUE)
  names(Wcb)[length(Wcb)] <- name_block
  Wcb
}


#' Construction of neighbourhood matrix
#' @noRd
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

#' Construction of precision matrix from neighbourhood matrix
#' @noRd
make_graph_precision <- function(neig.map, n) {
  deg <- lengths(neig.map) -
    vapply(neig.map, function(x) sum(x == 0), integer(1))
  idx <- rep.int(seq_len(n), deg)
  nbrs <- unlist(neig.map, use.names = FALSE)

  ind_island <- nbrs == 0

  # Remove island indices from off-diagonal entries
  idx <- idx[!ind_island]
  nbrs <- nbrs[!ind_island]

  Matrix::sparseMatrix(
    i = c(seq_len(n), idx),
    j = c(seq_len(n), nbrs),
    x = c(deg, rep.int(-1, length(nbrs))),
    dims = c(n, n)
  )
}

#' Transformation of Leroux range parameter
#' @noRd
leroux_rho <- function(v) {
  exp_v <- exp(v)
  exp_v / (1 + exp_v)
}

#' Construction of Leroux precision matrix
#' @noRd
make_leroux_precision <- function(log_tau, logit_rho, n, Rn_term) {
  rho <- leroux_rho(logit_rho)
  exp(log_tau) * (Matrix::Diagonal(n = n, x = 1 - rho) + rho * Rn_term)
}


#' Align order of regions in data and map
#' @noRd
prepare_spatial_data <- function(map, region, id_col) {

  region_levels <- levels(factor(region))
  map_ids <- as.character(map[[id_col]])

  if (anyDuplicated(map_ids)) {
    stop("Map contains duplicated region IDs.")
  }

  if (!setequal(map_ids, region_levels)) {
    stop(
      "Map region IDs do not match the levels of the regional variable."
    )
  }

  map <- map[match(region_levels, map_ids), ]

  if (!identical(as.character(map[[id_col]]), region_levels)) {
    stop("Failed to align map with regional variable.")
  }

  map
}
