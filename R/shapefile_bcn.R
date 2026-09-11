#' Spatial boundaries of the 10 Barcelona districts
#'
#' Spatial boundaries of the 10 administrative districts of Barcelona,
#' provided as an \code{sf} object. The dataset is used alongside the
#' simulated \code{barcelona_temp} dataset to facilitate spatial
#' visualisation and examples in the \code{DLNMLPS} package.
#'
#' @format An \code{sf} object with 10 features and 2 columns:
#' \describe{
#'   \item{District}{District identifier, ranging from 01 to 10.}
#'   \item{geometry}{Spatial geometry of each Barcelona district,
#'   represented as \code{POLYGON} or \code{MULTIPOLYGON} features.}
#' }
#'
#' The geometries use the ETRS89 / UTM zone 31N coordinate reference
#' system (EPSG:25831).
#'
#' @details
#' The district boundaries correspond to the 10 administrative
#' districts of Barcelona and can be linked to the \code{District}
#' variable in \code{barcelona_temp}.
#'
#' Unlike \code{barcelona_temp}, which contains simulated data, the
#' spatial boundaries in this dataset originate from the Barcelona
#' City Council Open Data BCN portal.
#'
#' @source
#' Barcelona City Council, Open Data BCN.
#' \emph{Districtes i barris}.
#' \url{https://opendata-ajuntament.barcelona.cat/data/dataset/20170706-districtes-barris/resource/11851135-6919-4dcb-91ed-821e5e87a428}
#'
#' @seealso
#' \code{\link{barcelona_temp}}
#'
#' @examples
#' data(shapefile_bcn)
#'
#' library(sf)
#' plot(shapefile_bcn["District"])
"shapefile_bcn"
