#' Simulated temperature and mortality data for Barcelona districts
#'
#' A simulated dataset containing daily temperature, population,
#' deprivation, and mortality data for 10 districts of Barcelona during
#' the summer period from 2007 to 2011.
#'
#' The dataset is provided solely for illustrative and demonstration
#' purposes, in particular for demonstrating the methods implemented
#' in the \code{DLNMLPS} package. The data are simulated and do not
#' represent actual observations from Barcelona.
#'
#' @format A data frame with 6100 observations and 9 variables:
#' \describe{
#'   \item{date}{Date of observation. Data cover 1 June to 30 September
#'   for each year from 2007 to 2011.}
#'   \item{Tmax}{Daily maximum temperature in degrees Celsius.}
#'   \item{Tmin}{Daily minimum temperature in degrees Celsius.}
#'   \item{Tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{District}{District identifier, ranging from 01 to 10.}
#'   \item{pop}{Population size of the district.}
#'   \item{depriv}{Scaled deprivation index of the district.}
#'   \item{dtot}{Daily mortality count in the district.}
#' }
#'
#' @details
#' The dataset represents the 10 Barcelona districts
#' observed daily during the summer period (1 June--30 September) over
#' five years (2007--2011). Temperature variables vary by date, while
#' population size and deprivation index are district-level
#' characteristics.
#'
#' The data were generated to provide a realistic example for
#' demonstrating distributed lag non-linear modelling and related
#' methods. They should not be interpreted as real epidemiological or
#' climatological observations.
#'
#' @source Simulated data generated for the \code{DLNMLPS} package.
#'
#' @examples
#' data(barcelona_temp)
#'
#' head(barcelona_temp)
#'
#' # Inspect the temperature variables
#' summary(barcelona_temp[, c("Tmax", "Tmin", "Tmean")])
#'
#' # Number of districts
#' unique(barcelona_temp$District)
#'
#' # Plot mean temperature over time
#' plot(barcelona_temp$date, barcelona_temp$Tmean,
#'      xlab = "Date", ylab = "Mean temperature (°C)",
#'      pch = 16, cex = 0.5)
"barcelona_temp"
