#' @name northplatte
#' @aliases northplatte
#' @title Exemplary sounding dataset - sample from LBF North Platte (WMO ID: 72562) - 03 July 1999, 00:00 UTC
#'
#' @description The object contains pre-downloaded sounding dataset from University of Wyoming sounding database.
#' Dataset can be downloaded with the following syntax: 
#' northplatte = get_sounding(wmo_id = 72562, yy = 1999, mm = 7, dd = 3, hh = 00)
#'
#' @importFrom utils data
#' 
#' @format A data frame with 71 rows and 6 variables as described in `get_sounding()`
#'
#' \describe{
#'   \item{pressure}{pressure [hPa]}
#'   \item{altitude}{altitude [m]}
#'   \item{temp}{temperature [degree Celsius]}
#'   \item{dpt}{dew point temperature [degree Celsius]}
#'   \item{wd}{wind direction [azimuth as degrees]}
#'   \item{ws}{wind speed [knots]}
#' }
#' 
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @docType data
#' @keywords datasets thunder
#' @usage data("northplatte")
#' @examples
#' data(northplatte)
#' head(northplatte)
"northplatte"
