#' @name chanhassen
#' @aliases chanhassen
#' @title Examplary sounding dataset - sample from Chanhassen (WMO ID: 72649) - 10 May 2001, 00:00 UTC
#'
#' @description The object contains pre-downloaded sounding dataset from University of Wyoming sounding database.
#' Dataset can be downloaded with the following syntax: 
#' chanhassen = get_sounding(wmo_id = 72649, yy = 2001, mm = 5, dd = 10, hh = 00)
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
#' @usage data("schanhassen")
#' @examples
#' data(chanhassen)
#' head(chanhassen)
"chanhassen"
