#' @name sounding_vienna
#' @aliases sounding_vienna
#' @title Examplary sounding dataset - sample from Vienna (WMO ID: 11035) - 23 August 2011, 1200 UTC
#'
#' @description The object contains pre-downloaded sounding dataset from University of Wyoming sounding database.
#' Dataset can be downloaded with the following syntax: 
#' demo_dataset = get_sounding(wmo_id = 11035, yy = 2011, mm = 8, dd = 23, hh = 12)        
#'
#' @importFrom utils data
#' 
#' @format A data frame with 88 rows and 6 variables as described in `get_sounding()`
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
#' @usage data("sounding_vienna")
#' @examples
#' data(sounding_vienna)
#' head(sounding_vienna)
"sounding_vienna"
