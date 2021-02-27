#' @name sounding_wien
#' @aliases sounding_wien
#' @title Examplary sounding dataset - sample from Wien (Austria) - 2011/08/23, 1200 UTC
#'
#' @description The object contains pre-downloaded sounding dataset as obtained for Vienna rawinsonde station
#' from http://weather.uwyo.edu/upperair/ with the use of the R climate package.
#' The snapshot was taken on 2021/02/27 with the following syntax: 
#' demo_dataset = get_sounding(wmo_id = 11035, yy = 2011, mm = 8, dd = 23, hh = 12)        
#'
#' @importFrom utils data
#' @importFrom climate sounding_wyoming
#' 
#' @format A data frame with 88 rows and 11 variables as described in `get_sounding()`
#'
#' \describe{
#'   \item{pressure}{Air pressure in hPa}
#'   \item{altitude}{Height in metres}
#'   \item{temp}{Air temperutre in degree Celsius}
#'   \item{dpt}{Dew point temperutre in degree Celsius}
#'   \item{rh}{Relative humidity in percent}
#'   \item{mixr}{Mixing ratio in g/kg}
#'   \item{wd}{Wind direction in degrees}
#'   \item{ws}{Wind speed in knots}
#'   \item{thta}{K}
#'   \item{thte}{K}
#'   \item{thtv}{K}
#' }
#' 
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @docType data
#' @keywords datasets thunder
#' @usage data("sounding_wien")
#' @examples
#' data(sounding_wien)
#' head(sounding_wien)
"sounding_wien"