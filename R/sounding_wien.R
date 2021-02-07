#' @name sounding_wien
#' @aliases sounding_wien
#' @title Examplary sounding dataset - sample from Vienna rawinsonde station (2011/08/23, 1200 UTC)
#'
#' @description The object contains pre-downloaded sounding dataset as obtained for Vienna rawinsonde station
#' from http://weather.uwyo.edu/upperair/ with the use of the R climate package.
#' The snapshot was taken on 2021/01/01 with the following syntax: 
#' demo_dataset = climate::sounding_wyoming(wmo_id = 11035, yy = 2011, mm = 8, dd = 23, hh=12)[[1]] 
#'
#' @importFrom utils data
#' 
#' @format A data frame with 87 rows and 11 variables:
#' \describe{
#'   \item{PRES}{Air pressure in hPa}
#'   \item{HGHT}{Height in metres}
#'   \item{TEMP}{Air temperutre in degree Celsius}
#'   \item{DWPT}{Dew point temperutre in degree Celsius}
#'   \item{RELH}{Relative humidity in percent}
#'   \item{MIXR}{Mixing ratio in g/kg}
#'   \item{DRCT}{Wind direction in degrees}
#'   \item{SKNT}{Wind speed in knots}
#'   \item{THTA}{K}
#'   \item{THTE}{K}
#'   \item{THTV}{K}
#' }
#'
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @docType data
#' @keywords datasets sounding
#' @usage data("sounding_wien")
#' @examples
#' data(sounding_wien)
#' head(sounding_wien)
"sounding_wien"