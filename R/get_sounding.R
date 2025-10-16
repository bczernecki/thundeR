#' Download rawinsonde measurement
#' 
#' Download rawinsonde measurement from sounding database of the University of Wyoming in a form convenient to use with thundeR package.
#' In case of problems with downloading the chosen dataset the url is checked 5 times in 5-second intervals.
#' 
#' @param wmo_id international WMO station code (e.g. 11035 for Vienna)
#' @param yy year - single number (e.g. 2010)
#' @param mm month - single number (e.g. 5)
#' @param dd day - single number (e.g. 23)
#' @param hh hour - single number (e.g. 0)
#' @param metadata - logical, whether to return metadata of downloaded sounding; default FALSE
#' 
#' @return Returns two lists with values described at: weather.uwyo.edu ; The first list contains:
#' \enumerate{
#'  \item pressure - pressure [hPa]
#'  \item altitude - altitude [meters]
#'  \item temp - temperature [degree Celsius]
#'  \item dpt - dew point temperature [degree Celsius]
#'  \item wd - wind direction [azimuth in degrees]
#'  \item ws - wind speed [knots]
#'  }
#'  If metadata = TRUE then retrieved data is wrapped into a second list containing available metadata 
#'
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @export
#'
#' @examples 
#' \donttest{
#' # download rawinsonde profile from Vienna (WMO ID: 11035) for 23 August 2011 1200 UTC:
#' 
#'   profile = get_sounding(wmo_id = 11035, 
#'                          yy = 2011,
#'                          mm = 8, 
#'                          dd = 23, 
#'                          hh = 12)
#'   head(profile)
#'   
#' }

get_sounding = function(x) {
  return(x*x)
}
