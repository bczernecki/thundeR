#' Download rawinsonde measurement
#' 
#' Download rawinsonde measurement from sounding database of the University of Wyoming. 
#' This is a wrapper for `climate::sounding_wyoming()` customized for the use with thundeR package.
#' 
#' @param wmo_id international WMO station code (e.g. 11035 for Vienna)
#' @param yy year - single number (e.g. 2010)
#' @param mm month - single number denoting month (e.g. 5)
#' @param dd day - single number denoting day (e.g. 23)
#' @param hh hour - single number denoting hour of measurement, for most stations it will 12 or 00 UTC (e.g. 0).
#' @param metadata - logical, whether to return metadata of downloaded sounding; default FALSE
#' 
#' @importFrom climate sounding_wyoming
#' @return Returns two lists with values described at: weather.uwyo.edu ; The first list contains:
#' \enumerate{
#'  \item pressure - pressure [hPa]
#'  \item altitude - altitude [metres above sea level]
#'  \item temp - temperature [degree C]
#'  \item dpt - dew point [degree C]
#'  \item wd - wind direction [azimuth as degrees]
#'  \item ws - wind speed [kn]
#'  }
#'  If metadata = TRUE then retrieved data is wrapped to a list with the second element containing metadata 
#'
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @export
#'
#' @examples 
#' \donttest{
#' # download rawinsonde profile from Leba station (WMO ID: 12120) for 20 August 2010 1200 UTC:
#' 
#'   profile = get_sounding(wmo_id = 12120, 
#'                          yy = 2010,
#'                          mm = 8, 
#'                          dd = 20, 
#'                          hh = 12)
#'   head(profile)
#'   
#' }

get_sounding = function(wmo_id, yy, mm, dd, hh, metadata = FALSE){

  # clipping to define max_hght
  sounding_data = climate::sounding_wyoming(wmo_id, yy, mm, dd, hh)
  
  colnames(sounding_data[[1]]) = c("pressure", "altitude", "temp", "dpt",
                                   "rh", "mixr", "wd", "ws", "thta", "thte", "thtv")
  
  sounding_data[[1]] = sounding_data[[1]][,c("pressure", "altitude", "temp", "dpt","wd", "ws")]
                                            
  sounding_data[[1]] = na.omit(sounding_data[[1]])
                                            
  if(!metadata){
    sounding_data = sounding_data[[1]]    
  }
                                            return(sounding_data)
}



        
