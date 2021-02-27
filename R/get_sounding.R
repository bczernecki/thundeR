#' Download observational sounding dataset
#' 
#' Download rawinsonde dataset for vertical profile of atmosphere from University of Wyoming's database. 
#' This is a wrapper for `climate::sounding_wyoming()` customized for easier use within thundeR,
#' such as customizing headers to align to naming convention used or dealing with missing data.
#' 
#' @param wmo_id international WMO station code (World Meteorological Organization ID)
#' @param yy year - single number
#' @param mm month - single number denoting month
#' @param dd day - single number denoting day
#' @param hh hour - single number denoting initial hour of sounding; for most stations this measurement is done twice a day (i.e. at 12 and 00 UTC), sporadically 4 times a day
#' @param metadata - logical, whether to return metadata of downloaded sounding; default FALSE
#' 
#' @importFrom climate sounding_wyoming
#' @return Returns two lists with values described at: weather.uwyo.edu ; The first list contains:
#' \enumerate{
#'  \item pressure - Air Pressure (hPa)
#'  \item altitude - Altitude (metres AMSL)
#'  \item temp - Temperature (degree C)
#'  \item dpt - Dew point (C)
#'  \item rh - Relative humidity (%)
#'  \item mixr - Mixing ratio (g/kg)
#'  \item wd - Wind direction (deg)
#'  \item ws - Wind speed (knots)
#'  \item thta - THETA-A (K)
#'  \item thte - Equivalent Potential Temperature (THETA-E) (K)
#'  \item thtv - THETA-V (K)
#'  }
#'  If metadata = TRUE then retrieved data is wrapped to a list with the second one containing metadata 
#'
#' @source http://weather.uwyo.edu/upperair/sounding.html
#' @export
#'
#' @examples 
#' \donttest{
#' # generate the date to download a random dataset for 12120 station (LEBA, PL):
#' 
#'   profile = get_sounding(wmo_id = 12120, 
#'                          yy = sample(2000:2019,1),
#'                          mm = sample(1:12,1), 
#'                          dd = sample(1:20,1), 
#'                          hh = 0)
#'   head(profile)
#'   
#' }

get_sounding = function(wmo_id, yy, mm, dd, hh, metadata = FALSE){
        
        # clipping to define max_hght
        sounding_data = climate::sounding_wyoming(wmo_id, yy, mm, dd, hh)
        
        colnames(sounding_data[[1]]) = c("pressure", "altitude", "temp", "dpt",
                                         "rh", "mixr", "wd", "ws", "thta", "thte", "thtv")
        
        sounding_data[[1]] = sounding_data[[1]][-unique(c(which(is.na(sounding_data[[1]]$pressure)),
          which(is.na(sounding_data[[1]]$temp)))),]
        
        if(!metadata){
            sounding_data = sounding_data[[1]]    
        }
        
        return(sounding_data)
}


        