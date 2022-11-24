#' Sounding export
#' 
#' Internal package function for exporting interpolated profile with 5 m steps
#' 
#' 
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m) altitude [meters]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' 
#' @importFrom climate sounding_wyoming
#' @export
#' 
#' @return Data frame of computed values for visualizing parcel trajectories
#'  \enumerate{
#'   \item pressure pressure [hPa]
#'   \item altitude altitude [m]
#'   \item temp temperature [degree Celsius]
#'   \item tempV virtual temperature [degree Celsius]
#'   \item dpt dew point temperature [degree Celsius]
#'   \item wd wind direction [azimuth in degrees]
#'   \item ws wind speed [knots]
#'   \item MU temperature for most unstable CAPE trajectory [degree Celsius]
#'   \item SB temperature for surface based CAPE trajectory [degree Celsius]
#'   \item ML temperature for mixed layer CAPE trajectory [degree Celsius]
#' }
#' 
#' @examples 
#' data("sounding_vienna")
#' attach(sounding_vienna)
#' skewt_plot(close_par = FALSE)
#' output = sounding_export(pressure, altitude, temp, dpt, wd, ws)
#' skewt_lines(output$dpt, output$pressure, col = 'forestgreen', lwd = 2.5)
#' skewt_lines(output$temp,output$pressure, col = 'red', lwd = 2.5)
#' skewt_lines(output$MU,output$pressure, col = "orange", lty = 1, lwd = 2)
#' skewt_lines(output$tempV,output$pressure, col = "red3", lty = 3, lwd = 1.5)


sounding_export = function(pressure, altitude, temp, dpt, wd, ws, accuracy = 3, interpolate_step = 5,
                          meanlayer_bottom_top = c(0,500), storm_motion = c(999,999,999)){

  parametry = sounding_default(pressure = pressure, altitude = altitude, 
                               temp = temp, dpt = dpt, wd = wd, ws = ws, 
                               export_profile = 1, accuracy = accuracy, interpolate_step = interpolate_step,
                              meanlayer_bottom_top = meanlayer_bottom_top, storm_motion = storm_motion)
  


return(parametry)
  
}
