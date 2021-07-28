#' Save `sounding_layout` to a graphical file
#' 
#' Auxiliary function to `sounding_plot` that plots a composite \
#' of Skew-T, hodograph and selected convective parameters \
#' on a single layout and saves as graphical file.
#' 
#' 
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as a surface, i.e h = 0 m)
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction in degrees [azimuth in degrees]
#' @param ws wind speed [knots]
#' @param title title to be added in the layout's header
#' @param parcel parcel tracing on Skew-T for "MU", "ML" or "SB" parcel
#' @param filename output file name with extension indicating file format (e.g. "my_plot.png" or "my_plot.svg")
#' @param ... other arguments that can be used with `sounding_plot` or other graphic arguments
#' @export
#' @import aiRthermo
#' @importFrom tools file_ext
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom grDevices svg
#' 
#' @return graphical file with Skew-T and hodograph layout
#' 
#' @examples
#' \donttest{
#' data("sounding_vienna")
#' sounding_save(filename = "Vienna.png", 
#'               sounding_vienna$pressure, sounding_vienna$altitude, sounding_vienna$temp, 
#'               sounding_vienna$dpt, sounding_vienna$wd, sounding_vienna$ws, 
#'               parcel = "MU", 
#'               title = "Vienna - 23 August 2011, 12:00 UTC")
#' } 
#'

sounding_save = function(pressure, altitude, temp, dpt, wd, ws,
                         title = "", parcel = "MU", filename, ...){
  
  convert = FALSE
  ptop = 100 
  
  stopifnot(length(filename) < 4)
  
  if(tools::file_ext(filename) == "png"){
    grDevices::png(filename = filename, width = 2000, height = 1200, res = 200)
    sounding_plot(pressure, altitude, temp, 
                  dpt, wd, ws, title = title, ...)
    grDevices::dev.off()
  }
  
  
  if(tools::file_ext(filename) == "svg"){
    grDevices::svg(filename = filename, width = 20, height = 12, pointsize=24)
    sounding_plot(pressure, altitude, temp, 
                  dpt, wd, ws, title = title)
    grDevices::dev.off()
  }
  
  
}
