#' Save `sounding_layout` to a graphical file
#' 
#' Auxiliary function to `sounding_plot` that plots a composite \
#' of Skew-T, hodograph and selected convective parameters \
#' on a single layout and saves as graphical file.
#' 
#' 
#' @param pressure - air pressure [hPa]
#' @param altitude - altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0m)
#' @param temp - air temperature [degree Celsius]
#' @param dpt - dew point temperature [degree Celsius]
#' @param wd - wind direction in degrees [azimuth in degrees]
#' @param ws - wind speed [knots]
#' @param title to be added in the layout's header
#' @param parcel tracing on Skew-T for "MU", "ML" or "SB" parcel
#' @param filename output file name with extension indicating file format (e.g. "my_plot.png" or "my_plot.svg")
#' @param ... other arguments that can be used with `sounding_plot` or other graphic arguments
#' @export
#' @import aiRthermo
#' @importFrom tools file_ext
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom grDevices svg
#' 
#' @examples
#' data("sounding_wien")
#' attach(sounding_wien)
#' sounding_save(filename = "myfile.png", 
#'               pressure, altitude, temp, dpt, wd, ws, parcel = "MU", 
#'               title = "Wien - 2011/08/23, 12:00 UTC")
#' 
#'

sounding_save = function(pressure, altitude, temp, dpt, wd, ws,
                         title = "", parcel = "MU", filename, ...){
  
  convert = FALSE
  ptop = 100 
  interpolate = TRUE
  
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
