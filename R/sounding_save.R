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
#' @param parcel parcel tracing on Skew-T for "MU", "ML" or "SB" parcel, "none" for no parcel line.
#' @param buoyancy_polygon logical, plotting area of parcel's positive (yellow) or negative (red) buoyancy (default  = TRUE)
#' @param SRH_polygon draws polygon for storm-relative helicity, available options are "0500m", "01km", "03km", "36km", "none", "03km" used as default
#' @param DCAPE draws downdraft parcel and polygon of downdraft's negative buoyancy (default = FALSE) 
#' @param filename output file name with extension indicating file format (e.g. "my_plot.png" or "my_plot.svg")
#' @param max_speed range of the hodograph to be drawn, 25 m/s used as default
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
#' attach(sounding_vienna)
#' sounding_save(filename = "Vienna.png", 
#'              pressure, altitude, temp, dpt, wd, ws, parcel = "MU", 
#'              title = "Vienna - 23 August 2011, 12:00 UTC")
#' } 
#'

sounding_save = function(pressure, altitude, temp, dpt, wd, ws,
                         title = "", parcel = "MU", max_speed = 25, buoyancy_polygon = TRUE,
                         SRH_polygon = "03km", DCAPE = FALSE, filename, ...) {
  
  convert = FALSE
  ptop = 100 
  
  stopifnot(length(filename) < 4)
  
  if (tools::file_ext(filename) == "png") {
    grDevices::png(filename = filename, width = 2000, height = 1200, res = 200)
    sounding_plot(pressure, altitude, temp, dpt, wd, ws, 
                  title, parcel, max_speed, buoyancy_polygon, SRH_polygon, DCAPE, ...)
    grDevices::dev.off()
  }

  if (tools::file_ext(filename) == "svg") {
    grDevices::svg(filename = filename, width = 20, height = 12, pointsize = 24)
    sounding_plot(pressure, altitude, temp, dpt, wd, ws, 
                  title, parcel, max_speed, buoyancy_polygon, SRH_polygon, DCAPE, ...)
    grDevices::dev.off()
  }
   
}
