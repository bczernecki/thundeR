#' Save `sounding_layout` to a graphical file
#' 
#' Auxiliary function to `sounding_plot` that plots a composite \
#' of Skew-T, hodograph and other relavant rawindsonde-derived profiles \
#' on a single layout and saves it into graphical file. \
#' Useful if graphical window is too small.
#' 
#' 
#' @param pressure - air pressure (hPa)
#' @param altitude - in metres
#' @param temp - air temperature (degree Celsius)
#' @param dpt - dew point temperature (degree Celsius)
#' @param wd - wind direction in degrees (0-360)
#' @param ws - wind speed in m/s
#' @param convert logical. Whether to convert wind speed from knots to m/s (default FALSE)
#' @param ptop Pressure top level to be used for plotting wind speed. Valid options should be < 200 hPa (100 by default)
#' @param interpolate logical, draw wind barbs only at interpolated altitudes with 1 km interval (default = TRUE)  instead of all wind barbs for a given input dataset
#' @param title title to be added in the layout's header
#' @param parcel tracing for "MU", "ML" or "SB" parcel
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
#' sounding_wien = na.omit(sounding_wien)
#' attach(sounding_wien)
#' sounding_save(filename = "myfile.png", 
#'               PRES, HGHT, TEMP, DWPT, DRCT, SKNT, parcel = "MU", 
#'               title = "Wien - 2011/08/23, 12:00 UTC")
#' 
#'

sounding_save = function(pressure, altitude, temp, dpt, wd, ws,
                         convert = FALSE, ptop = 100, interpolate = TRUE,
                         title = "", parcel = "MU", filename, ...){
        
        stopifnot(length(filename) < 4)
        
        if(tools::file_ext(filename) == "png"){
                grDevices::png(filename = filename, width = 980, height = 600, pointsize = 17)
                sounding_plot(pressure, altitude, temp, 
                              dpt, wd, ws, title = title, ...)
                grDevices::dev.off()
        }
        
        
        if(tools::file_ext(filename) == "svg"){
                grDevices::svg(filename = filename, width = 9.80, height = 6.00, pointsize = 12)
                sounding_plot(pressure, altitude, temp, 
                              dpt, wd, ws, title = title)
                grDevices::dev.off()
        }
        
                
}