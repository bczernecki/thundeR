#' Add line to a Skew-T diagram
#' 
#' Calculate X and Y coordinates for lines to be drawn on Skew-T diagram;
#' Draw any line on Skew-T diagram using temperature and pressure as coordinates
#'
#' 
#' @param temp coordinates to be used based on air temperature vector
#' @param pressure coordinates to be used base on air pressure vector
#' @param ptop upper limit of drawn trajectory (default: 100 hPa); use only if a line goes beyond the drawing area
#' @param ... other graphical parameters that can be passed to `lines()` function, such as `lwd`, `lty`, `col`, etc.
#' @export
#' 
#' @return adds line on a pre-defined Skew-T plot
#' 
#' @examples 
#' 
#' # take a sample sounding profile:
#' data("sounding_vienna")
#' attach(sounding_vienna)
#' 
#' # draw empty Skew-T plot:
#' skewt_plot(temp_stripes = TRUE)
#' 
#' # draw line for dew-point temperature:
#' skewt_lines(dpt, pressure, type = 'l', col = 'forestgreen', lwd = 2.5) 
#' # draw line for air temperature:
#' skewt_lines(temp, pressure, type = 'l', col='red', lwd = 2.5) 
#' 

skewt_lines = function (temp, pressure, ptop = 100, ...) {
  ind = pressure >= ptop
  v = skewty(pressure[ind]) # extra checks for NA coded as -99
  u = skewtx(temp[ind], v)
  u[u<(-27.75) | u>(26.25)] <- NA # to avoid going above plotting area (works smoothly when data every 5m is provided)
  graphics::lines(u, v,  ...)
}
