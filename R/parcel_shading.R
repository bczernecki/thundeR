#' Calculate X and Y coordinates for lines to be drawn on Skew-T diagram
#' 
#' Draw any line on Skew-T diagram using temperature and pressure as coordinates
#'
#' 
#' @param temp1 coordinates to be used based on air temperature vector
#' @param temp2 coordinates to be used based on air temperature vector
#' @param pressure1 coordinates to be used base on air pressure vector
#' @param pressure2 coordinates to be used base on air pressure vector
#' @param ptop upper limit of drawn trajectory (default: 100 hPa); use only if a line goes beyond the drawing area
#' @param ... other graphical parameters that can be passed to `graphics::lines()` function, such as `lwd`, `lty`, `col`, etc.
#' @noRd
#' 

parcel_shading <- function (temp1, pressure1, 
                            temp2, pressure2, 
                            ptop = 100, ...) {
        
        ind <- pressure1 >= ptop
        v <- skewty(pressure1[ind]) # extra checks for NA coded as -99
        u <- skewtx(temp1[ind], v)
        
        graphics::lines(u, v,  ...)
}
