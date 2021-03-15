#' Plot vertical wind speed profile
#' 
#' Function for plotting wind speed profile. 
#' Can be launched as standalone function or coupled with pre-drawn Skew-T diagram.
#'
#' @import stats
#' @import utils
#' @importFrom dplyr left_join
#' @importFrom grDevices colorRampPalette
#'
#' @param pressure pressure [hPa] 
#' @param ws wind speed [knots]
#' @param ptop pressure top level to be used for plotting wind speed. Valid options should be < 200 hPa (100 by default)
#' @param yaxs logic. Whether to add labels to heights on Y lab
#' @param ... extra graphic arguments
#' @export
#' 
#' @return vertical wind profile plot
#' 
#' @examples 
#' # load examplary dataset:
#' data("sounding_vienna")
#' attach(sounding_vienna)
#' sounding_wind(pressure = pressure, ws = ws, yaxs = TRUE)

sounding_wind <- function(pressure, ws, ptop = 100, yaxs = TRUE, ...){
        
        if(ptop > 200) {
                stop("\nptop argument needs to be set < 200 (hPa)!")
        }
        
        # convert wind speed from knots to m/s
        ws <- ws * 0.51444

        
        # define plotting area limits:
        ymax <- skewty(1050)
        #ymin <- skewty(50)
        ymin <- skewty(ptop)
        xmin <- 0
        xmax <- 80
        
        ws_units <- seq(0, xmax, by = 10)
        
        xc <- c(xmin, xmin, xmax, xmax, xmin)
        yc <- c(ymin, ymax, ymax, ymin, ymin)
        plot(xc, yc, type = "l", axes = FALSE, xlab = "", ylab = "", lwd = 1)
        axis(1, at = ws_units, labels = ws_units, pos = ymax)
        axis(3, at = ws_units, labels = NA, pos = ymin)
        mtext(side = 1, line = 1.2, paste("Wind speed (m/s)"), cex=0.8)
        
        segments(x0 = ws_units, y0 = ymax, x1 = ws_units, y1 = ymin, lwd = 0.5, col = "black", lty = 3)
        prs <- c(1050, 1000, 850, 700, 500, 400, 300, seq(from = 200, to = ptop, by = -50))
        Npressure <- length(prs)
        
        ypos <- skewty(prs[2:Npressure])
        
        if(yaxs){
                axis(2, at = ypos, labels = prs[2:Npressure], pos = xmin, padj = 1)
                # commented label for X-axis
                mtext(side = 2, line = 1.5, "pressure (hPa)", padj = 2, cex=0.8)
                
        } else {
                axis(2, at = ypos, labels = NA, pos = xmin, padj = 1)
        }
        
        
        # end of drawing diagram in square
        
        # draw pressuresure heights:
        y <- skewty(pressure)
        ind = y < ymin # clipping to max visible height
        x <- rep(xmin, length(y))
        
        
        # points(x = ws[ind], y = y[ind],  cex = 1, pch = 4, col = "black")
        # segments(x0 = 0.5, y0 = y, x1 = ws, y1 = y, col = data$cols, lwd = 5)
        
        
        data = approx(y = ws, x = y, n = 300)
        data = data.frame(x = data$y, y = data$x)
        #points(data$x, data$y, col = "blue")
        
        head(data)
        
        cols = c("white", "#c7ebfd", "#5689bc", "#62448e", "#007122", "#00ae00", "#b9d400", "#ffe300",
                 "#ffaf00", "#ff4400", "#b70000", "#73006a", "#f659f7", "#b1b1b1", "#fe7277")
        
        str(cols)
        
        cols = data.frame(x1 = 1:80, 
                          cols = paste0(colorRampPalette(cols)(80), "90"), 
                          stringsAsFactors = FALSE)
        
        
        data$x1 = round(data$x)
        
        data = dplyr::left_join(data, cols)
        
       
        # clipping data beyond the ptop
        data = data[data$y < ymin, ]
        
        #segments(x0 = 0.5, y0 = data$y, x1 = data$x, y1 = data$y, col = data$cols, lwd = 5)
        
        #lines(x = df$ws, y = skewty(df$pressuresure), lwd = 5, col = "white")
        
        
        # wariant no. 2 z wypelnieniem gradientem
        head(data)
        
        # mean Y height
        my = mean(diff(data$y))
        
        # gradient = function(x, y, col){
        #         y = y + 0.5*my
        #         polygon(x = c(0, x, x, 0), y = c(y - my, y - my, y, y), col = col, border = NA)
        # }
        
        # for (i in 1:nrow(data)){
        #         print(i)
        #         gradient(data$x[i]-0.1, data$y[i], col = data$cols[i])
        # }
        
        # improve layout
        lines(data$x, data$y, lwd = 2.5)
        #segments(x0 = ws_units, y0 = ymax, x1 = ws_units, y1 = ymin, lwd = 0.5, col = "black", lty = 3)
        #points(x = ws[ind], y = y[ind],  cex = 1, pch = 19, col = "black")
        lines(x = ws[ind], y = y[ind],  cex = 1, pch = 19, col = "black")
        
        
        # approach no. 2 for checking purposes
        for (i in 1:nrow(cols)){
                polygon(x = c(cols$x1[i]-1, cols$x1[i], cols$x1[i], cols$x1[i]-1),
                        y = c(ymax, ymax, ymin, ymin), 
                        col = cols$cols[i], border = NA)
        }
        
        polygon(x = c(xmax, data$x[1], data$x, data$x[length(data$x)], xmax), 
                y = c(ymax, ymax, data$y, ymin, ymin), 
                col = "white", border = NA)
        
        lines(data$x, data$y, lwd = 2.5, col = "#00000080")
        segments(x0 = ws_units, y0 = ymax, x1 = ws_units, y1 = ymin, lwd = 0.5, col = "black", lty = 3)
        points(x = ws[ind], y = y[ind],  cex = 0.8, pch = 4, col = "#00000095")
        #lines(x = ws[ind], y = y[ind],  cex = 1, pch = 19, col = "black")
        
        
}

