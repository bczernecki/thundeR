#' Wind barbs
#' 
#' A function to plot a wind barb. This is a modified version of `station.symbol` function from the RadioSonde package. 
#' Currently wind barbs are supported up to 190 knots. 
#' 
#' @import stats
#' @import utils
#'
#' @export
#' @param cx x coordinates on a plot
#' @param cy y coordinates on a plot
#' @param direction wind direction (0-360 degrees)
#' @param speed wind speed in knots
#' @param cex symbol size. Default 1
#' 
#' @examples 
#' plot(1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", frame = FALSE)
#' windbarbs(cx = 1, cy = 1, direction = 90, speed = 1, cex = 5)
#' 
#' # multiplot
#' par(mfrow=c(5,4), mar = c(1,1,1,1))
#' for (i in 19:38){
#'   sc = 5
#'   plot(0:2, xaxt = 'n', yaxt = 'n', type = "n", xlab = "", ylab = "")
#'   text(1.4,1, i*sc, cex = 1.5)
#'   windbarbs(cx = 2, cy = 1, direction = 60, speed = i*sc, cex = 3)
#' }


windbarbs = function (cx, cy, direction, speed = NA, cex = 1) {
        ns = length(cx)
        if (length(cy) != ns) {
                stop("x and y coordinates should have the same length")
        msg = "ALL VARIABLES SHOULD HAVE SAME LENGTH AS COORDINATES, OR BE MISSING!!!"
        }
        
        if (ns > 1) {
                if (length(direction) == 1)
                        if (!is.na(direction))
                                stop(msg)
                if (length(speed) == 1)
                        if (!is.na(speed))
                                stop(msg)

                if (length(direction) > 1 & length(direction) != ns)
                        stop(msg)
                if (length(speed) > 1 & length(speed) != ns)
                        stop(msg)
        }
        if (ns > 3) {
                mspd = mean(as.numeric(speed), na.rm = TRUE)
                rspd = sd(as.numeric(speed), na.rm = TRUE)
        }
        tpar = par()
        size = tpar$csi
        scalex = (tpar$usr[2] - tpar$usr[1])/(tpar$pin[1])
        scaley = (tpar$usr[4] - tpar$usr[3])/(tpar$pin[2])
        scalex = (cex * (scalex * size))/4.5
        scaley = (cex * (scaley * size))/4.5
        for (i in 1:ns) {
                spdcolor = "black"
                pscolor = "blue"
                x = cx[i]
                y = cy[i]
                if (is.na(x) | is.na(y))
                        next
                spd = speed[i]
                
                if (!is.na(spd)) {
                        xs = if (spd > 0) {
                                X1 = 0
                                X2 = 0
                                Y1 = 0
                                Y2 = 5
                                if (spd >= 5 & spd < 10) {
                                        X1 = c(X1, 0)
                                        X2 = c(X2, 1)
                                        Y1 = c(Y1, 5)
                                        Y2 = c(Y2, 5)
                                }
                                if (spd >= 10 & spd < 15) {
                                        X1 = c(X1, 0)
                                        X2 = c(X2, 2)
                                        Y1 = c(Y1, 5)
                                        Y2 = c(Y2, 5)
                                }
                                if (spd >= 15 & spd < 20) {
                                        X1 = c(X1, 0, 0)
                                        X2 = c(X2, 1, 2)
                                        Y1 = c(Y1, 4, 5)
                                        Y2 = c(Y2, 4, 5)
                                }
                                if (spd >= 20 & spd < 25) {
                                        X1 = c(X1, 0, 0)
                                        X2 = c(X2, 2, 2)
                                        Y1 = c(Y1, 4, 5)
                                        Y2 = c(Y2, 4, 5)
                                }
                                if (spd >= 25 & spd < 30) {
                                        X1 = c(X1, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2)
                                        Y1 = c(Y1, 3, 4, 5)
                                        Y2 = c(Y2, 3, 4, 5)
                                }
                                if (spd >= 30 & spd < 35) {
                                        X1 = c(X1, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2)
                                        Y1 = c(Y1, 3, 4, 5)
                                        Y2 = c(Y2, 3, 4, 5)
                                }
                                if (spd >= 35 & spd < 40) {
                                        X1 = c(X1, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2)
                                        Y1 = c(Y1, 2, 3, 4, 5)
                                        Y2 = c(Y2, 2, 3, 4, 5)
                                }
                                if (spd >= 40 & spd < 45) {
                                        X1 = c(X1, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 2, 3, 4, 5)
                                        Y2 = c(Y2, 2, 3, 4, 5)
                                }
                                if (spd >= 45 & spd < 50) {
                                        X1 = c(X1, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1, 2, 3, 4, 5)
                                        Y2 = c(Y2, 1, 2, 3, 4, 5)
                                }
                                if (spd >= 50 & spd < 55) {
                                        X1 = c(X1, 0, 0)
                                        X2 = c(X2, 2, 2)
                                        Y1 = c(Y1, 4, 5)
                                        Y2 = c(Y2, 4.5, 4.5)
                                }
                                if (spd >= 55 & spd < 60) {
                                        X1 = c(X1, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2)
                                        Y1 = c(Y1, 3, 4, 5)
                                        Y2 = c(Y2, 3, 4.5, 4.5)
                                }
                                if (spd >= 60 & spd < 65) {
                                        X1 = c(X1, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2)
                                        Y1 = c(Y1, 3, 4, 5)
                                        Y2 = c(Y2, 3, 4.5, 4.5)
                                }
                                
                                # appending here:
                                
                                if (spd >= 65 & spd < 70) {
                                        X1 = c(X1, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 1, 2, 2, 2)
                                        Y1 = c(Y1, 0, 2, 3, 4, 5)
                                        Y2 = c(Y2, 0, 2, 3, 4.5, 4.5)
                                }
                                
                                if (spd >= 70 & spd < 75) {
                                        X1 = c(X1, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 2, 2, 2, 2)
                                        Y1 = c(Y1, 0, 2, 3, 4, 5)
                                        Y2 = c(Y2, 0, 2, 3, 4.5, 4.5)
                                }
                                
                                if (spd >= 75 & spd < 80) {
                                        X1 = c(X1, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1, 2, 3, 4, 5)
                                        Y2 = c(Y2, 1, 2, 3, 4.5, 4.5)
                                }
                                
                                if (spd >= 80 & spd < 85) {
                                        X1 = c(X1, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 2.5, 3, 3.5, 4, 5)
                                        Y2 = c(Y2, 2.5, 3, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 85 & spd < 90) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 2, 2.5, 3, 3.5, 4, 5)
                                        Y2 = c(Y2, 2, 2.5, 3, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 90 & spd < 95) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 2, 2.5, 3, 3.5, 4, 5)
                                        Y2 = c(Y2, 2, 2.5, 3, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 95 & spd < 100) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 3.5, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 100 & spd < 105) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 0, 0, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 105 & spd < 110) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 0, 1, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 110 & spd < 115) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 0, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 115 & spd < 120) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 1, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 120 & spd < 125) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 125 & spd < 130) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 130 & spd < 135) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 135 & spd < 140) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 140 & spd < 145) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 145 & spd < 150) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 0.5, 1.0, 1.5, 2, 2.5, 3, 4, 4, 5)
                                        Y2 = c(Y2, 0.5, 1.0, 1.5, 2, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 150 & spd < 155) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 0, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 155 & spd < 160) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 1, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                
                                if (spd >= 160 & spd < 165) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 165 & spd < 170) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 0, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 170 & spd < 175) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 175 & spd < 180 ) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 180 & spd < 185 ) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 1, 2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 0.5, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 0.5, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                if (spd >= 185  ) {
                                        X1 = c(X1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                                        X2 = c(X2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
                                        Y1 = c(Y1, 0.5, 1.0, 1.5, 2, 3, 3, 4, 4, 5)
                                        Y2 = c(Y2, 0.5, 1.0, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 4.5)
                                }
                                
                                
                                # currently wind barbs are only supported up to 190 knots
                                if (spd > 190) {
                                        text = paste("Currently wind barbs are supported only up to 190 knots. Your wind speed is", spd)
                                        message(text)
                                        }
                                
                                
                                dir = (direction[i]/360) * 2 * pi
                                rot = cbind(c(cos(dir), -sin(dir)), c(sin(dir),
                                                                       cos(dir)))
                                S1 = rbind(X1, Y1)
                                S2 = rbind(X2, Y2)
                                S1 = rot %*% S1
                                S2 = rot %*% S2
                                S1 = S1 * c(scalex, scaley) + c(x, y)
                                S2 = S2 * c(scalex, scaley) + c(x, y)
                        }
                        
                        if (spd > 0) {
                                if (ns > 3) {
                                        if (spd > mspd + 3 * rspd | spd < mspd - 3 *
                                            rspd)
                                                spdcolor = "black"
                                }
                                
                                # fill the last barb:
                                if(spd >= 50){
                                        polygon(x = c(S1[1,ncol(S1)],
                                                      S2[1,(ncol(S2)-1)],
                                                      S1[1,(ncol(S1)-1)]),
                                                y = c(S1[2,ncol(S1)],
                                                      S2[2,(ncol(S2)-1)],
                                                      S1[2,(ncol(S1)-1)]),
                                                col = "black"
                                        )
                                }
                                
                                if(spd >= 100){
                                        polygon(x = c(S1[1,ncol(S1)-1],
                                                      S2[1,(ncol(S2)-2)],
                                                      S1[1,(ncol(S1)-3)]),
                                                y = c(S1[2,ncol(S1)-1],
                                                      S2[2,(ncol(S2)-2)],
                                                      S1[2,(ncol(S1)-3)]),
                                                col = "black"
                                        )
                                }
                                
                                if(spd >= 150){
                                        kol = ifelse(spd > 190, "red", "black")
                                        polygon(x = c(S1[1,ncol(S1)-3],
                                                      S2[1,(ncol(S2)-4)],
                                                      S1[1,(ncol(S1)-5)]),
                                                y = c(S1[2,ncol(S1)-3],
                                                      S2[2,(ncol(S2)-4)],
                                                      S1[2,(ncol(S1)-5)]),
                                                col = kol
                                        )
                                }
                                
                                segments(S1[1, ], S1[2, ], S2[1, ], S2[2, ],
                                         col = spdcolor, lwd = 1)
                        }
                }
        }
}

