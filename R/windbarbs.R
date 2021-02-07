#' Wind barbs
#' 
#' A slightly modified version of station.symbol function from the RadioSonde package.
#' @import stats
#' @import utils
#'
#' @export
#' @param cx x coordinates on a plot
#' @param cy y coordinates on a plot
#' @param direction wind direction (0-360 degrees)
#' @param speed wind speed in m/s
#' @param cex symbole size. Default 1


windbarbs <- function (cx, cy, direction, speed = NA, cex = 1) {
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
        tpar <- par()
        size <- tpar$csi
        scalex <- (tpar$usr[2] - tpar$usr[1])/(tpar$pin[1])
        scaley <- (tpar$usr[4] - tpar$usr[3])/(tpar$pin[2])
        scalex <- (cex * (scalex * size))/4.5
        scaley <- (cex * (scaley * size))/4.5
        for (i in 1:ns) {
                spdcolor = "black"
                pscolor = "blue"
                x = cx[i]
                y = cy[i]
                if (is.na(x) | is.na(y))
                        next
                spd = speed[i]
                
                if (!is.na(spd)) {
                        xs <- if (spd > 0) {
                                X1 <- 0
                                X2 <- 0
                                Y1 <- 0
                                Y2 <- 5
                                if (spd >= 5 & spd < 10) {
                                        X1 <- c(X1, 0)
                                        X2 <- c(X2, 1)
                                        Y1 <- c(Y1, 5)
                                        Y2 <- c(Y2, 5)
                                }
                                if (spd >= 10 & spd < 15) {
                                        X1 <- c(X1, 0)
                                        X2 <- c(X2, 2)
                                        Y1 <- c(Y1, 5)
                                        Y2 <- c(Y2, 5)
                                }
                                if (spd >= 15 & spd < 20) {
                                        X1 <- c(X1, 0, 0)
                                        X2 <- c(X2, 1, 2)
                                        Y1 <- c(Y1, 4, 5)
                                        Y2 <- c(Y2, 4, 5)
                                }
                                if (spd >= 20 & spd < 25) {
                                        X1 <- c(X1, 0, 0)
                                        X2 <- c(X2, 2, 2)
                                        Y1 <- c(Y1, 4, 5)
                                        Y2 <- c(Y2, 4, 5)
                                }
                                if (spd >= 25 & spd < 30) {
                                        X1 <- c(X1, 0, 0, 0)
                                        X2 <- c(X2, 1, 2, 2)
                                        Y1 <- c(Y1, 3, 4, 5)
                                        Y2 <- c(Y2, 3, 4, 5)
                                }
                                if (spd >= 30 & spd < 35) {
                                        X1 <- c(X1, 0, 0, 0)
                                        X2 <- c(X2, 2, 2, 2)
                                        Y1 <- c(Y1, 3, 4, 5)
                                        Y2 <- c(Y2, 3, 4, 5)
                                }
                                if (spd >= 35 & spd < 40) {
                                        X1 <- c(X1, 0, 0, 0, 0)
                                        X2 <- c(X2, 1, 2, 2, 2)
                                        Y1 <- c(Y1, 2, 3, 4, 5)
                                        Y2 <- c(Y2, 2, 3, 4, 5)
                                }
                                if (spd >= 40 & spd < 45) {
                                        X1 <- c(X1, 0, 0, 0, 0)
                                        X2 <- c(X2, 2, 2, 2, 2)
                                        Y1 <- c(Y1, 2, 3, 4, 5)
                                        Y2 <- c(Y2, 2, 3, 4, 5)
                                }
                                if (spd >= 45 & spd < 50) {
                                        X1 <- c(X1, 0, 0, 0, 0, 0)
                                        X2 <- c(X2, 1, 2, 2, 2, 2)
                                        Y1 <- c(Y1, 1, 2, 3, 4, 5)
                                        Y2 <- c(Y2, 1, 2, 3, 4, 5)
                                }
                                if (spd >= 50 & spd < 55) {
                                        X1 <- c(X1, 0, 0)
                                        X2 <- c(X2, 2, 2)
                                        Y1 <- c(Y1, 4, 5)
                                        Y2 <- c(Y2, 4.5, 4.5)
                                }
                                if (spd >= 55 & spd < 60) {
                                        X1 <- c(X1, 0, 0, 0)
                                        X2 <- c(X2, 1, 2, 2)
                                        Y1 <- c(Y1, 3, 4, 5)
                                        Y2 <- c(Y2, 3, 4.5, 4.5)
                                }
                                if (spd >= 60 & spd < 65) {
                                        X1 <- c(X1, 0, 0, 0)
                                        X2 <- c(X2, 2, 2, 2)
                                        Y1 <- c(Y1, 3, 4, 5)
                                        Y2 <- c(Y2, 3, 4.5, 4.5)
                                }
                                dir <- (direction[i]/360) * 2 * pi
                                rot <- cbind(c(cos(dir), -sin(dir)), c(sin(dir),
                                                                       cos(dir)))
                                S1 <- rbind(X1, Y1)
                                S2 <- rbind(X2, Y2)
                                S1 <- rot %*% S1
                                S2 <- rot %*% S2
                                S1 <- S1 * c(scalex, scaley) + c(x, y)
                                S2 <- S2 * c(scalex, scaley) + c(x, y)
                        }
                        if (spd > 0) {
                                if (ns > 3) {
                                        if (spd > mspd + 3 * rspd | spd < mspd - 3 *
                                            rspd)
                                                spdcolor = "black"
                                }
                                segments(S1[1, ], S1[2, ], S2[1, ], S2[2, ],
                                         col = spdcolor, lwd = 1.2)
                        }
                }
        }
}

