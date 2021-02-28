#' Plot empty Skew-T diagram
#' 
#' Function for plotting a customized version of the Skew-T diagram. Please note that drawing Skew-T may require increasing size or modifying 
#' aspect ratio of plotting window in order to provide readable results.
#'
#' @import graphics
#' 
#' @param ptop Pressure top level to be used for plotting diagram. Valid options: 200, 150, 100 (default) and 50 hPa
#' @param temp_stripes logical, whether to draw color stripes for isotherms
#' @param isoterms_col color to be used for drawing dry isoterms
#' @param mixing_ratio_col color to be used for drawing mixing ratio isolines and labels. If set to NA or empty string isolines are not drawn
#' @param dry_adiabats_col color to be used for drawing dry adiabats. If set to NA or not provided drawing lines skipped
#' @param moist_adiabats_col color to be used for drawing moist adiabats. If set to NA or not provided drawing lines skipped
#' @param deg45 whether to preserve 45 degrees for diagonal isolines on Skew-T diagram regardless ploting window aspect ratio. [Logical, default: FALSE]
#' @param isotherm0 whether to deliminate 0 degree Celsius isother[Logical, default: TRUE]
#' @param ... additional (mostly graphical) parameters to be passed
#' @export
#' 
#' @examples 
#' skewt_plot(ptop = 100)
#' 
#' skewt_plot(ptop = 150, temp_stripes = FALSE) # na color stripes for temperature
#' 
#' skewt_plot(ptop = 100)
#' title("Your title")
#' mtext('WMO ID: 11035, 2011-08-23 1200 UTC', padj = -0.5, col = "white")
#' data("sounding_wien")
#' attach(sounding_wien)
#' # dataset obtained from Wyoming University's sounding webpage:
#' # climate::sounding_wyoming(wmo_id = 11035, yy = 2011, mm = 8, dd = 23, hh=12)[[1]]
#' 
#' output <- sounding_export(PRES, HGHT, TEMP, DWPT, DRCT, SKNT)
#' skewt_lines(output$dpt, output$pressure,type='l',col='forestgreen',lwd = 2.5)
#' skewt_lines(output$temp,output$pressure,type='l',col='red', lwd = 2.5)
#' skewt_lines(output$MU,output$pressure, col = "orange", lty = 1, lwd = 2)
#' skewt_lines(output$tempV,output$pressure, col = "red3", lty = 3, lwd = 1.5)
#' 

skewt_plot = function(ptop = 100, 
                      isoterms_col = "#d8be9b",
                      temp_stripes = FALSE, 
                      mixing_ratio_col = "#8470FF90", 
                      dry_adiabats_col = "#d6878750",
                      moist_adiabats_col = "#00FF0095", 
                      deg45 = FALSE,
                      isotherm0 = TRUE,
                          ...){
  
  if(deg45){
    par(pty = "s") # preserve correct aspect ratio
  }
  
  #par(pty = "s")
  # marginesy w ukladzie: 
  #c(bottom, left, top, right) 
  #par(mar = c(2, 1.5, 1 ,6))
  ymax <- skewty(1050)
  #ymin <- skewty(50)
  ymin <- skewty(ptop)
  xmin <- skewtx(-50, skewty(1050))
  # przesuwanie rysowanego zakresu na wykresie:
  xmax <- skewtx(48.3, skewty(995))
  
  xc <- c(xmin, xmin, xmax, xmax, xmin)
  yc <- c(ymin, ymax, ymax, ymin, ymin)
  
  plot(xc, yc, type = "l", axes = FALSE, xlab = "", ylab = "", lwd = 1) 
  
  # par("usr")
  # [1] -29.828975  23.667431  -2.424056  37.790678
  # abline(v = c(-29.82, 23.66))
  ypos <- skewty(1050)
  degc <- seq(-50, 50, by = 10)
  axis(1, at = skewtx(degc, ypos), labels = seq(-50, 50, by = 10), pos = ymax, cex.axis=0.65, padj=-0.15, tck=-0.01)
  mtext(side = 1, line = 0, expression(paste("Temperature [\u00b0C]")), cex=0.65)
  
  pres <- c(1050, 1000, 850, 700, 500, 300, 200, 100)
  NPRES <- length(pres)
  xpl <- rep(xmin, times = NPRES)
  xpr <- c(xmax, xmax, xmax, xmax, skewtx(20, skewty(500)))
  
  #abline(h = y)
  ypos <- skewty(pres[2:NPRES])
  axis(2, las = 1, at = ypos, labels = pres[2:NPRES], pos = xmin, cex.axis = 0.65, lwd=0)
  mtext(side = 2, line = 1.3, "Pressure [hPa]", padj = 2, cex = 0.65)
  # end of drawing diagram in square
  
  # dry isotherms
  kinkx <- skewtx(10.5, skewty(400))
  temp <- seq(from = -150, to = 60, by = 10)
  NTEMP <- length(temp[temp < 60])
  lendt <- rep(1050, NTEMP) # lower limit in hPa
  
  #if(ptop == 150) lendt[1:8] <- c(222, 303, 414, 565, 770, 1050, 1050, 1050)
  lendt[1:11] <- c(49, 63, 87, 118, 163, 222, 303, 414, 565, 770, 1050)
  
  inds <- seq(1, length(temp))[(temp > -50 & temp < 60)]
  exponent <- (127.182 - (kinkx - 0.54 * temp[inds])/0.90692)/44.061
  #rendt <- rep(50, NTEMP)
  
  rendt <- rep(ptop, NTEMP) # change for 100 or 150hpa
  rendt[inds] <- 10^exponent # upper limit
  
  # check limit of upper lines:
  rendt[which(rendt < ptop)] = ptop # clip coordinates up to ptop
  lendt[which(lendt < ptop)] = ptop
  yl <- skewty(rendt)
  xl <- skewtx(temp[temp < 60], yl)
  yr <- skewty(lendt)
  xr <- skewtx(temp[temp < 60], yr)
  
  segments(xl, yl, xr, yr, col = isoterms_col, lwd = 0.8)
  
  if(temp_stripes){
    strt <- ifelse(ptop == 150, 1, 2)
    for (i in seq(strt, length(xl), by = 2)){
      polygon(x = c(xl[i], xr[i], xr[i+1], xl[i+1]), 
              y = c(yl[i], yr[i], yr[i+1], yl[i+1]), 
              border = NA,
              col = "#f0e8f475")
    }
  }
  
  if(isotherm0){
    inds = which(temp == 0)
    segments(xl[inds], yl[inds], xr[inds], yr[inds], col = "blue3", lwd = 1, lty = 3)
    inds = which(temp == -20)
    segments(xl[inds], yl[inds], xr[inds], yr[inds], col = "blue3", lwd = 1, lty = 3)
  }
  # upper labels for temperature:
  #ind <- (round(xl - max(xl),2) != 0) & (xl > xmin)
  
  #text(xl, yl, temp)
  
  #text(xl[ind] - 0.3, yl[ind] + 1, 
  #     labels = paste(" ", as.character(temp[ind])), col = isoterms_col, cex = 0.9)
  # right labels for temperature:
  #ind <- (round(xl - max(xl),2) != 0) 
  #end of drawing isotherms
  
  # mixing ratio
  temp1050 <- c(256.65, 265.25, 274.45, 284.40, 295.1)-273.15
  temp600 <- c(250.15, 258.25, 266.90, 276.25, 286.25)-273.15
  
  yl <- skewty((rep(600, times = length(temp600))))
  xr <- skewtx(temp1050, skewty(rep(1050, times = length(temp1050))))
  yr <- skewty((rep(1050, times = length(temp1050))))
  xl <- skewtx(temp600, skewty(rep(600, times = length(temp600))))
  
  if(!is.na(mixing_ratio_col) || mixing_ratio_col != ""){
    segments(xl, yl, xr, yr, col = mixing_ratio_col, lwd = 0.8, lty = 1)
    text(xl, yl + 0.75, labels = c(1, 2, 4, 8, 16), col = "#8470FF90", adj = 0.5, cex = 0.6)
  }
  
  # dry adiabats (theta):
  theta <- seq(from = -50, to = 250, by = 10)
  NTHETA <- length(theta)
  lendth <- rep(ptop, times = NTHETA) # set limit to 50 or 100 hPa
  lendth[1:8] <- c(950, 690, 500, 360, 250, 170, 110, 60)
  
  # correction for ptop
  if((!is.na(dry_adiabats_col) || dry_adiabats_col != "")){
    lendth[lendth < ptop] <- ptop
    rendth <- rep(1050, times = NTHETA)
    for (itheta in 1:NTHETA) {
      p <- seq(from = lendth[itheta], to = rendth[itheta], length = 200)
      sy <- skewty(p)
      dry <- tda(theta[itheta], p)
      sx <- skewtx(dry, sy)
      sx <- sx[sx <= xmax]
      if(length(sx)){
        sy=sy[1:length(sx)]
        lines(sx, sy, lty = 1, col = dry_adiabats_col, lwd = 0.7)
      }
    }
  }
  # end of dry adiabats
  
  
  # beginning of moist adiabats:
  if(!is.na(moist_adiabats_col) || moist_adiabats_col != ""){
    p <- seq(from = 1050, to = ptop, by = -2)
    npts <- length(p)
    sy <- skewty(p)
    sx <- double(length = npts)
    
    pseudo <- c(34,28, 21, 14, 5, -5,-18,-32,-45)
    NPSEUDO <- length(pseudo)
    holdx <- matrix(0, nrow = npts, ncol = NPSEUDO)
    holdy <- matrix(0, nrow = npts, ncol = NPSEUDO)
    for (ipseudo in 1:NPSEUDO) {
      for (ilen in 1:npts) {
        moist <- satlft(pseudo[ipseudo], p[ilen])
        sx[ilen] <- skewtx(moist, sy[ilen])
      }
      inds <- (sx < xmin)
      sx[inds] <- NA
      sy[inds] <- NA
      holdx[, ipseudo] <- sx
      holdy[, ipseudo] <- sy
    }
    
    for (ipseudo in 1:NPSEUDO) {
      sx <- holdx[, ipseudo]
      sy <- holdy[, ipseudo]
      lines(sx, sy, lty = 3, col = moist_adiabats_col)
      #text(sx[length(sx)]+0.5, sy[length(sy)]-0.5, as.character(pseudo[ipseudo]), 
      #     col = moist_adiabats_col, adj = 0.5, cex = 0.75)
    }
  } # end of moist adiabats
  
  # pp <- recordPlot()
  # devtools::use_data(pp)
  # #saveRDS(object = pp, file = "data/pp.rds") # opcjonalne zapisanie do wykorzystania potem
  # rm(pp)
  # 
  # print(pp)
  
  # draw pressure heights:
  y <- skewty(pres)
  segments(-27.85, y, 26, y, col = "black", lwd = 0.25, lty = 1)
  
  
}

