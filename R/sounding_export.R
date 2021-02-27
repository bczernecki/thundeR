#' Sounding export
#' 
#' Internal function to be used for re-calculate parcel trajectories that can be analysed or drawn
#' on Skew-T diagrams
#' 
#' 
#' @param pressure - air pressure (hPa)
#' @param altitude - in metres
#' @param temp - air temperature (degree Celsius)
#' @param dpt - dew point temperature (degree Celsius)
#' @param wd - wind direction in degrees (0-360)
#' @param ws - wind speed in [m/s or knots / TODO / TODISCUSS]
#' 
#' @importFrom climate sounding_wyoming
#' @export
#' 
#' @examples 
#' data("sounding_wien")
#' attach(sounding_wien)
#' skewt_plot()
#' output <- sounding_export(pressure, altitude, temp, dpt, wd, ws)
#' skewt_lines(output$dpt, output$pressure, col = 'forestgreen',lwd = 2.5)
#' skewt_lines(output$temp,output$pressure, col = 'red', lwd = 2.5)
#' skewt_lines(output$MU,output$pressure, col = "orange", lty = 1, lwd = 2)
#' skewt_lines(output$tempV,output$pressure, col = "red3", lty = 3, lwd = 1.5)

sounding_export = function(pressure, altitude, temp, dpt, wd, ws){

  parametry = sounding_compute(pressure = pressure, altitude = altitude, 
                               temp = temp, dpt = dpt, wd = wd, ws = ws, 
                               export_profile = 1, accuracy = 3)
  
  LP = max(which(!is.na(names(parametry)))) # no. of parameters 
  ###
  pozMU = parametry[LP+1]
  MUs = parametry[LP+2]
  MUw = parametry[(LP+3):(LP+2+pozMU)]
  ###
  pozALL = parametry[LP+pozMU+3]
  SBw = parametry[(LP+pozMU+4):(LP+pozMU+pozALL+3)]
  MLw = parametry[(LP+pozMU+pozALL+5):(LP+pozMU+(2*pozALL)+4)]
  Pw = parametry[(LP+pozMU+(2*pozALL)+6):(LP+pozMU+(3*pozALL)+5)]
  Hw <- parametry[(LP+pozMU+(3*pozALL)+7):(LP+pozMU+(4*pozALL)+6)]
  Tw <- parametry[(LP+pozMU+(4*pozALL)+8):(LP+pozMU+(5*pozALL)+7)]
  TDw <- parametry[(LP+pozMU+(5*pozALL)+9):(LP+pozMU+(6*pozALL)+8)]
  WDw <- parametry[(LP+pozMU+(6*pozALL)+10):(LP+pozMU+(7*pozALL)+9)]
  WSw <- parametry[(LP+pozMU+(7*pozALL)+11):(LP+pozMU+(8*pozALL)+10)]
  TVw <- parametry[(LP+pozMU+(8*pozALL)+12):(LP+pozMU+(9*pozALL)+11)]
  MUw <- c(rep(NA,pozALL-pozMU),MUw) # correction for elevated MU
  ###
  res <- data.frame(pressure = Pw,
                    altitude = Hw,
                    temp = Tw, 
                    tempV = TVw, 
                    dpt = TDw,
                    wd = WDw,
                    ws = WSw,
                    MU = MUw,
                    SB = SBw,
                    ML = MLw)
  
  return(res)
  
}
