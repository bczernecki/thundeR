#' Sounding export
#' 
#' Internal package function for exporting interpolated profile with 5 m steps
#' 
#' 
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m) altitude [meters]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' 
#' @importFrom climate sounding_wyoming
#' @export
#' 
#' @return Data frame of computed values for visualizing parcel trajectories
#'  \enumerate{
#'   \item pressure pressure [hPa]
#'   \item altitude altitude [m]
#'   \item temp temperature [degree Celsius]
#'   \item tempV virtual temperature [degree Celsius]
#'   \item dpt dew point temperature [degree Celsius]
#'   \item wd wind direction [azimuth in degrees]
#'   \item ws wind speed [knots]
#'   \item MU temperature for most unstable CAPE trajectory [degree Celsius]
#'   \item SB temperature for surface based CAPE trajectory [degree Celsius]
#'   \item ML temperature for mixed layer CAPE trajectory [degree Celsius]
#' }
#' 
#' @examples 
#' data("sounding_vienna")
#' attach(sounding_vienna)
#' skewt_plot(close_par = FALSE)
#' output = sounding_export(pressure, altitude, temp, dpt, wd, ws)
#' skewt_lines(output$dpt, output$pressure, col = 'forestgreen', lwd = 2.5)
#' skewt_lines(output$temp,output$pressure, col = 'red', lwd = 2.5)
#' skewt_lines(output$MU,output$pressure, col = "orange", lty = 1, lwd = 2)
#' skewt_lines(output$tempV,output$pressure, col = "red3", lty = 3, lwd = 1.5)

sounding_export = function(pressure, altitude, temp, dpt, wd, ws){

  parametry = sounding_default(pressure = pressure, altitude = altitude, 
                               temp = temp, dpt = dpt, wd = wd, ws = ws, 
                               export_profile = 1, accuracy = 3)
  
  LP = length(sounding_default(pressure = pressure, altitude = altitude, 
                               temp = temp, dpt = dpt, wd = wd, ws = ws, 
                               export_profile = 0, accuracy = 1)) # no. of parameters 
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
DNw <- parametry[(LP+pozMU+(8*pozALL)+12):(LP+pozMU+(9*pozALL)+11)]
TVw <- parametry[(LP+pozMU+(9*pozALL)+13):(LP+pozMU+(10*pozALL)+12)]
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
                  ML = MLw,
                  DN = DNw)

res$DN[res$altitude-res$altitude[1] >= 4000] <- NA

return(res)
  
}
