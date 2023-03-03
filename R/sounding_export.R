#' Sounding export
#'
#' Internal package function for exporting interpolated profile with 5 m (or user-defined) steps
#'
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m) altitude [meters]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' @param accuracy accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
#' @param interpolate_step interpolation step to be used for vertical interpolation. Valid only if `accuracy` is set to 3 (default is 5 m)
#' @param meanlayer_bottom_top (optional) vector of length 2 for bottom and top heights used for computing parcel starting parameters; default: 0, 500
#' @param storm_motion (optional) for moving storms only - one can define 
#' wind speed and wind directions (TODO: units!!!) that will be used to compute adjusted SRH parameters
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
#' skewt_lines(output$dpt, output$pressure, col = "forestgreen", lwd = 2.5)
#' skewt_lines(output$temp, output$pressure, col = "red", lwd = 2.5)
#' skewt_lines(output$MU, output$pressure, col = "orange", lty = 1, lwd = 2)
#' skewt_lines(output$tempV, output$pressure, col = "red3", lty = 3, lwd = 1.5)
sounding_export = function(pressure, altitude, temp, dpt, wd, ws, accuracy = 3,
                           interpolate_step = 5,
                           meanlayer_bottom_top = c(0, 500),
                           storm_motion = c(999, 999)) {
  
  parametry = sounding_default(
    pressure = pressure, altitude = altitude,
    temp = temp, dpt = dpt, wd = wd, ws = ws,
    export_profile = 1, accuracy = accuracy, interpolate_step = interpolate_step,
    meanlayer_bottom_top = meanlayer_bottom_top, storm_motion = storm_motion
  )

  LP = length(sounding_default(
    pressure = pressure, altitude = altitude,
    temp = temp, dpt = dpt, wd = wd, ws = ws,
    export_profile = 0, accuracy = 1, interpolate_step = interpolate_step,
    meanlayer_bottom_top = meanlayer_bottom_top, storm_motion = storm_motion
  )) # no. of parameters

  ###
  pozMU = parametry[LP + 1]
  MUs = parametry[LP + 2]
  MUw = parametry[(LP + 3):(LP + 2 + pozMU)]
  ###
  pozALL = parametry[LP + pozMU + 3]
  SBw = parametry[(LP + pozMU + 4):(LP + pozMU + pozALL + 3)]
  pozML = parametry[(LP + pozMU + pozALL + 4)]
  MLs = parametry[(LP + pozMU + pozALL + 5)]
  MLw = parametry[(LP + pozMU + pozALL + 6):(LP + pozMU + pozALL + pozML + 5)]
  Pw = parametry[(LP + pozMU + pozML + (pozALL) + 7):(LP + pozMU + pozML + (2 * pozALL) + 6)]
  Hw = parametry[(LP + pozMU + pozML + (2 * pozALL) + 8):(LP + pozMU + pozML + (3 * pozALL) + 7)]
  Tw = parametry[(LP + pozMU + pozML + (3 * pozALL) + 9):(LP + pozMU + pozML + (4 * pozALL) + 8)]
  TDw = parametry[(LP + pozMU + pozML + (4 * pozALL) + 10):(LP + pozMU + pozML + (5 * pozALL) + 9)]
  WDw = parametry[(LP + pozMU + pozML + (5 * pozALL) + 11):(LP + pozMU + pozML + (6 * pozALL) + 10)]
  WSw = parametry[(LP + pozMU + pozML + (6 * pozALL) + 12):(LP + pozMU + pozML + (7 * pozALL) + 11)]
  DNw = parametry[(LP + pozMU + pozML + (7 * pozALL) + 13):(LP + pozMU + pozML + (8 * pozALL) + 12)]
  TVw = parametry[(LP + pozMU + pozML + (8 * pozALL) + 14):(LP + pozMU + pozML + (9 * pozALL) + 13)]
  MUw = c(rep(NA, pozALL - pozMU), MUw) # correction for elevated MU
  MLw = c(rep(NA, pozALL - pozML), MLw) # correction for elevated ML

  ###
  res = data.frame(
    pressure = Pw,
    altitude = Hw,
    temp = Tw,
    tempV = TVw,
    dpt = TDw,
    wd = WDw,
    ws = WSw,
    MU = MUw,
    SB = SBw,
    ML = MLw,
    DN = DNw
  )

  res$DN[res$altitude - res$altitude[1] >= 4000] = NA

  return(res)
}
