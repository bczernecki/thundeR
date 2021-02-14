#' Calculate convective indices
#'
#' Function calculates presently over 100 atmospheric thermal-, kinematic- and composite indices 
#' (e.g. CAPE, CIN, Lapse Rates, EHI, STP, SCP, etc...)
#' related mostly to severe weather phenomena (Thunderstorms, Large Hail, Tornadoes, etc...)
#' 
#' Wrapper for a generic C++ based function. Please consider using a slightly faster 
#' `?sounding_default()` for computing large databases (i.e. at least over 1,000,000 soundings).
#'
#' @export
#' 
#' @return 
#' \enumerate{
#'  \item MU_CAPE
#'  \item MU_03km_CAPE
#'  \item MU_HGL_CAPE
#'  \item MU_CIN
#'  \item MU_LCL_HGT
#'  \item MU_LFC_HGT
#'  \item MU_EL_HGT
#'  \item MU_LI
#'  \item MU_WMAX
#'  \item MU_EL_TEMP
#'  \item MU_LCL_TEMP
#'  \item MU_LFC_TEMP
#'  \item MU_MIXR
#'  \item SB_CAPE
#'  \item SB_03km_CAPE
#'  \item SB_HGL_CAPE
#'  \item SB_CIN
#'  \item SB_LCL_HGT
#'  \item SB_LFC_HGT
#'  \item SB_EL_HGT
#'  \item SB_LI
#'  \item SB_WMAX
#'  \item SB_EL_TEMP
#'  \item SB_LCL_TEMP
#'  \item SB_LFC_TEMP
#'  \item SB_MIXR
#'  \item ML_CAPE
#'  \item ML_03km_CAPE
#'  \item ML_HGL_CAPE
#'  \item ML_CIN
#'  \item ML_LCL_HGT
#'  \item ML_LFC_HGT
#'  \item ML_EL_HGT
#'  \item ML_LI
#'  \item ML_WMAX
#'  \item ML_EL_TEMP
#'  \item ML_LCL_TEMP
#'  \item ML_LFC_TEMP
#'  \item ML_MIXR
#'  \item LR_01km
#'  \item LR_03km
#'  \item LR_24km
#'  \item LR_36km
#'  \item LR_500700hPa
#'  \item LR_500800hPa
#'  \item FRZG_HGT
#'  \item FRZG_wetbulb_HGT
#'  \item HGT_max_thetae_03km
#'  \item HGT_min_thetae_04km
#'  \item Delta_thetae
#'  \item DCAPE
#'  \item Cold_Pool_Strength
#'  \item Wind_Index  
#'  \item PRCP_WATER
#'  \item Moisture_Flux_02km
#'  \item RH_02km
#'  \item RH_25km
#'  \item RH_HGL
#'  \item BS_01km
#'  \item BS_02km
#'  \item BS_03km
#'  \item BS_06km
#'  \item BS_08km
#'  \item BS_36km
#'  \item BS_18km
#'  \item BS_EFF_MU
#'  \item BS_EFF_SB
#'  \item BS_EFF_ML
#'  \item BS_SFC_to_HGL
#'  \item BS_MU_LFC_to_HGL
#'  \item BS_SB_LFC_to_HGL
#'  \item BS_ML_LFC_to_HGL
#'  \item MW_01km
#'  \item MW_02km
#'  \item MW_06km
#'  \item MW_13km
#'  \item SRH_100m_RM
#'  \item SRH_500m_RM
#'  \item SRH_1km_RM
#'  \item SRH_3km_RM
#'  \item SRH_100m_LM
#'  \item SRH_500m_LM
#'  \item SRH_1km_LM
#'  \item SRH_3km_LM
#'  \item K_Index
#'  \item Showalter_Index
#'  \item TotalTotals_Index
#'  \item SWEAT_Index
#'  \item STP
#'  \item STP_new
#'  \item SCP
#'  \item SCP_new
#'  \item SHIP
#'  \item DCP
#'  \item MU_WMAXSHEAR
#'  \item SB_WMAXSHEAR
#'  \item ML_WMAXSHEAR
#'  \item MU_EFF_WMAXSHEAR
#'  \item SB_EFF_WMAXSHEAR
#'  \item ML_EFF_WMAXSHEAR
#' }
#'
#' @param pressure pressure [hPa]
#' @param altitude altitude [metres]
#' @param temp air temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [degrees]
#' @param ws wind speed [knots or m/s - TODO/TODISCUSS]
#' @param export_profile  whether to export interpolated levels for drawing on Skew-T diagram. Binary [0 - default, no extra data for drawing exported, 1 - extracting data to be used ]
#' @param accuracy how accurately integrate the positive/negative area. Valid options (1 - default, 2 - fast implementation, 3 - very accurate)
#' @export 
#' @examples
#' pressure <- c(1000, 855, 700, 500, 300, 100, 10)
#' altitude <- c(0, 1500, 2500, 6000, 8500, 12000, 25000)
#' temp <- c(25, 10, 0, -15, -30, -50, -92)
#' dpt <- c(20, 5, -5, -30, -55, -80, -99)
#' wd <- c(0, 90, 135, 180, 270, 350, 0)
#' ws <- c(5, 10, 20, 30, 40, 5, 0)
#' options(scipen = 999) # change formatting
#' sounding_compute(pressure, altitude, temp, dpt, wd, ws)


sounding_compute = function(pressure, altitude, temp, dpt, wd, ws, 
                            export_profile = 0, accuracy = 1){
  
  tmp = sounding_default(pressure, altitude, temp, dpt, wd, ws, export_profile, accuracy)
  
  names(tmp) = c(
  "MU_CAPE",
  "MU_03km_CAPE",
  "MU_HGL_CAPE",
  "MU_CIN",
  "MU_LCL_HGT",
  "MU_LFC_HGT",
  "MU_EL_HGT",
  "MU_LI",
  "MU_WMAX",
  "MU_EL_TEMP",
  "MU_LCL_TEMP",
  "MU_LFC_TEMP",
  "MU_MIXR",
  "SB_CAPE",
  "SB_03km_CAPE",
  "SB_HGL_CAPE",
  "SB_CIN",
  "SB_LCL_HGT",
  "SB_LFC_HGT",
  "SB_EL_HGT",
  "SB_LI",
  "SB_WMAX",
  "SB_EL_TEMP",
  "SB_LCL_TEMP",
  "SB_LFC_TEMP",
  "SB_MIXR",
  "ML_CAPE",
  "ML_03km_CAPE",
  "ML_HGL_CAPE",
  "ML_CIN",
  "ML_LCL_HGT",
  "ML_LFC_HGT",
  "ML_EL_HGT",
  "ML_LI",
  "ML_WMAX",
  "ML_EL_TEMP",
  "ML_LCL_TEMP",
  "ML_LFC_TEMP",
  "ML_MIXR",
  "LR_01km",
  "LR_03km",
  "LR_24km",
  "LR_36km",
  "LR_500700hPa",
  "LR_500800hPa",
  "FRZG_HGT",
  "FRZG_wetbulb_HGT",
  "HGT_max_thetae_03km",
  "HGT_min_thetae_04km",
  "Delta_thetae",
  "DCAPE",
  "Cold_Pool_Strength",
  "Wind_Index  ",
  "PRCP_WATER",
  "Moisture_Flux_02km",
  "RH_02km",
  "RH_25km",
  "RH_HGL",
  "BS_01km",
  "BS_02km",
  "BS_03km",
  "BS_06km",
  "BS_08km",
  "BS_36km",
  "BS_18km",
  "BS_EFF_MU",
  "BS_EFF_SB",
  "BS_EFF_ML",
  "BS_SFC_to_HGL",
  "BS_MU_LFC_to_HGL",
  "BS_SB_LFC_to_HGL",
  "BS_ML_LFC_to_HGL",
  "MW_01km",
  "MW_02km",
  "MW_06km",
  "MW_13km",
  "SRH_100m_RM",
  "SRH_500m_RM",
  "SRH_1km_RM",
  "SRH_3km_RM",
  "SRH_100m_LM",
  "SRH_500m_LM",
  "SRH_1km_LM",
  "SRH_3km_LM",
  "K_Index",
  "Showalter_Index",
  "TotalTotals_Index",
  "SWEAT_Index",
  "STP",
  "STP_new",
  "SCP",
  "SCP_new",
  "SHIP",
  "DCP",
  "MU_WMAXSHEAR",
  "SB_WMAXSHEAR",
  "ML_WMAXSHEAR",
  "MU_EFF_WMAXSHEAR",
  "SB_EFF_WMAXSHEAR",
  "ML_EFF_WMAXSHEAR")
       
  return(tmp)
}
