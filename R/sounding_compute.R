#' Calculate convective parameters
#'
#' A core function for calculating convective parameters commonly used in the operational prediction of severe convective storms. Returns a vector of parameters. 
#'
#' @export
#' 
#' @return Named vector of 100+ convective indices
#' \enumerate{
#'  \item MU_CAPE 
#'  \item MU_CAPE_M10_fraction 
#'  \item MU_02km_CAPE 
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
#'  \item MU_CAPE_500
#'  \item MU_CIN_500
#'  \item MU_LI_500
#'  \item SB_CAPE 
#'  \item SB_CAPE_M10_fraction 
#'  \item SB_02km_CAPE 
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
#'  \item ML_CAPE_M10_fraction 
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
#'  \item LR_0500m 
#'  \item LR_01km 
#'  \item LR_02km 
#'  \item LR_03km 
#'  \item LR_04km  
#'  \item LR_06km 
#'  \item LR_16km 
#'  \item LR_26km
#'  \item LR_24km 
#'  \item LR_36km 
#'  \item LR_26km_MAX 
#'  \item LR_500700hPa 
#'  \item LR_500800hPa 
#'  \item LR_600800hPa 
#'  \item FRZG_HGT 
#'  \item FRZG_wetbulb_HGT 
#'  \item HGT_max_thetae_03km 
#'  \item HGT_min_thetae_04km 
#'  \item Delta_thetae 
#'  \item Delta_thetae_min04km 
#'  \item DCAPE 
#'  \item Cold_Pool_Strength 
#'  \item Wind_Index 
#'  \item PRCP_WATER 
#'  \item Moisture_Flux_02km 
#'  \item RH_02km 
#'  \item RH_14km 
#'  \item RH_25km 
#'  \item RH_36km 
#'  \item RH_HGL 
#'  \item BS_01km 
#'  \item BS_02km 
#'  \item BS_03km 
#'  \item BS_06km 
#'  \item BS_08km 
#'  \item BS_36km 
#'  \item BS_26km 
#'  \item BS_16km 
#'  \item BS_18km 
#'  \item BS_EFF_MU 
#'  \item BS_EFF_SB 
#'  \item BS_EFF_ML 
#'  \item BS_SFC_to_M10 
#'  \item BS_1km_to_M10 
#'  \item BS_2km_to_M10 
#'  \item BS_MU_LFC_to_M10 
#'  \item BS_SB_LFC_to_M10 
#'  \item BS_ML_LFC_to_M10
#'  \item BS_MW02_to_SM 
#'  \item BS_MW02_to_RM 
#'  \item BS_MW02_to_LM 
#'  \item BS_HGL_to_SM 
#'  \item BS_HGL_to_RM 
#'  \item BS_HGL_to_LM 
#'  \item MW_01km 
#'  \item MW_02km 
#'  \item MW_03km 
#'  \item MW_06km 
#'  \item MW_13km 
#'  \item SRH_100m_RM 
#'  \item SRH_250m_RM 
#'  \item SRH_500m_RM 
#'  \item SRH_1km_RM 
#'  \item SRH_3km_RM 
#'  \item SRH_36km_RM 
#'  \item SRH_100m_LM 
#'  \item SRH_250m_LM 
#'  \item SRH_500m_LM 
#'  \item SRH_1km_LM 
#'  \item SRH_3km_LM 
#'  \item SRH_36km_LM 
#'  \item Bunkers_RM_A 
#'  \item Bunkers_RM_M 
#'  \item Bunkers_LM_A 
#'  \item Bunkers_LM_M 
#'  \item Bunkers_MW_A 
#'  \item Bunkers_MW_M 
#'  \item Corfidi_downwind_A 
#'  \item Corfidi_downwind_M 
#'  \item Corfidi_upwind_A 
#'  \item Corfidi_upwind_M 
#'  \item K_Index 
#'  \item Showalter_Index 
#'  \item TotalTotals_Index 
#'  \item SWEAT_Index 
#'  \item STP_fix 
#'  \item STP_new 
#'  \item SCP_fix 
#'  \item SCP_new 
#'  \item SHIP 
#'  \item HSI 
#'  \item DCP 
#'  \item MU_WMAXSHEAR 
#'  \item SB_WMAXSHEAR 
#'  \item ML_WMAXSHEAR 
#'  \item MU_EFF_WMAXSHEAR 
#'  \item SB_EFF_WMAXSHEAR 
#'  \item ML_EFF_WMAXSHEAR 
#'  \item EHI_500m 
#'  \item EHI_01km 
#'  \item EHI_03km

#' }
#'
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m)  altitude [metres]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' @param accuracy accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
#' @export 
#' @examples
#' pressure = c(1000, 855, 700, 500, 300, 100, 10)
#' altitude = c(0, 1500, 2500, 6000, 8500, 12000, 25000)
#' temp = c(25, 10, 0, -15, -30, -50, -92)
#' dpt = c(20, 5, -5, -30, -55, -80, -99)
#' wd = c(0, 90, 135, 180, 270, 350, 0)
#' ws = c(5, 10, 20, 30, 40, 5, 0)
#' accuracy = 2 
#' options(digits = 2) #change output formatting precision 
#' sounding_compute(pressure, altitude, temp, dpt, wd, ws, accuracy)

sounding_compute = function(pressure, altitude, temp, dpt, wd, ws, accuracy = 2){
  
  export_profile = 0 
  
  tmp = sounding_default(pressure, altitude, temp, dpt, wd, ws, export_profile, accuracy)
  
  names(tmp) = c(
"MU_CAPE",
"MU_CAPE_M10_fraction",
"MU_02km_CAPE",
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
"MU_CAPE_500", 
"MU_CIN_500", 
"MU_LI_500", 
"SB_CAPE",
"SB_CAPE_M10_fraction",
"SB_02km_CAPE",
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
"ML_CAPE_M10_fraction",
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
"LR_0500m",
"LR_01km",
"LR_02km",
"LR_03km",
"LR_04km ",
"LR_06km",
"LR_16km",
"LR_26km", 
"LR_24km",
"LR_36km",
"LR_26km_MAX",
"LR_500700hPa",
"LR_500800hPa",
"LR_600800hPa",
"FRZG_HGT",
"FRZG_wetbulb_HGT",
"HGT_max_thetae_03km",
"HGT_min_thetae_04km",
"Delta_thetae",
"Delta_thetae_min04km",
"DCAPE",
"Cold_Pool_Strength",
"Wind_Index",
"PRCP_WATER",
"Moisture_Flux_02km",
"RH_02km",
"RH_14km",
"RH_25km",
"RH_36km",
"RH_HGL",
"BS_01km",
"BS_02km",
"BS_03km",
"BS_06km",
"BS_08km",
"BS_36km",
"BS_26km",
"BS_16km",
"BS_18km",
"BS_EFF_MU",
"BS_EFF_SB",
"BS_EFF_ML",
"BS_SFC_to_M10",
"BS_1km_to_M10",
"BS_2km_to_M10",
"BS_MU_LFC_to_M10",
"BS_SB_LFC_to_M10",
"BS_ML_LFC_to_M10", 
"BS_MW02_to_SM",
"BS_MW02_to_RM",
"BS_MW02_to_LM",
"BS_HGL_to_SM",
"BS_HGL_to_RM",
"BS_HGL_to_LM",
"MW_01km",
"MW_02km",
"MW_03km",
"MW_06km",
"MW_13km",
"SRH_100m_RM",
"SRH_250m_RM",
"SRH_500m_RM",
"SRH_1km_RM",
"SRH_3km_RM",
"SRH_36km_RM",
"SRH_100m_LM",
"SRH_250m_LM",
"SRH_500m_LM",
"SRH_1km_LM",
"SRH_3km_LM",
"SRH_36km_LM",
"Bunkers_RM_A",
"Bunkers_RM_M",
"Bunkers_LM_A",
"Bunkers_LM_M",
"Bunkers_MW_A",
"Bunkers_MW_M",
"Corfidi_downwind_A",
"Corfidi_downwind_M",
"Corfidi_upwind_A",
"Corfidi_upwind_M",
"K_Index",
"Showalter_Index",
"TotalTotals_Index",
"SWEAT_Index",
"STP_fix",
"STP_new",
"SCP_fix",
"SCP_new",
"SHIP",
"HSI",
"DCP",
"MU_WMAXSHEAR",
"SB_WMAXSHEAR",
"ML_WMAXSHEAR",
"MU_EFF_WMAXSHEAR",
"SB_EFF_WMAXSHEAR",
"ML_EFF_WMAXSHEAR",
"EHI_500m",
"EHI_01km",
"EHI_03km")
       
  return(tmp)
}
