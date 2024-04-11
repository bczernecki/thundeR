#' Calculate convective parameters
#'
#' A core function for calculating convective parameters commonly used in the operational
#' prediction of severe convective storms. Returns a vector of parameters. 
#' 
#' \enumerate{
#'  \item MU_CAPE
#'  \item MU_CAPE_M10
#'  \item MU_CAPE_M10_PT
#'  \item MU_CAPE_3km
#'  \item MU_CAPE_HGL
#'  \item MU_buoyancy
#'  \item MU_buoyancy_M10
#'  \item MU_CIN
#'  \item MU_LCL_HGT
#'  \item MU_LFC_HGT
#'  \item MU_EL_HGT
#'  \item MU_LI
#'  \item MU_LI_M10
#'  \item MU_WMAX
#'  \item MU_EL_TEMP
#'  \item MU_LCL_TEMP
#'  \item MU_LFC_TEMP
#'  \item MU_MIXR
#'  \item MU_cold_cloud
#'  \item MU_warm_cloud
#'  \item MU_equal_layer
#'  \item MU_etilde
#'  \item MU_updraft_radius
#'  \item MU_ECAPE
#'  \item MU_ECAPE_HGL
#'  \item MU_ECAPE_M10
#'  \item MU_EWMAX
#'  \item MU_ECAPE_3km
#'  \item MUML_CAPE
#'  \item MUML_CAPE_M10
#'  \item MUML_CAPE_M10_PT
#'  \item MUML_CAPE_3km
#'  \item MUML_CAPE_HGL
#'  \item MUML_buoyancy
#'  \item MUML_buoyancy_M10
#'  \item MUML_CIN
#'  \item MUML_LCL_HGT
#'  \item MUML_LFC_HGT
#'  \item MUML_EL_HGT
#'  \item MUML_LI
#'  \item MUML_LI_M10
#'  \item MUML_WMAX
#'  \item MUML_EL_TEMP
#'  \item MUML_LCL_TEMP
#'  \item MUML_LFC_TEMP
#'  \item MUML_MIXR
#'  \item MUML_cold_cloud
#'  \item MUML_warm_cloud
#'  \item MUML_equal_layer
#'  \item MUML_etilde
#'  \item MUML_updraft_radius
#'  \item MUML_ECAPE
#'  \item MUML_ECAPE_HGL
#'  \item MUML_ECAPE_M10
#'  \item MUML_EWMAX
#'  \item MUML_ECAPE_3km
#'  \item SB_CAPE
#'  \item SB_CAPE_M10
#'  \item SB_CAPE_M10_PT
#'  \item SB_CAPE_3km
#'  \item SB_CAPE_HGL
#'  \item SB_buoyancy
#'  \item SB_buoyancy_M10
#'  \item SB_CIN
#'  \item SB_LCL_HGT
#'  \item SB_LFC_HGT
#'  \item SB_EL_HGT
#'  \item SB_LI
#'  \item SB_LI_M10
#'  \item SB_WMAX
#'  \item SB_EL_TEMP
#'  \item SB_LCL_TEMP
#'  \item SB_LFC_TEMP
#'  \item SB_MIXR
#'  \item SB_cold_cloud
#'  \item SB_warm_cloud
#'  \item SB_equal_layer
#'  \item SB_etilde
#'  \item SB_updraft_radius
#'  \item SB_ECAPE
#'  \item SB_ECAPE_HGL
#'  \item SB_ECAPE_M10
#'  \item SB_EWMAX
#'  \item SB_ECAPE_3km
#'  \item ML_CAPE
#'  \item ML_CAPE_M10
#'  \item ML_CAPE_M10_PT
#'  \item ML_CAPE_3km
#'  \item ML_CAPE_HGL
#'  \item ML_buoyancy
#'  \item ML_buoyancy_M10
#'  \item ML_CIN
#'  \item ML_LCL_HGT
#'  \item ML_LFC_HGT
#'  \item ML_EL_HGT
#'  \item ML_LI
#'  \item ML_LI_M10
#'  \item ML_WMAX
#'  \item ML_EL_TEMP
#'  \item ML_LCL_TEMP
#'  \item ML_LFC_TEMP
#'  \item ML_MIXR
#'  \item ML_cold_cloud
#'  \item ML_warm_cloud
#'  \item ML_equal_layer
#'  \item ML_etilde
#'  \item ML_updraft_radius
#'  \item ML_ECAPE
#'  \item ML_ECAPE_HGL
#'  \item ML_ECAPE_M10
#'  \item ML_EWMAX
#'  \item ML_ECAPE_3km
#'  \item MU500_CAPE
#'  \item MU500_CAPE_M10
#'  \item MU500_CAPE_M10_PT
#'  \item MU500_CAPE_3km
#'  \item MU500_CAPE_HGL
#'  \item MU500_CIN
#'  \item MU500_LI
#'  \item MU500_LI_M10
#'  \item MU500_buoyancy
#'  \item MU500_buoyancy_M10
#'  \item MU500_etilde
#'  \item MU500_radius
#'  \item MU500_ECAPE
#'  \item MU500_ECAPE_HGL
#'  \item MU500_ECAPE_M10
#'  \item MU500_EWMAX
#'  \item MU500_ECAPE_3km
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
#'  \item RH_01km 
#'  \item RH_02km 
#'  \item RH_14km 
#'  \item RH_25km 
#'  \item RH_36km 
#'  \item RH_HGL 
#'  \item FRZG_HGT 
#'  \item M10_HGT
#'  \item FRZG_wetbulb_HGT 
#'  \item MU_HGT
#'  \item MUML_HGT
#'  \item HGT_min_thetae_04km 
#'  \item Delta_thetae 
#'  \item Delta_thetae_min04km 
#'  \item Thetae_01km 
#'  \item Thetae_02km 
#'  \item DCAPE 
#'  \item Cold_Pool_Strength 
#'  \item PRCP_WATER 
#'  \item Moisture_Flux
#'  \item Moisture_Flux_SR
#'  \item Moisture_Flux_SR_eff
#'  \item BS_0500m
#'  \item BS_01km 
#'  \item BS_02km 
#'  \item BS_03km 
#'  \item BS_06km 
#'  \item BS_08km 
#'  \item BS_36km 
#'  \item BS_13km 
#'  \item BS_16km 
#'  \item BS_18km 
#'  \item BS_14km 
#'  \item BS_25km 
#'  \item BS_eff_MU 
#'  \item BS_eff_MUML
#'  \item BS_eff_SB 
#'  \item BS_eff_ML 
#'  \item BS_sfc_to_M10 
#'  \item BS_1km_to_M10 
#'  \item BS_2km_to_M10 
#'  \item BS_MU_LCL_to_M10 
#'  \item BS_MUML_LCL_to_M10 
#'  \item BS_SB_LCL_to_M10 
#'  \item BS_ML_LCL_to_M10
#'  \item MW_SR_0500m_RM
#'  \item MW_SR_01km_RM
#'  \item MW_SR_02km_RM
#'  \item MW_SR_03km_RM
#'  \item MW_SR_HGL_RM
#'  \item MW_SR_0500m_RM_eff
#'  \item MW_SRVM_0500m_RM
#'  \item MW_SRVM_01km_RM
#'  \item MW_SRVM_03km_RM
#'  \item MW_SR_0500m_LM
#'  \item MW_SR_01km_LM
#'  \item MW_SR_02km_LM
#'  \item MW_SR_03km_LM
#'  \item MW_SR_HGL_LM
#'  \item MW_SR_0500m_LM_eff
#'  \item MW_SRVM_0500m_LM
#'  \item MW_SRVM_01km_LM
#'  \item MW_SRVM_03km_LM
#'  \item MW_SR_0500m_MW
#'  \item MW_SR_01km_MW
#'  \item MW_SR_02km_MW
#'  \item MW_SR_03km_MW
#'  \item MW_SR_HGL_MW
#'  \item MW_SR_0500m_MW_eff
#'  \item Peters_SR_inflow
#'  \item Peters_SR_inflow_eff
#'  \item MW_0500m
#'  \item MW_01km 
#'  \item MW_02km 
#'  \item MW_03km 
#'  \item MW_06km 
#'  \item MW_13km 
#'  \item SRH_0100m_RM 
#'  \item SRH_0250m_RM 
#'  \item SRH_0500m_RM 
#'  \item SRH_01km_RM 
#'  \item SRH_03km_RM 
#'  \item SRH_36km_RM 
#'  \item SRH_0100m_LM 
#'  \item SRH_0250m_LM 
#'  \item SRH_0500m_LM 
#'  \item SRH_01km_LM 
#'  \item SRH_03km_LM 
#'  \item SRH_36km_LM
#'  \item SRH_01km_MW
#'  \item SRH_03km_MW
#'  \item SRH_36km_MW
#'  \item SRH_01km_RM_eff
#'  \item SRH_01km_LM_eff
#'  \item SRH_01km_MW_eff
#'  \item SRH_03km_RM_eff
#'  \item SRH_03km_LM_eff
#'  \item SRH_03km_MW_eff
#'  \item SV_0500m_RM
#'  \item SV_01km_RM
#'  \item SV_03km_RM
#'  \item SV_0500m_LM
#'  \item SV_01km_LM
#'  \item SV_03km_LM
#'  \item SV_FRA_0500m_RM
#'  \item SV_FRA_01km_RM
#'  \item SV_FRA_03km_RM
#'  \item SV_FRA_0500m_LM
#'  \item SV_FRA_01km_LM
#'  \item SV_FRA_03km_LM
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
#'  \item Peters_vector_A
#'  \item Peters_vector_M
#'  \item Peters_vector_eff_A
#'  \item Peters_vector_eff_M
#'  \item K_Index 
#'  \item Showalter_Index 
#'  \item TotalTotals_Index 
#'  \item SWEAT_Index 
#'  \item Wind_Index 
#'  \item STP_fix_RM
#'  \item STP_new_RM 
#'  \item STP_fix_LM 
#'  \item STP_new_LM 
#'  \item SCP_fix_RM 
#'  \item SCP_new_RM 
#'  \item SCP_fix_LM 
#'  \item SCP_new_LM 
#'  \item SHIP 
#'  \item HSI 
#'  \item HSIv2 
#'  \item DCP 
#'  \item MU_WMAXSHEAR 
#'  \item MUML_WMAXSHEAR
#'  \item SB_WMAXSHEAR 
#'  \item ML_WMAXSHEAR 
#'  \item MU_EFF_EWMAXSHEAR 
#'  \item MUML_EFF_EWMAXSHEAR
#'  \item SB_EFF_EWMAXSHEAR 
#'  \item ML_EFF_EWMAXSHEAR 
#'  \item MU_EFF_EWMAXSHEAR_HGL 
#'  \item MUML_EFF_EWMAXSHEAR_HGL
#'  \item SB_EFF_EWMAXSHEAR_HGL
#'  \item ML_EFF_EWMAXSHEAR_HGL
#'  \item EHI_500m_RM 
#'  \item EHI_01km_RM 
#'  \item EHI_03km_RM
#'  \item EHI_500m_LM 
#'  \item EHI_01km_LM
#'  \item EHI_03km_LM
#'  \item SHERBS3
#'  \item SHERBE
#'  \item SHERBS3_v2
#'  \item SHERBE_v2
#'  \item DEI
#'  \item DEI_eff
#'  \item TIP
#'  \item THTE_LR03
#'  \item THTE_LR04
#'  \item THTE_LR13
#'  \item THTE_LR14
#'  \item THTE_LR5_eff
#'  \item THTE_LR4_eff

#' }
#'
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m)  altitude [metres]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' @param accuracy accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
#' @param interpolate_step interpolation step to be used for vertical interpolation. Valid only if `accuracy` is set to 3 (default is 5 m)
#' @param meanlayer_bottom_top (optional) vector of length 2 for bottom and top heights used for computing parcel starting parameters; default: 0, 500
#' @param storm_motion (optional) for moving storms only - one can define vector of length two with
#' wind speed (m/s) and wind directions (degrees) that will be used to compute adjusted SRH parameters
#' @return Named vector of 200+ convective indices
#' @export 
#' @examples
#' old_options = options(scipen = 99) 
#' pressure = c(1000, 855, 700, 500, 300, 100, 10)
#' altitude = c(0, 1500, 2500, 6000, 8500, 12000, 25000)
#' temp = c(25, 10, 0, -15, -30, -50, -92)
#' dpt = c(20, 5, -5, -30, -55, -80, -99)
#' wd = c(0, 90, 135, 180, 270, 350, 0)
#' ws = c(5, 10, 20, 30, 40, 5, 0)
#' accuracy = 2
#' sounding_compute(pressure, altitude, temp, dpt, wd, ws, accuracy)
#' options(old_options) 

sounding_compute = function(pressure, altitude, temp, dpt, wd, ws,
                            accuracy = 2,
                            interpolate_step = 5,
                            meanlayer_bottom_top = c(0, 500),
                            storm_motion = c(999, 999)) {
  
  export_profile = 0 
  
  if (sum(storm_motion == c(999, 999)) != 2) {
    ws2 = storm_motion[1]
    wd2 = storm_motion[2]
    storm_motion[2] = round(ws2 * sin(wd2 * pi/180), 2)
    storm_motion[1] = round(ws2 * cos(wd2 * pi/180), 2)
    storm_motion[3] = 0
  } else {
    storm_motion[3] = 999
  }

  tmp = thunder::sounding_default(pressure, altitude, temp, dpt, wd, ws,
                                  export_profile, accuracy, interpolate_step,
                                  meanlayer_bottom_top, storm_motion)
  
  names(tmp) = c(
"MU_CAPE",
"MU_CAPE_M10",
"MU_CAPE_M10_PT",
"MU_CAPE_3km",
"MU_CAPE_HGL",
"MU_buoyancy",
"MU_buoyancy_M10",
"MU_CIN",
"MU_LCL_HGT",
"MU_LFC_HGT",
"MU_EL_HGT",
"MU_LI",
"MU_LI_M10",
"MU_WMAX",
"MU_EL_TEMP",
"MU_LCL_TEMP",
"MU_LFC_TEMP",
"MU_MIXR",
"MU_cold_cloud",
"MU_warm_cloud",
"MU_equal_layer",
"MU_etilde",
"MU_updraft_radius",
"MU_ECAPE",
"MU_ECAPE_HGL",
"MU_ECAPE_M10",
"MU_EWMAX",
"MU_ECAPE_3km",
"MUML_CAPE",
"MUML_CAPE_M10",
"MUML_CAPE_M10_PT",
"MUML_CAPE_3km",
"MUML_CAPE_HGL",
"MUML_buoyancy",
"MUML_buoyancy_M10",
"MUML_CIN",
"MUML_LCL_HGT",
"MUML_LFC_HGT",
"MUML_EL_HGT",
"MUML_LI",
"MUML_LI_M10",
"MUML_WMAX",
"MUML_EL_TEMP",
"MUML_LCL_TEMP",
"MUML_LFC_TEMP",
"MUML_MIXR",
"MUML_cold_cloud",
"MUML_warm_cloud",
"MUML_equal_layer",
"MUML_etilde",
"MUML_updraft_radius",
"MUML_ECAPE",
"MUML_ECAPE_HGL",
"MUML_ECAPE_M10",
"MUML_EWMAX",
"MUML_ECAPE_3km",
"SB_CAPE",
"SB_CAPE_M10",
"SB_CAPE_M10_PT",
"SB_CAPE_3km",
"SB_CAPE_HGL",
"SB_buoyancy",
"SB_buoyancy_M10",
"SB_CIN",
"SB_LCL_HGT",
"SB_LFC_HGT",
"SB_EL_HGT",
"SB_LI",
"SB_LI_M10",
"SB_WMAX",
"SB_EL_TEMP",
"SB_LCL_TEMP",
"SB_LFC_TEMP",
"SB_MIXR",
"SB_cold_cloud",
"SB_warm_cloud",
"SB_equal_layer",
"SB_etilde",
"SB_updraft_radius",
"SB_ECAPE",
"SB_ECAPE_HGL",
"SB_ECAPE_M10",
"SB_EWMAX",
"SB_ECAPE_3km",
"ML_CAPE",
"ML_CAPE_M10",
"ML_CAPE_M10_PT",
"ML_CAPE_3km",
"ML_CAPE_HGL",
"ML_buoyancy",
"ML_buoyancy_M10",
"ML_CIN",
"ML_LCL_HGT",
"ML_LFC_HGT",
"ML_EL_HGT",
"ML_LI",
"ML_LI_M10",
"ML_WMAX",
"ML_EL_TEMP",
"ML_LCL_TEMP",
"ML_LFC_TEMP",
"ML_MIXR",
"ML_cold_cloud",
"ML_warm_cloud",
"ML_equal_layer",
"ML_etilde",
"ML_updraft_radius",
"ML_ECAPE",
"ML_ECAPE_HGL",
"ML_ECAPE_M10",
"ML_EWMAX",
"ML_ECAPE_3km",
"MU500_CAPE",
"MU500_CAPE_M10",
"MU500_CAPE_M10_PT",
"MU500_CAPE_3km",
"MU500_CAPE_HGL",
"MU500_CIN",
"MU500_LI",
"MU500_LI_M10",
"MU500_buoyancy",
"MU500_buoyancy_M10",
"MU500_etilde",
"MU500_radius",
"MU500_ECAPE",
"MU500_ECAPE_HGL",
"MU500_ECAPE_M10",
"MU500_EWMAX",
"MU500_ECAPE_3km",
"LR_0500m",
"LR_01km",
"LR_02km",
"LR_03km",
"LR_04km",
"LR_06km",
"LR_16km",
"LR_26km",
"LR_24km",
"LR_36km",
"LR_26km_MAX",
"LR_500700hPa",
"LR_500800hPa",
"LR_600800hPa",
"RH_01km",
"RH_02km",
"RH_14km",
"RH_25km",
"RH_36km",
"RH_HGL",
"FRZG_HGT",
"M10_HGT",
"FRZG_wetbulb_HGT",
"MU_HGT",
"MUML_HGT",
"HGT_min_thetae_04km",
"Delta_thetae",
"Delta_thetae_min04km",
"Thetae_01km",
"Thetae_02km",
"DCAPE",
"Cold_Pool_Strength",
"PRCP_WATER",
"Moisture_Flux",
"Moisture_Flux_SR",
"Moisture_Flux_SR_eff",
"BS_0500m",
"BS_01km",
"BS_02km",
"BS_03km",
"BS_06km",
"BS_08km",
"BS_36km",
"BS_13km",
"BS_16km",
"BS_18km",
"BS_14km",
"BS_25km",
"BS_eff_MU",
"BS_eff_MUML",
"BS_eff_SB",
"BS_eff_ML",
"BS_sfc_to_M10",
"BS_1km_to_M10",
"BS_2km_to_M10",
"BS_MU_LCL_to_M10",
"BS_MUML_LCL_to_M10",
"BS_SB_LCL_to_M10",
"BS_ML_LCL_to_M10",
"MW_SR_0500m_RM",
"MW_SR_01km_RM",
"MW_SR_02km_RM",
"MW_SR_03km_RM",
"MW_SR_HGL_RM",
"MW_SR_0500m_RM_eff",
"MW_SRVM_0500m_RM",
"MW_SRVM_01km_RM",
"MW_SRVM_03km_RM",
"MW_SR_0500m_LM",
"MW_SR_01km_LM",
"MW_SR_02km_LM",
"MW_SR_03km_LM",
"MW_SR_HGL_LM",
"MW_SR_0500m_LM_eff",
"MW_SRVM_0500m_LM",
"MW_SRVM_01km_LM",
"MW_SRVM_03km_LM",
"MW_SR_0500m_MW",
"MW_SR_01km_MW",
"MW_SR_02km_MW",
"MW_SR_03km_MW",
"MW_SR_HGL_MW",
"MW_SR_0500m_MW_eff",
"Peters_SR_inflow",
"Peters_SR_inflow_eff",
"MW_0500m",
"MW_01km",
"MW_02km",
"MW_03km",
"MW_06km",
"MW_13km",
"SRH_0100m_RM",
"SRH_0250m_RM",
"SRH_0500m_RM",
"SRH_01km_RM",
"SRH_03km_RM",
"SRH_36km_RM",
"SRH_0100m_LM",
"SRH_0250m_LM",
"SRH_0500m_LM",
"SRH_01km_LM",
"SRH_03km_LM",
"SRH_36km_LM",
"SRH_01km_MW",
"SRH_03km_MW",
"SRH_36km_MW",
"SRH_01km_RM_eff",
"SRH_01km_LM_eff",
"SRH_01km_MW_eff",
"SRH_03km_RM_eff",
"SRH_03km_LM_eff",
"SRH_03km_MW_eff",
"SV_0500m_RM",
"SV_01km_RM",
"SV_03km_RM",
"SV_0500m_LM",
"SV_01km_LM",
"SV_03km_LM",
"SV_FRA_0500m_RM",
"SV_FRA_01km_RM",
"SV_FRA_03km_RM",
"SV_FRA_0500m_LM",
"SV_FRA_01km_LM",
"SV_FRA_03km_LM",
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
"Peters_vector_A",
"Peters_vector_M",
"Peters_vector_eff_A",
"Peters_vector_eff_M",
"K_Index",
"Showalter_Index",
"TotalTotals_Index",
"SWEAT_Index",
"Wind_Index",
"STP_fix_RM",
"STP_new_RM",
"STP_fix_LM",
"STP_new_LM",
"SCP_fix_RM",
"SCP_new_RM",
"SCP_fix_LM",
"SCP_new_LM",
"SHIP",
"HSI",
"HSIv2",
"DCP",
"MU_WMAXSHEAR",
"MUML_WMAXSHEAR",
"SB_WMAXSHEAR",
"ML_WMAXSHEAR",
"MU_EFF_EWMAXSHEAR",
"MUML_EFF_EWMAXSHEAR",
"SB_EFF_EWMAXSHEAR",
"ML_EFF_EWMAXSHEAR",
"MU_EFF_EWMAXSHEAR_HGL",
"MUML_EFF_EWMAXSHEAR_HGL",
"SB_EFF_EWMAXSHEAR_HGL",
"ML_EFF_EWMAXSHEAR_HGL",
"EHI_500m_RM",
"EHI_01km_RM",
"EHI_03km_RM",
"EHI_500m_LM",
"EHI_01km_LM",
"EHI_03km_LM",
"SHERBS3",
"SHERBE",
"SHERBS3_v2",
"SHERBE_v2",
"DEI",
"DEI_eff",
"TIP",
"THTE_LR03",
"THTE_LR04",
"THTE_LR13",
"THTE_LR14",
"THTE_LR5_eff",
"THTE_LR4_eff")
  
  return(tmp)
}
