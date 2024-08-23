#' Calculate convective parameters
#'
#' A core function for calculating convective parameters commonly used in the operational
#' prediction of severe convective storms. Returns a vector of parameters. 
#' 
#' \enumerate{
#'  \item 	SB_CAPE
#'  \item 	SB_CAPE_M10
#'  \item 	SB_CAPE_M10_PT
#'  \item 	SB_CAPE_3km
#'  \item 	SB_CAPE_HGL
#'  \item 	SB_buoy
#'  \item 	SB_buoy_M10
#'  \item 	SB_buoy_3km
#'  \item 	SB_LI
#'  \item 	SB_LI_M10
#'  \item 	SB_CIN
#'  \item 	SB_LCL_hgt
#'  \item 	SB_LFC_hgt
#'  \item 	SB_EL_hgt
#'  \item 	SB_LCL_tmp
#'  \item 	SB_LFC_tmp
#'  \item 	SB_EL_tmp
#'  \item 	SB_cold_cloud
#'  \item 	SB_warm_cloud
#'  \item 	SB_equal_layer
#'  \item 	SB_MIXR
#'  \item 	SB_WMAXSHEAR
#'  \item 	SB_E_CAPE
#'  \item 	SB_E_CAPE_M10
#'  \item 	SB_E_CAPE_3km
#'  \item 	SB_E_CAPE_HGL
#'  \item 	SB_E_buoy
#'  \item 	SB_E_buoy_M10
#'  \item 	SB_E_buoy_3km
#'  \item 	SB_E_LI
#'  \item 	SB_E_tilde 
#'  \item 	SB_E_radius
#'  \item 	SB_E_WMAXSHEAR
#'  \item 	SB_E_WMAXSHEAR_HGL
#'  \item 	SB_E_WMAXSHEAR_3km
#'  \item 	ML_CAPE
#'  \item 	ML_CAPE_M10
#'  \item 	ML_CAPE_M10_PT
#'  \item 	ML_CAPE_3km
#'  \item 	ML_CAPE_HGL
#'  \item 	ML_buoy
#'  \item 	ML_buoy_M10
#'  \item 	ML_buoy_3km
#'  \item 	ML_LI
#'  \item 	ML_LI_M10
#'  \item 	ML_CIN
#'  \item 	ML_LCL_hgt
#'  \item 	ML_LFC_hgt
#'  \item 	ML_EL_hgt
#'  \item 	ML_LCL_tmp
#'  \item 	ML_LFC_tmp
#'  \item 	ML_EL_tmp
#'  \item 	ML_cold_cloud
#'  \item 	ML_warm_cloud
#'  \item 	ML_equal_layer
#'  \item 	ML_MIXR
#'  \item 	ML_WMAXSHEAR
#'  \item 	ML_E_CAPE
#'  \item 	ML_E_CAPE_M10
#'  \item 	ML_E_CAPE_3km
#'  \item 	ML_E_CAPE_HGL
#'  \item 	ML_E_buoy
#'  \item 	ML_E_buoy_M10
#'  \item 	ML_E_buoy_3km
#'  \item 	ML_E_LI
#'  \item 	ML_E_tilde 
#'  \item 	ML_E_radius
#'  \item 	ML_E_WMAXSHEAR
#'  \item 	ML_E_WMAXSHEAR_HGL
#'  \item 	ML_E_WMAXSHEAR_3km
#'  \item 	MU_CAPE
#'  \item 	MU_CAPE_M10
#'  \item 	MU_CAPE_M10_PT
#'  \item 	MU_CAPE_3km
#'  \item 	MU_CAPE_HGL
#'  \item 	MU_buoy
#'  \item 	MU_buoy_M10
#'  \item 	MU_buoy_3km
#'  \item 	MU_LI
#'  \item 	MU_LI_M10
#'  \item 	MU_CIN
#'  \item 	MU_LCL_hgt
#'  \item 	MU_LFC_hgt
#'  \item 	MU_EL_hgt
#'  \item 	MU_LCL_tmp
#'  \item 	MU_LFC_tmp
#'  \item 	MU_EL_tmp
#'  \item 	MU_cold_cloud
#'  \item 	MU_warm_cloud
#'  \item 	MU_equal_layer
#'  \item 	MU_MIXR
#'  \item 	MU_WMAXSHEAR
#'  \item 	MU_E_CAPE
#'  \item 	MU_E_CAPE_M10
#'  \item 	MU_E_CAPE_3km
#'  \item 	MU_E_CAPE_HGL
#'  \item 	MU_E_buoy
#'  \item 	MU_E_buoy_M10
#'  \item 	MU_E_buoy_3km
#'  \item 	MU_E_LI
#'  \item 	MU_E_tilde 
#'  \item 	MU_E_radius
#'  \item 	MU_E_WMAXSHEAR
#'  \item 	MU_E_WMAXSHEAR_HGL
#'  \item 	MU_E_WMAXSHEAR_3km
#'  \item 	MUML_CAPE
#'  \item 	MUML_CAPE_M10
#'  \item 	MUML_CAPE_M10_PT
#'  \item 	MUML_CAPE_3km
#'  \item 	MUML_CAPE_HGL
#'  \item 	MUML_buoy
#'  \item 	MUML_buoy_M10
#'  \item 	MUML_buoy_3km
#'  \item 	MUML_LI
#'  \item 	MUML_LI_M10
#'  \item 	MUML_CIN
#'  \item 	MUML_LCL_hgt
#'  \item 	MUML_LFC_hgt
#'  \item 	MUML_EL_hgt
#'  \item 	MUML_LCL_tmp
#'  \item 	MUML_LFC_tmp
#'  \item 	MUML_EL_tmp
#'  \item 	MUML_cold_cloud
#'  \item 	MUML_warm_cloud
#'  \item 	MUML_equal_layer
#'  \item 	MUML_MIXR
#'  \item 	MUML_WMAXSHEAR
#'  \item 	MUML_E_CAPE
#'  \item 	MUML_E_CAPE_M10
#'  \item 	MUML_E_CAPE_3km
#'  \item 	MUML_E_CAPE_HGL
#'  \item 	MUML_E_buoy
#'  \item 	MUML_E_buoy_M10
#'  \item 	MUML_E_buoy_3km
#'  \item 	MUML_E_LI
#'  \item 	MUML_E_tilde 
#'  \item 	MUML_E_radius
#'  \item 	MUML_E_WMAXSHEAR
#'  \item 	MUML_E_WMAXSHEAR_HGL
#'  \item 	MUML_E_WMAXSHEAR_3km
#'  \item 	MU5_CAPE
#'  \item 	MU5_CAPE_M10
#'  \item 	MU5_CAPE_M10_PT
#'  \item 	MU5_CAPE_3km
#'  \item 	MU5_CAPE_HGL
#'  \item 	MU5_buoy
#'  \item 	MU5_buoy_M10
#'  \item 	MU5_buoy_3km
#'  \item 	MU5_LI
#'  \item 	MU5_LI_M10
#'  \item 	MU5_CIN
#'  \item 	MU5_E_CAPE
#'  \item 	MU5_E_CAPE_M10
#'  \item 	MU5_E_CAPE_3km
#'  \item 	MU5_E_CAPE_HGL
#'  \item 	MU5_E_buoy
#'  \item 	MU5_E_buoy_M10
#'  \item 	MU5_E_buoy_3km
#'  \item 	MU5_E_LI
#'  \item 	MU5_E_tilde 
#'  \item 	MU5_E_radius
#'  \item 	LR_0500m
#'  \item 	LR_01km
#'  \item 	LR_02km
#'  \item 	LR_03km
#'  \item 	LR_04km
#'  \item 	LR_06km
#'  \item 	LR_16km
#'  \item 	LR_24km
#'  \item 	LR_26km
#'  \item 	LR_36km
#'  \item 	LR_26km_max
#'  \item 	LR_500700
#'  \item 	LR_500800
#'  \item 	LR_600800
#'  \item 	HGT_ISO_0
#'  \item 	HGT_ISO_0_wetbulb
#'  \item 	HGT_ISO_M10
#'  \item 	HGT_MU
#'  \item 	HGT_MUML
#'  \item 	HGT_THETAE_MIN
#'  \item 	THETAE_delta
#'  \item 	THETAE_delta_4km
#'  \item 	THETAE_LCL_M10
#'  \item 	THETAE_MU_M10
#'  \item 	THETAE_01km
#'  \item 	THETAE_02km
#'  \item 	THETAE_LR_03km
#'  \item 	THETAE_LR_13km
#'  \item 	THETAE_LR_14km
#'  \item 	DCAPE
#'  \item 	Cold_Pool_Strength
#'  \item 	RH_01km
#'  \item 	RH_02km
#'  \item 	RH_14km
#'  \item 	RH_25km
#'  \item 	RH_36km
#'  \item 	RH_HGL
#'  \item 	RH_500850
#'  \item 	PRCP_WATER
#'  \item 	PRCP_WATER_eff
#'  \item 	Moisture_Flux_02km
#'  \item 	Moisture_Flux_SR
#'  \item 	Moisture_Flux_SR_eff
#'  \item 	MW_0500m
#'  \item 	MW_01km
#'  \item 	MW_02km
#'  \item 	MW_03km
#'  \item 	MW_06km
#'  \item 	MW_13km
#'  \item 	BS_0500m
#'  \item 	BS_01km
#'  \item 	BS_02km
#'  \item 	BS_03km
#'  \item 	BS_06km
#'  \item 	BS_08km
#'  \item 	BS_36km
#'  \item 	BS_13km
#'  \item 	BS_14km
#'  \item 	BS_16km
#'  \item 	BS_18km
#'  \item 	BS_25km
#'  \item 	BS_eff_SB
#'  \item 	BS_eff_ML
#'  \item 	BS_eff_MU
#'  \item 	BS_eff_MUML
#'  \item 	BS_eff_M500
#'  \item 	BS_0km_M10
#'  \item 	BS_1km_M10
#'  \item 	BS_2km_M10
#'  \item 	BS_ML_LCL_M10
#'  \item 	BS_MU_LCL_M10
#'  \item 	BS_MUML_LCL_M10
#'  \item 	BS_06km_smoothness
#'  \item 	BS_16km_smoothness
#'  \item 	SRW_0500m_RM
#'  \item 	SRW_0500m_LM
#'  \item 	SRW_0500m_MW
#'  \item 	SRW_0500m_CBV
#'  \item 	SRW_01km_RM
#'  \item 	SRW_01km_LM
#'  \item 	SRW_01km_MW
#'  \item 	SRW_02km_RM
#'  \item 	SRW_02km_LM
#'  \item 	SRW_02km_MW
#'  \item 	SRW_03km_RM
#'  \item 	SRW_03km_LM
#'  \item 	SRW_03km_MW
#'  \item 	SRW_16km_RM
#'  \item 	SRW_16km_LM
#'  \item 	SRW_36km_RM
#'  \item 	SRW_36km_LM
#'  \item 	SRW_HGL_RM
#'  \item 	SRW_HGL_LM
#'  \item 	SRW_HGL_MW
#'  \item 	SRW_eff_RM
#'  \item 	SRW_eff_LM
#'  \item 	SRW_eff_MW
#'  \item 	SRW_eff_CBV
#'  \item 	SRW_VM_0500m_RM
#'  \item 	SRW_VM_0500m_LM
#'  \item 	SRW_VM_01km_RM
#'  \item 	SRW_VM_01km_LM
#'  \item 	SRW_VM_03km_RM
#'  \item 	SRW_VM_03km_LM
#'  \item 	SRW_VM_16km_RM
#'  \item 	SRW_VM_16km_LM
#'  \item 	SRW_VM_36km_RM
#'  \item 	SRW_VM_36km_LM
#'  \item 	Ventilation_13km
#'  \item 	Ventilation_16km
#'  \item 	Ventilation_36km
#'  \item 	Ventilation_69km
#'  \item 	Ventilation_912km
#'  \item 	SRH_0100m_RM
#'  \item 	SRH_0100m_LM
#'  \item 	SRH_0250m_RM
#'  \item 	SRH_0250m_LM
#'  \item 	SRH_0500m_RM
#'  \item 	SRH_0500m_LM
#'  \item 	SRH_01km_RM
#'  \item 	SRH_01km_LM
#'  \item 	SRH_01km_MW
#'  \item 	SRH_03km_RM
#'  \item 	SRH_03km_LM
#'  \item 	SRH_03km_MW
#'  \item 	SRH_36km_RM
#'  \item 	SRH_36km_LM
#'  \item 	SRH_36km_MW
#'  \item 	SRH_13km_RM
#'  \item 	SRH_13km_LM
#'  \item 	SRH_13km_MW
#'  \item 	SRH_16km_RM
#'  \item 	SRH_16km_LM
#'  \item 	SRH_16km_MW
#'  \item 	SRH_eff_1km_RM
#'  \item 	SRH_eff_1km_LM
#'  \item 	SRH_eff_1km_MW
#'  \item 	SRH_eff_3km_RM
#'  \item 	SRH_eff_3km_LM
#'  \item 	SRH_eff_3km_MW
#'  \item 	SV_0100m_RM
#'  \item 	SV_0100m_LM
#'  \item 	SV_0500m_RM
#'  \item 	SV_0500m_LM
#'  \item 	SV_01km_RM
#'  \item 	SV_01km_LM
#'  \item 	SV_03km_RM
#'  \item 	SV_03km_LM
#'  \item 	SV_0100m_RM_fra
#'  \item 	SV_0100m_LM_fra
#'  \item 	SV_0500m_RM_fra
#'  \item 	SV_0500m_LM_fra
#'  \item 	SV_01km_RM_fra
#'  \item 	SV_01km_LM_fra
#'  \item 	SV_03km_RM_fra
#'  \item 	SV_03km_LM_fra
#'  \item 	Bunkers_RM_A
#'  \item 	Bunkers_RM_M
#'  \item 	Bunkers_LM_A
#'  \item 	Bunkers_LM_M
#'  \item 	Bunkers_MW_A
#'  \item 	Bunkers_MW_M
#'  \item 	CBV_A
#'  \item 	CBV_M
#'  \item 	Corfidi_downwind_A
#'  \item 	Corfidi_downwind_M
#'  \item 	Corfidi_upwind_A
#'  \item 	Corfidi_upwind_M
#'  \item 	K_Index
#'  \item 	TotalTotals_Index
#'  \item 	SWEAT_Index
#'  \item 	Wind_Index
#'  \item 	STP_fix_RM
#'  \item 	STP_fix_LM
#'  \item 	STP_eff_RM
#'  \item 	STP_eff_LM
#'  \item 	SCP_fix_RM
#'  \item 	SCP_fix_LM
#'  \item 	SCP_eff_RM
#'  \item 	SCP_eff_LM
#'  \item 	SHIP
#'  \item 	HSI
#'  \item 	HSI_mod
#'  \item 	DCP
#'  \item 	DCP_eff
#'  \item 	EHI_0500m_RM
#'  \item 	EHI_0500m_LM
#'  \item 	EHI_01km_RM
#'  \item 	EHI_01km_LM
#'  \item 	EHI_03km_RM
#'  \item 	EHI_03km_LM
#'  \item 	SHERBS3
#'  \item 	SHERBE
#'  \item 	SHERBS3_mod
#'  \item 	SHERBE_mod
#'  \item 	DEI
#'  \item 	DEI_eff
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
#' @return Named vector of 300+ convective indices
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
"SB_CAPE",
"SB_CAPE_M10",
"SB_CAPE_M10_PT",
"SB_CAPE_3km",
"SB_CAPE_HGL",
"SB_buoy",
"SB_buoy_M10",
"SB_buoy_3km",
"SB_LI",
"SB_LI_M10",
"SB_CIN",
"SB_LCL_hgt",
"SB_LFC_hgt",
"SB_EL_hgt",
"SB_LCL_tmp",
"SB_LFC_tmp",
"SB_EL_tmp",
"SB_cold_cloud",
"SB_warm_cloud",
"SB_equal_layer",
"SB_MIXR",
"SB_WMAXSHEAR",
"SB_E_CAPE",
"SB_E_CAPE_M10",
"SB_E_CAPE_3km",
"SB_E_CAPE_HGL",
"SB_E_buoy",
"SB_E_buoy_M10",
"SB_E_buoy_3km",
"SB_E_LI",
"SB_E_tilde ",
"SB_E_radius",
"SB_E_WMAXSHEAR",
"SB_E_WMAXSHEAR_HGL",
"SB_E_WMAXSHEAR_3km",
"ML_CAPE",
"ML_CAPE_M10",
"ML_CAPE_M10_PT",
"ML_CAPE_3km",
"ML_CAPE_HGL",
"ML_buoy",
"ML_buoy_M10",
"ML_buoy_3km",
"ML_LI",
"ML_LI_M10",
"ML_CIN",
"ML_LCL_hgt",
"ML_LFC_hgt",
"ML_EL_hgt",
"ML_LCL_tmp",
"ML_LFC_tmp",
"ML_EL_tmp",
"ML_cold_cloud",
"ML_warm_cloud",
"ML_equal_layer",
"ML_MIXR",
"ML_WMAXSHEAR",
"ML_E_CAPE",
"ML_E_CAPE_M10",
"ML_E_CAPE_3km",
"ML_E_CAPE_HGL",
"ML_E_buoy",
"ML_E_buoy_M10",
"ML_E_buoy_3km",
"ML_E_LI",
"ML_E_tilde ",
"ML_E_radius",
"ML_E_WMAXSHEAR",
"ML_E_WMAXSHEAR_HGL",
"ML_E_WMAXSHEAR_3km",
"MU_CAPE",
"MU_CAPE_M10",
"MU_CAPE_M10_PT",
"MU_CAPE_3km",
"MU_CAPE_HGL",
"MU_buoy",
"MU_buoy_M10",
"MU_buoy_3km",
"MU_LI",
"MU_LI_M10",
"MU_CIN",
"MU_LCL_hgt",
"MU_LFC_hgt",
"MU_EL_hgt",
"MU_LCL_tmp",
"MU_LFC_tmp",
"MU_EL_tmp",
"MU_cold_cloud",
"MU_warm_cloud",
"MU_equal_layer",
"MU_MIXR",
"MU_WMAXSHEAR",
"MU_E_CAPE",
"MU_E_CAPE_M10",
"MU_E_CAPE_3km",
"MU_E_CAPE_HGL",
"MU_E_buoy",
"MU_E_buoy_M10",
"MU_E_buoy_3km",
"MU_E_LI",
"MU_E_tilde ",
"MU_E_radius",
"MU_E_WMAXSHEAR",
"MU_E_WMAXSHEAR_HGL",
"MU_E_WMAXSHEAR_3km",
"MUML_CAPE",
"MUML_CAPE_M10",
"MUML_CAPE_M10_PT",
"MUML_CAPE_3km",
"MUML_CAPE_HGL",
"MUML_buoy",
"MUML_buoy_M10",
"MUML_buoy_3km",
"MUML_LI",
"MUML_LI_M10",
"MUML_CIN",
"MUML_LCL_hgt",
"MUML_LFC_hgt",
"MUML_EL_hgt",
"MUML_LCL_tmp",
"MUML_LFC_tmp",
"MUML_EL_tmp",
"MUML_cold_cloud",
"MUML_warm_cloud",
"MUML_equal_layer",
"MUML_MIXR",
"MUML_WMAXSHEAR",
"MUML_E_CAPE",
"MUML_E_CAPE_M10",
"MUML_E_CAPE_3km",
"MUML_E_CAPE_HGL",
"MUML_E_buoy",
"MUML_E_buoy_M10",
"MUML_E_buoy_3km",
"MUML_E_LI",
"MUML_E_tilde ",
"MUML_E_radius",
"MUML_E_WMAXSHEAR",
"MUML_E_WMAXSHEAR_HGL",
"MUML_E_WMAXSHEAR_3km",
"MU5_CAPE",
"MU5_CAPE_M10",
"MU5_CAPE_M10_PT",
"MU5_CAPE_3km",
"MU5_CAPE_HGL",
"MU5_buoy",
"MU5_buoy_M10",
"MU5_buoy_3km",
"MU5_LI",
"MU5_LI_M10",
"MU5_CIN",
"MU5_E_CAPE",
"MU5_E_CAPE_M10",
"MU5_E_CAPE_3km",
"MU5_E_CAPE_HGL",
"MU5_E_buoy",
"MU5_E_buoy_M10",
"MU5_E_buoy_3km",
"MU5_E_LI",
"MU5_E_tilde ",
"MU5_E_radius",
"LR_0500m",
"LR_01km",
"LR_02km",
"LR_03km",
"LR_04km",
"LR_06km",
"LR_16km",
"LR_24km",
"LR_26km",
"LR_36km",
"LR_26km_max",
"LR_500700",
"LR_500800",
"LR_600800",
"HGT_ISO_0",
"HGT_ISO_0_wetbulb",
"HGT_ISO_M10",
"HGT_MU",
"HGT_MUML",
"HGT_THETAE_MIN",
"THETAE_delta",
"THETAE_delta_4km",
"THETAE_LCL_M10",
"THETAE_MU_M10",
"THETAE_01km",
"THETAE_02km",
"THETAE_LR_03km",
"THETAE_LR_13km",
"THETAE_LR_14km",
"DCAPE",
"Cold_Pool_Strength",
"RH_01km",
"RH_02km",
"RH_14km",
"RH_25km",
"RH_36km",
"RH_HGL",
"RH_500850",
"PRCP_WATER",
"PRCP_WATER_eff",
"Moisture_Flux_02km",
"Moisture_Flux_SR",
"Moisture_Flux_SR_eff",
"MW_0500m",
"MW_01km",
"MW_02km",
"MW_03km",
"MW_06km",
"MW_13km",
"BS_0500m",
"BS_01km",
"BS_02km",
"BS_03km",
"BS_06km",
"BS_08km",
"BS_36km",
"BS_13km",
"BS_14km",
"BS_16km",
"BS_18km",
"BS_25km",
"BS_eff_SB",
"BS_eff_ML",
"BS_eff_MU",
"BS_eff_MUML",
"BS_eff_M500",
"BS_0km_M10",
"BS_1km_M10",
"BS_2km_M10",
"BS_ML_LCL_M10",
"BS_MU_LCL_M10",
"BS_MUML_LCL_M10",
"BS_06km_smoothness",
"BS_16km_smoothness",
"SRW_0500m_RM",
"SRW_0500m_LM",
"SRW_0500m_MW",
"SRW_0500m_CBV",
"SRW_01km_RM",
"SRW_01km_LM",
"SRW_01km_MW",
"SRW_02km_RM",
"SRW_02km_LM",
"SRW_02km_MW",
"SRW_03km_RM",
"SRW_03km_LM",
"SRW_03km_MW",
"SRW_16km_RM",
"SRW_16km_LM",
"SRW_36km_RM",
"SRW_36km_LM",
"SRW_HGL_RM",
"SRW_HGL_LM",
"SRW_HGL_MW",
"SRW_eff_RM",
"SRW_eff_LM",
"SRW_eff_MW",
"SRW_eff_CBV",
"SRW_VM_0500m_RM",
"SRW_VM_0500m_LM",
"SRW_VM_01km_RM",
"SRW_VM_01km_LM",
"SRW_VM_03km_RM",
"SRW_VM_03km_LM",
"SRW_VM_16km_RM",
"SRW_VM_16km_LM",
"SRW_VM_36km_RM",
"SRW_VM_36km_LM",
"Ventilation_13km",
"Ventilation_16km",
"Ventilation_36km",
"Ventilation_69km",
"Ventilation_912km",
"SRH_0100m_RM",
"SRH_0100m_LM",
"SRH_0250m_RM",
"SRH_0250m_LM",
"SRH_0500m_RM",
"SRH_0500m_LM",
"SRH_01km_RM",
"SRH_01km_LM",
"SRH_01km_MW",
"SRH_03km_RM",
"SRH_03km_LM",
"SRH_03km_MW",
"SRH_36km_RM",
"SRH_36km_LM",
"SRH_36km_MW",
"SRH_13km_RM",
"SRH_13km_LM",
"SRH_13km_MW",
"SRH_16km_RM",
"SRH_16km_LM",
"SRH_16km_MW",
"SRH_eff_1km_RM",
"SRH_eff_1km_LM",
"SRH_eff_1km_MW",
"SRH_eff_3km_RM",
"SRH_eff_3km_LM",
"SRH_eff_3km_MW",
"SV_0100m_RM",
"SV_0100m_LM",
"SV_0500m_RM",
"SV_0500m_LM",
"SV_01km_RM",
"SV_01km_LM",
"SV_03km_RM",
"SV_03km_LM",
"SV_0100m_RM_fra",
"SV_0100m_LM_fra",
"SV_0500m_RM_fra",
"SV_0500m_LM_fra",
"SV_01km_RM_fra",
"SV_01km_LM_fra",
"SV_03km_RM_fra",
"SV_03km_LM_fra",
"Bunkers_RM_A",
"Bunkers_RM_M",
"Bunkers_LM_A",
"Bunkers_LM_M",
"Bunkers_MW_A",
"Bunkers_MW_M",
"CBV_A",
"CBV_M",
"Corfidi_downwind_A",
"Corfidi_downwind_M",
"Corfidi_upwind_A",
"Corfidi_upwind_M",
"K_Index",
"TotalTotals_Index",
"SWEAT_Index",
"Wind_Index",
"STP_fix_RM",
"STP_fix_LM",
"STP_eff_RM",
"STP_eff_LM",
"SCP_fix_RM",
"SCP_fix_LM",
"SCP_eff_RM",
"SCP_eff_LM",
"SHIP",
"HSI",
"HSI_mod",
"DCP",
"DCP_eff",
"EHI_0500m_RM",
"EHI_0500m_LM",
"EHI_01km_RM",
"EHI_01km_LM",
"EHI_03km_RM",
"EHI_03km_LM",
"SHERBS3",
"SHERBE",
"SHERBS3_mod",
"SHERBE_mod",
"DEI",
"DEI_eff")
return(tmp)
}
