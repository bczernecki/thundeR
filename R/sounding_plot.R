#' Plot Skew-T, hodograph and convective indices on a single layout 
#' 
#' Function to plot a composite of Skew-T, hodograph and selected convective parameters on a single layout
#' 
#' @param pressure pressure [hPa]
#' @param altitude altitude [m] (can be above sea level or above ground level as function always consider first level as surface, i.e h = 0 m) - altitude [meters]
#' @param temp temperature [degree Celsius]
#' @param dpt dew point temperature [degree Celsius]
#' @param wd wind direction [azimuth in degrees]
#' @param ws wind speed [knots]
#' @param title title to be added in the layout's header
#' @param parcel parcel tracing on Skew-T for "MU", "ML" or "SB" parcel
#' @param max_speed range of the hodograph to be drawn, 25 m/s used as default
#' @param hazards logical, whether to add extra information about possibility of convective hazards given convective initiation (default  = FALSE)
#' @param ... extra graphic arguments
#' @export
#' @import aiRthermo
#' @import grDevices
#' 
#' @return Skew-T, hodograph and table with convective indices drawn on a pre-defined single layout
#' 
#' @examples
#' data("sounding_vienna")
#' sounding_vienna = na.omit(sounding_vienna)
#' attach(sounding_vienna)
#' sounding_plot(pressure, altitude, temp, dpt, wd, ws, 
#'               parcel = "MU", title = "Vienna - 23 August 2011, 12:00 UTC")
#' 

sounding_plot = function(pressure, altitude, temp, dpt, wd, ws,
                        title = "", parcel = "MU", max_speed = 25, hazards = FALSE, ...){

   convert = FALSE
   ptop = 100 
  
   dev_size = dev.size("in")
  if(dev_size[1] < 10 | dev_size[2] < 7.5){
    text = paste("Your display device is", dev_size[1], "x", dev_size[2], "in. \nIt is recommended to use at least 10 x 7.5 in. plotting window \nor consider saving the layout into file")
    message(text)
  }

  # restore old par settings on exit:
  oldpar = par(no.readonly = TRUE) 
  on.exit(par(oldpar))

  # new par settings:
  par(pty = "m")
  plot.new()
  par(fig = c(0.028, 0.51, 0.03, 0.95), 
      new = T, 
      mar = c(0, 0, 0, 0), 
      oma = c(0, 0, 0, 0),
      mgp = c(0,0.15,0))

  ####
  t_col = function(color, percent, name = NULL) {
    rgb.val = col2rgb(color)
    t.col = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 maxColorValue = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    invisible(t.col)
  }

  ####
  output = sounding_export(pressure, altitude, temp, dpt, wd, ws)
  output2 = sounding_export(pressure, altitude, temp, ifelse(dpt==-273, NA,dpt), wd, ws)
  RH = aiRthermo::dewpointdepression2rh(output$pressure*100, output$temp+273.15, output$temp-output$dpt, consts = export_constants())
  SH = aiRthermo::rh2shum(output$pressure*100, output$temp+273.15, RH, consts = export_constants())
  E = aiRthermo::q2e(output$pressure*100, SH, consts = export_constants())
  W = aiRthermo::e2w(E, output$pressure*100, consts = export_constants())
  TLCL = aiRthermo::boltonTLCL(output$temp+273.15, RH, consts = export_constants())
  THETAE = aiRthermo::equivalentPotentialTemperature(output$pressure*100, output$temp+273.15, W, TLCL, consts = export_constants())
  ####

  skewt_plot(close_par = FALSE)
  skewt_lines(output2$dpt,output2$pressure, col = t_col('forestgreen',10), lwd = 2, ptop = 100)
  skewt_lines(output$temp,output$pressure, col = t_col('red',10), lwd = 2, ptop = 100)
  skewt_lines(output$tempV,output$pressure, col = t_col("red3",0), lwd = 1, lty = 3, ptop = 100)

  text(20, 28, "Hail Growth\nLayer (HGL)",col="#8470FF90",cex=0.65,srt=56)

  ####
  parametry = sounding_compute(pressure, altitude, temp, dpt, wd, ws, accuracy = 3)
  LP = max(which(!is.na(names(parametry))))

  ###

  if(hazards){
    rect(10.75,37.5,26.1,44,col=rgb(255,255,255, maxColorValue = 255, alpha = 200),lwd=0.2)
    text(11.5,43.25, "Possible storm hazards:",col="black",cex=0.65, adj=c(0,1))
    
    WIND = NA
    HAIL = NA
    TORN = NA
    
    ###
    if((parametry[which(names(parametry[1:LP]) == "DCAPE")]>700 |
        parametry[which(names(parametry[1:LP]) == "Delta_thetae")]>20 |
        parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")]>450 |
        parametry[which(names(parametry[1:LP]) == "MW_13km")]>15 | 
        parametry[which(names(parametry[1:LP]) == "SCP_new")]>1) & 
       parametry[which(names(parametry[1:LP]) == "MU_CAPE")]>150){
      WIND = 1
    }
    
    if((parametry[which(names(parametry[1:LP]) == "DCAPE")]>1100 |
        parametry[which(names(parametry[1:LP]) == "Delta_thetae")]>35 |
        parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")]>1000 |
        parametry[which(names(parametry[1:LP]) == "MW_13km")]>25 |
        parametry[which(names(parametry[1:LP]) == "SCP_new")]>10) & 
       parametry[which(names(parametry[1:LP]) == "MU_CAPE")]>150){
      WIND = 2
    }
    
    if(is.na(WIND)){WIND = 0}
    if(WIND ==1 ){text(11.5,41.75, "- 25m/s+ severe wind",  font=1, col="blue",cex=0.67, adj=c(0,1))}
    if(WIND ==2 ){text(11.5,41.75, "- 32m/s+ severe wind",  font=1, col="blue",cex=0.67, adj=c(0,1))}
    
    ###
    if(parametry[which(names(parametry[1:LP]) == "SHIP")]>0.5 |
       parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")]>450){
      HAIL = 1
    }
    
    if(parametry[which(names(parametry[1:LP]) == "SHIP")]>1 |
       parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")]>1000){
      HAIL = 2
    }
    
    if(is.na(HAIL)){HAIL = 0}
    if(HAIL ==1 ){text(11.5,40.25, "- 2cm+ large hail",  font=1, col="forestgreen",cex=0.67, adj=c(0,1))}
    if(HAIL ==2 ){text(11.5,40.25, "- 5cm+ large hail",  font=1, col="forestgreen",cex=0.67, adj=c(0,1))}
    
    ###
    if(parametry[which(names(parametry[1:LP]) == "STP")]>0.3 |
       parametry[which(names(parametry[1:LP]) == "STP_new")]>0.3 |
       parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")]>1000){
      TORN = 1
    }
    
    if(parametry[which(names(parametry[1:LP]) == "STP")]>0.75 |
       parametry[which(names(parametry[1:LP]) == "STP_new")]>0.75){
      TORN = 2
    }
    
    if(is.na(TORN)){TORN = 0}
    if(TORN ==1 ){text(11.5,38.75, "- F0+ tornado",  font=1, col="red",cex=0.67, adj=c(0,1))}
    if(TORN ==2 ){text(11.5,38.75, "- F2+ tornado",  font=1, col="red",cex=0.67, adj=c(0,1))}
    
  }
  
  ###

  if(parcel=="ML"){
    vsb_lcl = parametry[which(names(parametry[1:LP]) == "ML_LCL_HGT")] + output$altitude[1]
    vsb_lfc = parametry[which(names(parametry[1:LP]) == "ML_LFC_HGT")] + output$altitude[1]
    vsb_el = parametry[which(names(parametry[1:LP]) == "ML_EL_HGT")] + output$altitude[1]
    vsb_eff = (parametry[which(names(parametry[1:LP]) == "ML_EL_HGT")]/2) + output$altitude[1]
    ind_lcl = which.min(abs(output$altitude - vsb_lcl))
    ind_lfc = which.min(abs(output$altitude - vsb_lfc))
    ind_el = which.min(abs(output$altitude - vsb_el))
    ind_eff = which.min(abs(output$altitude - vsb_eff))
    y_eff = skewty(output$pressure[ind_eff])
    x_eff = skewtx(output$ML[ind_eff],skewty(output$pressure[ind_eff]))
    y_el = skewty(output$pressure[ind_el])
    x_el = skewtx(output$ML[ind_el],skewty(output$pressure[ind_el]))
    y_lfc = skewty(output$pressure[ind_lfc])
    x_lfc = skewtx(output$ML[ind_lfc],skewty(output$pressure[ind_lfc]))
    y_lcl = skewty(output$pressure[ind_lcl])
    x_lcl = skewtx(output$ML[ind_lcl],skewty(output$pressure[ind_lcl]))
    v = skewty(c(output$pressure[ind_lfc:ind_el], output$pressure[ind_el:ind_lfc])) # extra checks for NA coded as -99
    u = skewtx(c(output$tempV[ind_lfc:ind_el], rev(output$ML[ind_lfc:ind_el])), v)
    skewt_lines(output$ML,output$pressure, col = "orange", lty = 1, lwd = 1, ptop = 100)
  }

  if(parcel=="MU"){
    vsb_lcl = parametry[which(names(parametry[1:LP]) == "MU_LCL_HGT")] + output$altitude[1]
    vsb_lfc = parametry[which(names(parametry[1:LP]) == "MU_LFC_HGT")] + output$altitude[1]
    vsb_el = parametry[which(names(parametry[1:LP]) == "MU_EL_HGT")] + output$altitude[1]
    vsb_eff = ((parametry[which(names(parametry[1:LP]) == "MU_EL_HGT")]-parametry[which(names(parametry[1:LP]) == "HGT_max_thetae_03km")])/2) + output$altitude[1]
    ind_lcl = which.min(abs(output$altitude - vsb_lcl))
    ind_lfc = which.min(abs(output$altitude - vsb_lfc))
    ind_el = which.min(abs(output$altitude - vsb_el))
    ind_eff = which.min(abs(output$altitude - vsb_eff))
    y_eff = skewty(output$pressure[ind_eff])
    x_eff = skewtx(output$MU[ind_eff],skewty(output$pressure[ind_eff]))
    y_el = skewty(output$pressure[ind_el])
    x_el = skewtx(output$MU[ind_el],skewty(output$pressure[ind_el]))
    y_lfc = skewty(output$pressure[ind_lfc])
    x_lfc = skewtx(output$MU[ind_lfc],skewty(output$pressure[ind_lfc]))
    y_lcl = skewty(output$pressure[ind_lcl])
    x_lcl = skewtx(output$MU[ind_lcl],skewty(output$pressure[ind_lcl]))
    v = skewty(c(output$pressure[ind_lfc:ind_el], output$pressure[ind_el:ind_lfc])) # extra checks for NA coded as -99
    u = skewtx(c(output$tempV[ind_lfc:ind_el], rev(output$MU[ind_lfc:ind_el])), v)
    skewt_lines(output$MU,output$pressure, col = "orange", lty = 1, lwd = 1, ptop = 100)
  }

  if(parcel=="SB"){
    vsb_lcl = parametry[which(names(parametry[1:LP]) == "SB_LCL_HGT")] + output$altitude[1]
    vsb_lfc = parametry[which(names(parametry[1:LP]) == "SB_LFC_HGT")] + output$altitude[1]
    vsb_el = parametry[which(names(parametry[1:LP]) == "SB_EL_HGT")] + output$altitude[1]
    vsb_eff = (parametry[which(names(parametry[1:LP]) == "SB_EL_HGT")]/2) + output$altitude[1]
    ind_lcl = which.min(abs(output$altitude - vsb_lcl))
    ind_lfc = which.min(abs(output$altitude - vsb_lfc))
    ind_el = which.min(abs(output$altitude - vsb_el))
    ind_eff = which.min(abs(output$altitude - vsb_eff))
    y_eff = skewty(output$pressure[ind_eff])
    x_eff = skewtx(output$SB[ind_eff],skewty(output$pressure[ind_eff]))
    y_el = skewty(output$pressure[ind_el])
    x_el = skewtx(output$SB[ind_el],skewty(output$pressure[ind_el]))
    y_lfc = skewty(output$pressure[ind_lfc])
    x_lfc = skewtx(output$SB[ind_lfc],skewty(output$pressure[ind_lfc]))
    y_lcl = skewty(output$pressure[ind_lcl])
    x_lcl = skewtx(output$SB[ind_lcl],skewty(output$pressure[ind_lcl]))
    v = skewty(c(output$pressure[ind_lfc:ind_el], output$pressure[ind_el:ind_lfc])) # extra checks for NA coded as -99
    u = skewtx(c(output$tempV[ind_lfc:ind_el], rev(output$SB[ind_lfc:ind_el])), v)
    skewt_lines(output$SB,output$pressure, col = "orange", lty = 1, lwd = 1, ptop = 100)
  }

  ###
  u = subset(u, v < 44)
  v = subset(v, v < 44)
  polygon(u, v, col = "#FFA50025", border = NA)

  ###
  altitude_to_pressure = function(altitude){
    skewty(output$pressure[which(output$altitude-output$altitude[1] == altitude)])
  }
  
  ###  
  if(parcel=="SB"){
    if(parametry[which(names(parametry[1:LP]) == "SB_CAPE")] > 0){
      #text(x_eff, y_eff, paste0("---- Effective"), pos = 4, cex = 0.62, col = "black")
    if(output$pressure[which(output$altitude-output$altitude[1] == (parametry[which(names(parametry[1:LP]) == "SB_EL_HGT")]))] > 100 & which(names(parametry[1:LP]) == "SB_EL_HGT")!=0){
      text(x_el, y_el, paste0("---- SB EL"), pos = 4, cex = 0.62, col = "black") }
      text(x_lcl, y_lcl, paste0("---- SB LCL"), pos = 4, cex = 0.62, col = "black") }
  }

  if(parcel=="MU"){
    if(parametry[which(names(parametry[1:LP]) == "MU_CAPE")] > 0){
      #text(x_eff, y_eff, paste0("---- Effective"), pos = 4, cex = 0.62, col = "black")
      if(output$pressure[which(output$altitude-output$altitude[1] == (parametry[which(names(parametry[1:LP]) == "MU_EL_HGT")]))] > 100 & which(names(parametry[1:LP]) == "MU_EL_HGT")!=0){
        text(x_el, y_el, paste0("---- MU EL"), pos = 4, cex = 0.62, col = "black") }
      text(x_lcl, y_lcl, paste0("---- MU LCL"), pos = 4, cex = 0.62, col = "black") }
  }
  
  if(parcel=="ML"){
    if(parametry[which(names(parametry[1:LP]) == "ML_CAPE")] > 0){
      #text(x_eff, y_eff, paste0("---- Effective"), pos = 4, cex = 0.62, col = "black")
      if(output$pressure[which(output$altitude-output$altitude[1] == (parametry[which(names(parametry[1:LP]) == "ML_EL_HGT")]))] > 100 & which(names(parametry[1:LP]) == "ML_EL_HGT")!=0){
        text(x_el, y_el, paste0("---- ML EL"), pos = 4, cex = 0.62, col = "black") }
      text(x_lcl, y_lcl, paste0("---- ML LCL"), pos = 4, cex = 0.62, col = "black")}
  }

  ###
  text(-29, altitude_to_pressure(0), paste0("--- Sfc (",output$altitude[1]," m) ---"), pos=4, cex = 0.65, col = "black")
  if(max(output$altitude-output$altitude[1]) > 1000){if(output$pressure[output$altitude-output$altitude[1]==1000] > 100){text(-29, altitude_to_pressure(1000), paste0("--- 1 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 2000){if(output$pressure[output$altitude-output$altitude[1]==2000] > 100){text(-29, altitude_to_pressure(2000), paste0("--- 2 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 3000){if(output$pressure[output$altitude-output$altitude[1]==3000] > 100){text(-29, altitude_to_pressure(3000), paste0("--- 3 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 4000){if(output$pressure[output$altitude-output$altitude[1]==4000] > 100){text(-29, altitude_to_pressure(4000), paste0("--- 4 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 5000){if(output$pressure[output$altitude-output$altitude[1]==5000] > 100){text(-29, altitude_to_pressure(5000), paste0("--- 5 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 6000){if(output$pressure[output$altitude-output$altitude[1]==6000] > 100){text(-29, altitude_to_pressure(6000), paste0("--- 6 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 7000){if(output$pressure[output$altitude-output$altitude[1]==7000] > 100){text(-29, altitude_to_pressure(7000), paste0("--- 7 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 8000){if(output$pressure[output$altitude-output$altitude[1]==8000] > 100){text(-29, altitude_to_pressure(8000), paste0("--- 8 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 9000){if(output$pressure[output$altitude-output$altitude[1]==9000] > 100){text(-29, altitude_to_pressure(9000), paste0("--- 9 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 10000){if(output$pressure[output$altitude-output$altitude[1]==10000] > 100){text(-29, altitude_to_pressure(10000), paste0("--- 10 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 11000){if(output$pressure[output$altitude-output$altitude[1]==11000] > 100){text(-29, altitude_to_pressure(11000), paste0("--- 11 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 12000){if(output$pressure[output$altitude-output$altitude[1]==12000] > 100){text(-29, altitude_to_pressure(12000), paste0("--- 12 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 13000){if(output$pressure[output$altitude-output$altitude[1]==13000] > 100){text(-29, altitude_to_pressure(13000), paste0("--- 13 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 14000){if(output$pressure[output$altitude-output$altitude[1]==14000] > 100){text(-29, altitude_to_pressure(14000), paste0("--- 14 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 15000){if(output$pressure[output$altitude-output$altitude[1]==15000] > 100){text(-29, altitude_to_pressure(15000), paste0("--- 15 km"), pos=4, cex = 0.65, col = "black")}}
  if(max(output$altitude-output$altitude[1]) > 16000){if(output$pressure[output$altitude-output$altitude[1]==16000] > 100){text(-29, altitude_to_pressure(16000), paste0("--- 16 km"), pos=4, cex = 0.65, col = "black")}}

  ###
  par(fig = c(0.46, 0.57, 0.03, 0.95), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  over100 = output[["pressure"]] > 100
  sounding_barbs(pressure = output[["pressure"]][over100],
                 ws = output[["ws"]][over100], 
                 wd = output[["wd"]][over100], 
                 altitude = output[["altitude"]][over100], 
                 convert = FALSE, barb_cex = 0.8)
  
  #### draw table ####
  par(fig = c(0.54, 0.99, 0.063, 0.48), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))  
  box()

  ###
  FONTSIZE = 0.64
  
  ###
  text(0.04, 43.1, "MIXR", cex = FONTSIZE, adj = c(0,0))
  text(0.04, 41.1, "[g/kg]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.124, 43.1, "CAPE", cex = FONTSIZE, adj = c(0,0))
  text(0.124, 41.1, "[J/kg]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.21, 43.1, "CAPE03", cex = FONTSIZE, adj = c(0,0))
  text(0.21, 41.1, "[J/kg]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.325, 43.1, "CAPEHGL", cex = FONTSIZE, adj = c(0,0))
  text(0.325, 41.1, "[J/kg]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.46, 43.1, "CIN", cex = FONTSIZE, adj = c(0,0))
  text(0.46, 41.1, "[J/kg]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.54, 43.1, "LI", cex = FONTSIZE, adj = c(0,0))
  text(0.54, 41.1, "[K]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.60, 43.1, "LCL", cex = FONTSIZE, adj = c(0,0))
  text(0.60, 41.1, "[m]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.68, 43.1, "LFC", cex = FONTSIZE, adj = c(0,0))
  text(0.68, 41.1, "[m]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.76, 43.1, "EL", cex = FONTSIZE, adj = c(0,0))
  text(0.76, 41.1, "[m]", cex = FONTSIZE-0.075, adj = c(0,0))
  text(0.855, 43.1, "WMAXSHEAR", cex = FONTSIZE, adj = c(0,0))
  text(0.855, 41.1, "[m2/s2]", cex = FONTSIZE-0.075, adj = c(0,0))
  
  text(-0.016, 37.9, "SB", cex = FONTSIZE, adj = c(0,0))
  text(-0.021, 34.7, "MU", cex = FONTSIZE, adj = c(0,0))
  text(-0.016, 31.5, "ML", cex = FONTSIZE, adj = c(0,0))
  
  text(0.04, 37.9, sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "SB_MIXR")],digits=1)), cex = FONTSIZE,adj = c(0,0))
  text(0.04, 34.7, sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "MU_MIXR")],digits=1)), cex = FONTSIZE,adj = c(0,0))
  text(0.04, 31.5, sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "ML_MIXR")],digits=1)), cex = FONTSIZE,adj = c(0,0))
  
  text(0.125, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.125, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.125, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.21, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_03km_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.21, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_03km_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.21, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_03km_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.325, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_HGL_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.325, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_HGL_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.325, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_HGL_CAPE")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.46, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_CIN")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.46, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_CIN")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.46, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_CIN")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.54, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_LI")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.54, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_LI")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.54, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_LI")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.60, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_LCL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.60, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_LCL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.60, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_LCL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.68, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_LFC_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.68, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_LFC_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.68, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_LFC_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.76, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_EL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.76, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_EL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.76, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_EL_HGT")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.855, 37.9, round(parametry[which(names(parametry[1:LP]) == "SB_WMAXSHEAR")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.855, 34.7, round(parametry[which(names(parametry[1:LP]) == "MU_WMAXSHEAR")],digits=0), cex = FONTSIZE,adj = c(0,0))
  text(0.855, 31.5, round(parametry[which(names(parametry[1:LP]) == "ML_WMAXSHEAR")],digits=0), cex = FONTSIZE,adj = c(0,0))
  
  text(0.925, 37.9, paste0("(E ",round(parametry[which(names(parametry[1:LP]) == "SB_EFF_WMAXSHEAR")],digits=0),")"), cex = FONTSIZE,adj = c(0,0))
  text(0.925, 34.7, paste0("(E ",round(parametry[which(names(parametry[1:LP]) == "MU_EFF_WMAXSHEAR")],digits=0),")"), cex = FONTSIZE,adj = c(0,0))
  text(0.925, 31.5, paste0("(E ",round(parametry[which(names(parametry[1:LP]) == "ML_EFF_WMAXSHEAR")],digits=0),")"), cex = FONTSIZE,adj = c(0,0))
  
  abline(30.25,0)
  
  ###
  
  text(0.171, 27.55, "Bulk wind shear", cex = FONTSIZE, adj = c(1,0))
  text(0.171, 25.55, "[m/s]", cex = FONTSIZE-0.075, adj = c(1,0))
  
  text(0.11,22.55, "Sfc - 1 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.12,22.55, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_01km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,19.15, "Sfc - 3 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.12,19.15, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_03km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,15.75, "Sfc - 6 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.12,15.75, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_06km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,12.35, "Sfc - 8 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.12,12.35, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_08km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,8.95, "Sfc - HGL:", cex = FONTSIZE, adj = c(1,0))
  text(0.12,8.95, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_SFC_to_M10")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,5.55, "Effec. (SB):", cex = FONTSIZE, adj = c(1,0))
  text(0.12,5.55, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_EFF_SB")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,2.15, "Effec. (MU):", cex = FONTSIZE, adj = c(1,0))
  text(0.12,2.15, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_EFF_MU")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.11,-1.25, "Effec. (ML):", cex = FONTSIZE, adj = c(1,0))
  text(0.12,-1.25, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "BS_EFF_ML")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  segments(0.19, -3, 0.19, 30.22)

  ###

  text(0.4375, 27.55, "SRH RM", cex = FONTSIZE, adj = c(1,0))
  text(0.4375, 25.55, "[m2/s2]", cex = FONTSIZE-0.075, adj = c(1,0))
  
  text(0.55, 27.55, "SRH LM", cex = FONTSIZE, adj = c(1,0))
  text(0.55, 25.55, "[m2/s2]", cex = FONTSIZE-0.075, adj = c(1,0))
  
  text(0.345,22.55, "Sfc - 100 m:", cex = FONTSIZE, adj = c(1,0))
  text(0.36,22.55, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_100m_RM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  text(0.4725,22.55, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_100m_LM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  # 
  text(0.345,19.15, "Sfc - 500 m:", cex = FONTSIZE, adj = c(1,0))
  text(0.36,19.15, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_500m_RM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  text(0.4725,19.15, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_500m_LM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.345,15.75, "Sfc - 1 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.36,15.75, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_1km_RM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  text(0.4725,15.75, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_1km_LM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.345,12.35, "Sfc - 3 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.36,12.35, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_3km_RM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  text(0.4725,12.35, paste0(round(parametry[which(names(parametry[1:LP]) == "SRH_3km_LM")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  segments(0.565, 10.25, 0.565, 30.22)
  segments(0.19, 10.25, 0.6, 10.25)
  
  ###
  
  text(0.77, 27.55, "Mean wind", cex = FONTSIZE, adj = c(1,0))
  text(0.77, 25.55, "[m/s]", cex = FONTSIZE-0.075, adj = c(1,0))
  
  text(0.71,22.55, "Sfc - 1 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.72,22.55, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "MW_01km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.71,19.15, "Sfc - 2 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.72,19.15, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "MW_02km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.71,15.75, "1 - 3 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.72,15.75, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "MW_13km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.71,12.35, "Sfc - 6 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.72,12.35, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "MW_06km")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  ###
  
  text(1.025, 27.55, "Lapse rate", cex = FONTSIZE, adj = c(1,0))
  text(1.025, 25.55, "[K/km]", cex = FONTSIZE-0.075, adj = c(1,0))
  
  text(0.95,22.55, "Sfc - 1 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.96,22.55, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "LR_01km")]*-1,digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.95,19.15, "Sfc - 3 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.96,19.15, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "LR_03km")]*-1,digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.95,15.75, "3 - 6 km:", cex = FONTSIZE, adj = c(1,0))
  text(0.96,15.75, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "LR_36km")]*-1,digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.95,12.35, "500700 hPa:", cex = FONTSIZE, adj = c(1,0))
  text(0.96,12.35, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "LR_500700hPa")]*-1,digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  segments(0.7875, 10.25, 0.7875, 30.22)
  segments(0.6, 10.25, 1.2, 10.25)
  
  ####
  
  text(0.43,6.5, "Precip. water [mm]:", cex = FONTSIZE, adj = c(1,0))
  text(0.44,6.5, paste0(round(parametry[which(names(parametry[1:LP]) == "PRCP_WATER")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.43,2.9, "2 - 5 km RH [%]:", cex = FONTSIZE, adj = c(1,0))
  text(0.44,2.9, paste0(round(parametry[which(names(parametry[1:LP]) == "RH_25km")]*100,digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.43,-0.7, "Sfc - 2 km RH [%]:", cex = FONTSIZE, adj = c(1,0))
  text(0.44,-0.7, paste0(round(parametry[which(names(parametry[1:LP]) == "RH_02km")]*100,digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  segments(0.50, 10.25, 0.50, -3)
  
  ###
  
  text(0.78,6.3, "Moisture flux [g/s/m2]:", cex = FONTSIZE, adj = c(1,0))
  text(0.79,6.5, paste0(round(parametry[which(names(parametry[1:LP]) == "Moisture_Flux_02km")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.78,2.9, "4 km DCAPE [J/kg]:", cex = FONTSIZE, adj = c(1,0))
  text(0.79,2.9, paste0(round(parametry[which(names(parametry[1:LP]) == "DCAPE")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  text(0.78,-0.7, "4 km delta theta-e [K]:", cex = FONTSIZE, adj = c(1,0))
  text(0.79,-0.7, paste0(round(parametry[which(names(parametry[1:LP]) == "Delta_thetae")],digits=0),""), cex = FONTSIZE,adj = c(0,0))
  
  segments(0.865, 10.25, 0.865, -3)
  
  ###
  
  text(0.955,6.5, "SHIP:", cex = FONTSIZE, adj = c(1,0))
  text(0.965,6.5, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "SHIP")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.955,2.9, "SCP:", cex = FONTSIZE, adj = c(1,0))
  text(0.965,2.9, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "SCP_new")],digits=1))), cex = FONTSIZE,adj = c(0,0))
  
  text(0.955,-0.7, "STP:", cex = FONTSIZE, adj = c(1,0))
  text(0.965,-0.7, sprintf("%.1f",(round(parametry[which(names(parametry[1:LP]) == "STP_new")],digits=1))), cex = FONTSIZE,adj = c(0,0))

  ###
  par(fig=c(0.54, 0.99, 0, 0.06), new = TRUE, mar = c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
  mtext(expression(paste(bold("thundeR")," - rawinsonde processing tool for R v1.0 (2021)")),
        line = -1.2, cex=0.7)
  #mtext("www.rawinsonde.com",line = -1.8, cex=0.7)
  ####

  ### RH and Theta-E profiles ####
  par(fig = c(0.54, 0.68, 0.49, 0.69875), new = TRUE, 
      mar = c(0.55, 0, 0, 0), oma = c(0, 0, 0, 0))  
  plot(THETAE, output$altitude-output$altitude[1],xaxt="n",yaxt="n",xlab="",  
       type='l',xlim=c(245,385),ylim=c(100,7500),lwd=2,col='magenta')
  axis(4, at=seq(0,7000,1000), las = 1, padj=-0.35, hadj=1.25, xpd = TRUE, 
       labels=c("sfc ","1 km","2 km","3 km","4 km","5 km","6 km","7 km"), 
       cex.axis=0.56, tck=0.1, lwd=0.35)
  axis(1, at=seq(245,385,25), xpd = TRUE, padj=-1.1, 
       cex.axis=0.56, tck=-0.04, lwd=0.35)
  text(315,7100, "Theta-e [K]", cex = 0.65, col = "black")
  
  ####
  par(fig = c(0.54, 0.68, 0.70875, 0.9175), new = TRUE, 
      mar = c(0.55, 0, 0, 0), oma = c(0, 0, 0, 0))  
  plot(RH, output$altitude-output$altitude[1],xaxt="n",yaxt="n",xlab="",  
       type='l',xlim=c(0,119),ylim=c(100,7500),lwd=2,col='blue')
  axis(4, at=seq(0,7000,1000), las = 1, padj=-0.35, hadj=1.25, xpd = TRUE, 
       labels=c("sfc ","1 km","2 km","3 km","4 km","5 km","6 km","7 km"), 
       cex.axis=0.56, tck=0.1, lwd=0.35)
  axis(1, at=seq(0,119,20), xpd = TRUE, padj=-1.1, 
       cex.axis=0.56, tck=-0.04, lwd=0.35)
  text(60,7100, "RH [%]", cex = 0.65, col = "black")

  ####

  par(fig = c(0.69, 0.99, 0.49, 0.9175), new = TRUE, 
      mar = c(0, 0, 0, 0), oma=c(0, 0, 0, 0))

  sounding_hodograph(output$ws, output$wd, output$altitude-output$altitude[1], max_speed = max_speed, frame = FALSE, close_par = FALSE, ...)

  up = round(max_speed * -1.55  * cos(135 * pi/180), 2)
  vp = round(max_speed * 1.35 * sin(135 * pi/180), 2)
  points(up, vp, pch = 19, col = "white", cex = 2)
  text(up, vp, labels = "[m/s]", col = "gray20", cex = 0.75)
  
  RM_y = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_RM_M")] * cos(parametry[which(names(parametry[1:LP]) == "Bunkers_RM_A")] * pi/180), 2)
  RM_x = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_RM_M")] * sin(parametry[which(names(parametry[1:LP]) == "Bunkers_RM_A")] * pi/180), 2)  # LM_x = m_x1 - 7.5 * (sqrt( 1 / (1+slope2*slope2)))
  LM_y = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_LM_M")] * cos(parametry[which(names(parametry[1:LP]) == "Bunkers_LM_A")] * pi/180), 2)
  LM_x = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_LM_M")] * sin(parametry[which(names(parametry[1:LP]) == "Bunkers_LM_A")] * pi/180), 2)  # LM_x = m_x1 - 7.5 * (sqrt( 1 / (1+slope2*slope2)))
  MW_y = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_MW_M")] * cos(parametry[which(names(parametry[1:LP]) == "Bunkers_MW_A")] * pi/180), 2)
  MW_x = round(-parametry[which(names(parametry[1:LP]) == "Bunkers_MW_M")] * sin(parametry[which(names(parametry[1:LP]) == "Bunkers_MW_A")] * pi/180), 2)  # LM_x = m_x1 - 7.5 * (sqrt( 1 / (1+slope2*slope2)))

  vectors_color = rgb(128,128,128, maxColorValue = 255, alpha = 150)
  circle_color = rgb(255,255,255, maxColorValue = 255, alpha = 150)
  points(LM_x,LM_y, cex = 1.75, pch = 1, col=vectors_color)
  points(LM_x,LM_y, cex = 0.75, pch = 1, col='black')
  points(RM_x,RM_y, cex = 1.75, pch = 1, col=vectors_color)
  points(RM_x,RM_y, cex = 0.75, pch = 1, col='black')
  arrows(0,0,MW_x,MW_y, lty=1, length = 0.11, angle = 22, code = 2, col = rgb(75,75,75, maxColorValue = 255, alpha = 100), lwd = 2.75)
  text(LM_x,LM_y, font=2,"LM", col = vectors_color, adj = c(0.5, -1.5),cex = 0.5)
  text(RM_x,RM_y, font=2,"RM", col = vectors_color, adj = c(0.5, 2.5),cex = 0.5)

  ####

  max_speed = max_speed*1.25

  rect(round(max_speed * 1.38  * cos(135 * pi/180), 2), 
       round(max_speed * -0.72  * sin(135 * pi/180), 2), 
       round(max_speed * 0.02  * cos(135 * pi/180), 2),
       round(max_speed * -1.17  * sin(135 * pi/180), 2),
       lty = par("lty"), lwd = 0, col=rgb(255,255,255, maxColorValue = 255, alpha = 200))
  
  text(round(max_speed * 1.37  * cos(135 * pi/180), 2), 
       round(max_speed * -1.14 * sin(135 * pi/180), 2),  
       adj=c(0,0),
       labels = paste0("Right-moving: ",sprintf("%.0f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_RM_A")],digits = 0))," / ", sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_RM_M")],digits = 1))," m/s"), col = "black", cex = 0.66)
  text(round(max_speed * 1.37 * cos(135 * pi/180), 2), 
       round(max_speed * -0.99 * sin(135 * pi/180), 2),  
       adj=c(0,0),
       labels = paste0("Storm-motion: ",sprintf("%.0f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_MW_A")],digits = 0))," / ", sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_MW_M")],digits = 1))," m/s"), col = "black", cex = 0.66)
  text(round(max_speed * 1.37  * cos(135 * pi/180), 2), 
       round(max_speed * -0.84 * sin(135 * pi/180), 2),  
       adj=c(0,0),
       labels = paste0("Left-moving: ",sprintf("%.0f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_LM_A")],digits = 0))," / ", sprintf("%.1f",round(parametry[which(names(parametry[1:LP]) == "Bunkers_LM_M")],digits = 1))," m/s"), col = "black", cex = 0.66)

#box()

####

  title(title,outer=T,line=-1.5)

}
