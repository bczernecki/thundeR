#' Calculate X and Y-coordinates on Skew-T diagram
#'
#' @noRd
#' 

skewtx = function(temp, ycoord) {
        0.54 * temp + 0.90692 * ycoord
}

skewty = function(pres) {
        132.182 - 44.061 * log10(pres)
}


#' Calculate satlft
#' 
#' Taken from the R Radiosonde package. Originally from the "thermo.f" routines and technical documentation
#' National Oceanic and Atmospheric Administration 
#' HM Stone 1986 — Part A:--Program Information and Installation Procedure. 
#' Program Name; CONVECT.SV. AAL ID;. Revision No.; 01.00.
#'
#' @noRd
#' 

satlft = function(thw, p) {
        cta <- 273.15
        akap <- 0.28541
        pwrp <- (p/1000)^akap
        tone <- (thw + cta) * pwrp - cta
        eone <- wobf(tone) - wobf(thw)
        rate <- 1
        dlt <- 1
        while (abs(dlt) > 0.1) {
                ttwo <- tone - eone * rate
                pt <- (ttwo + cta)/pwrp - cta
                etwo <- pt + wobf(ttwo) - wobf(pt) - thw
                dlt <- etwo * rate
                rate <- (ttwo - tone)/(etwo - eone)
                tone <- ttwo
                eone <- etwo
        }
        ttwo - dlt
}

#' calculate TDA - for Skew-T diagram
#' 
#'
#' @noRd
#' 

tda = function(o, p) {
        ok <- o + 273.15
        tdak <- ok * ((p * 0.001)^0.286)
        tdak - 273.15
}

#' Calculate wobf
#' 
#' Taken from the R Radiosonde package. Originally from the "thermo.f" routines and technical documentation
#' National Oceanic and Atmospheric Administration 
#' HM Stone 1986 — Part A:--Program Information and Installation Procedure. 
#' Program Name; CONVECT.SV. AAL ID;. Revision No.; 01.00.
#'
#' @noRd
#' 
wobf <- function (temp) {
        x <- temp - 20
        if (x <= 0) {
         pol <- 1 + x * (-0.0088416605 + x * (0.00014714143 + 
            x * (-9.671989e-07 + x * (-3.2607217e-08 + x * (-3.8598073e-10)))))
        wbts <- 15.13/pol^4
        }
        else {
         pol <- 1 + x * (0.0036182989 + x * (-1.3603273e-05 + 
            x * (4.9618922e-07 + x * (-6.1059365e-09 + x * (3.9401551e-11 + 
                x * (-1.2588129e-13 + x * (1.668828e-16)))))))
         wbts <- 29.93/pol^4 + 0.96 * x - 14.8
        }
        return(wbts)
        }