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
#' Taken from the R Radiosonde package
#'
#' @noRd
#' @importFrom RadioSonde wobf
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
