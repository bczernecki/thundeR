Sounding - R Package
====================

<!-- badges: start -->
[![R-CMD-check](https://github.com/bczernecki/sounding/workflows/R-CMD-check/badge.svg)](https://github.com/bczernecki/sounding/actions)
<!-- badges: end -->
To be changed when CRAN comes...
[![CRAN status](https://www.r-pkg.org/badges/version/climate)](https://cran.r-project.org/package=climate)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/climate)](https://cran.r-project.org/package=climate)

`sounding` is a freeware R package for performing analyses on atmospheric sounding profiles.

The main core of this computational code is a highly optimized version of C++ code dedicated for calculating sounding derived indices related with atmospheric convections.

Installation
------------

The development version from github:

``` r
library(devtools);install_github("bczernecki/sounding", auth_token = "8caadca559389c2ead3a59edc7c9c1e1ff38eb94")
```

Usage
-----

``` r
library('sounding')
pressure <- c(1000, 855, 700, 500, 300, 100, 10) # pressure
altitude <- c(0, 1500, 2500, 6000, 8500, 12000, 25000) # altitude
temp <- c(25, 10, 0, -15, -30, -50, -92) # air temperature
dpt <- c(20, 5, -5, -30, -55, -80, -99) # dew point temperature
wd <- c(0, 90, 135, 180, 270, 350, 0) # wind direction
ws <- c(5, 10, 20, 30, 40, 5, 0) # wind speed

sounding(pressure, altitude, temp, dpt, wd, ws)
#>                 K_Index         Showalter_Index        Showalter_IndexV 
#>              25.0000000               3.6423438               8.0337671 
#>        MostUnstableCAPE      LLMostUnstableCAPE         MostUnstableCIN 
#>            1987.4458685             487.8597626               0.0000000 
#>         MostUnstableLCL         MostUnstableLFC          MostUnstableEL 
#>             800.0000000             800.0000000            8200.0000000 
#>          MostUnstableLI        MostUnstableVmax       VMostUnstableCAPE 
#>              -8.8967187              63.0467425            2685.1881933 
#>     VLLMostUnstableCAPE        VMostUnstableCIN        VMostUnstableLCL 
#>             725.4059002               0.0000000             800.0000000 
#>        VMostUnstableLFC         VMostUnstableEL         VMostUnstableLI 
#>             800.0000000            8200.0000000              -7.4590191 
#>       VMostUnstableVmax        SurfaceBasedCAPE      LLSurfaceBasedCAPE 
#>              73.2828519            1987.4458685             487.8597626 
#>         SurfaceBasedCIN         SurfaceBasedLCL         SurfaceBasedLFC 
#>               0.0000000             800.0000000             800.0000000 
#>          SurfaceBasedEL          SurfaceBasedLI        SurfaceBasedVmax 
#>            8200.0000000              -8.8967187              63.0467425 
#>       VSurfaceBasedCAPE     VLLSurfaceBasedCAPE        VSurfaceBasedCIN 
#>            2685.1881933             725.4059002               0.0000000 
#>        VSurfaceBasedLCL        VSurfaceBasedLFC         VSurfaceBasedEL 
#>             800.0000000             800.0000000            8200.0000000 
#>         VSurfaceBasedLI       VSurfaceBasedVmax           MeanLayerCAPE 
#>              -7.4590191              73.2828519            1358.5941690 
#>         LLMeanLayerCAPE            MeanLayerCIN            MeanLayerLCL 
#>             344.6297324               0.0000000            1000.0000000 
#>            MeanLayerLFC             MeanLayerEL             MeanLayerLI 
#>            1000.0000000            7700.0000000              -6.2014062 
#>           MeanLayerVmax          VMeanLayerCAPE        VLLMeanLayerCAPE 
#>              52.1266567            1857.4079325             512.0822525 
#>           VMeanLayerCIN           VMeanLayerLCL           VMeanLayerLFC 
#>               0.0000000            1000.0000000            1000.0000000 
#>            VMeanLayerEL            VMeanLayerLI          VMeanLayerVmax 
#>            7700.0000000              -4.1287940              60.9492893 
#>                   DCAPE                  VDCAPE             CrossTotals 
#>             765.4951174            1119.0097413              20.0000000 
#>          VerticalTotals             TotalTotals              SWEATIndex 
#>              25.0000000              45.0000000             110.0000000 
#>             LapseRate01             LapseRate24                 TQIndex 
#>             -10.0000000              -5.9183673              15.0000000 
#>              ZeroHeight       WetBulbZeroHeight                MLHeight 
#>            2500.0000000            2300.0000000             200.0000000 
#>           MinTHTEHeight                    BS01                    BS02 
#>            3700.0000000               3.8344392               8.7821084 
#>                    BS03                    BS36                    BS06 
...
```

Skew-T diagrams:

``` r
x <- data.frame(lev=c(1000, 855, 700, 500, 400, 300, 100, 10),
                HGT=c(0, 1500, 2500, 5300, 6700, 8500, 12000, 25000),
                TMP=c(25, 10, 0, -15, -30, -40, -50, -92),
                DPT = c(20, 5, -5, -30, -55, -65,-80, -99),
                WD = c(0, 90, 135, 180, 270, 290, 350, 0),
                WS = c(5, 10, 20, 30, 40, 45, 5, 0),
                lon=52,
                lat=19,
                date2=Sys.time())

#load("data/holdxy.Rdata")
load(system.file("data/holdxx.rda",  package = "sounding"))
load(system.file("data/holdyy.rda",  package = "sounding"))


plotuj(x, wczytane_holdxy=TRUE, etykiety_hgt = T, zaznacz_cape = T, ramka = T)
```

![](figure/test%20two-1.png)

Contributions
-------------

[Feel free to submit issues and enhancement requests.](https://github.com/bczernecki/sounding/issues)
