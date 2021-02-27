# thundeR

###### Rapid computation and visualisation of convective parameters from rawinsonde and NWP data <img src="man/figures/logo.png" align="right" width="180" /> 

<!-- badges: start --> 
[![R-CMD-check](https://github.com/bczernecki/thunder/workflows/R-CMD-check/badge.svg)](https://github.com/bczernecki/thunder/actions)
[![Codecov test coverage](https://codecov.io/gh/bczernecki/thunder/branch/devel/graph/badge.svg?token=JGZPB7RUFI)](https://codecov.io/gh/bczernecki/thunder)
[![CRAN status](https://www.r-pkg.org/badges/version/climate)](https://cran.r-project.org/package=climate)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/climate)](https://cran.r-project.org/package=climate)
<!-- badges: end -->

**`thundeR`** is a freeware R package and collection of functions for rapid computation and visualisation of convective parameters commonly used in the operational forecasting of severe convective storms. Core algorithm is based on C++ code implemented into R language via RCPP. This solution allows to compute over 100 thermodynamic and kinematic parameters in less than 0.02s per profile and process large datasets such as reanalyses or operational NWP models in a reasonable amount of time. Package has been developed since 2017 by research meteorologists specializing in severe convective storms and is constantly updated with new features.


## Online browser
------------

Online rawinsonde browser of **thundeR** package is available at www.rawinsonde.com 


## Installation
------------

The stable version can be installed from the CRAN repository:

``` r
install.packages("thunder")
```


The development version can be installed directly from the github repository:

``` r
library(remotes);install_github("bczernecki/thunder")
```

## Usage
-----

### Draw Skew-T diagram, hodograph and selected convective parameters on a single layout and save it to png file

``` r
data("sounding_wien") # load example dataset (rawinsonde observation):
pressure = sounding_wien$PRES # vector of pressure [hPa]
altitude = sounding_wien$HGHT # vector of altitude [m above ground level]
temp = sounding_wien$TEMP  # vector of air temperature [°C]
dpt = sounding_wien$DWPT # vector of dew point temperature [°C]
wd = sounding_wien$DRCT # vector of wind direction as azimuth [°]
ws = sounding_wien$SKNT # vector of wind speed [kn]
sounding_save(filename = "myfile.png", title = "Vienna - 23 August 2011 1200 UTC", pressure, altitude, temp, dpt, wd, ws)
```

![](inst/figures/Vienna_profile.png)



### Compute all convective indices in a form of numeric vector based on a sample vertical profile data:

``` r
library("thunder")
pressure <- c(1000, 855, 700, 500, 300, 100, 10) # pressure [hPa]
altitude <- c(0, 1500, 2500, 6000, 8500, 12000, 25000) # altitude [m above ground level]
temp <- c(25, 10, 0, -15, -30, -50, -92) # air temperature [°C]
dpt <- c(20, 5, -5, -30, -55, -80, -99) # dew point temperature [°C]
wd <- c(0, 90, 135, 180, 270, 350, 0) # wind direction [°]
ws <- c(5, 10, 20, 30, 40, 5, 0) # wind speed [kn]
accuracy <- 2 # accuracy of computations where 3 = high (slow), 2 = medium (recommended), 1 = low (fast)
options(digits = 2) # change output formatting precision 
sounding_compute(pressure, altitude, temp, dpt, wd, ws, accuracy)


#             MU_CAPE        MU_03km_CAPE         MU_HGL_CAPE              MU_CIN 
#             2195.41              575.59             1568.01                0.00 
#          MU_LCL_HGT          MU_LFC_HGT           MU_EL_HGT               MU_LI 
#              730.00              730.00             8200.00               -9.63 
#             MU_WMAX          MU_EL_TEMP         MU_LCL_TEMP         MU_LFC_TEMP 
#               66.26              -28.20               17.70               17.70 
#             MU_MIXR             SB_CAPE        SB_03km_CAPE         SB_HGL_CAPE 
#               14.88             2195.41              575.59             1568.01 
#              SB_CIN          SB_LCL_HGT          SB_LFC_HGT           SB_EL_HGT 
#                0.00              730.00              730.00             8200.00 
#               SB_LI             SB_WMAX          SB_EL_TEMP         SB_LCL_TEMP 
#               -9.63               66.26              -28.20               17.70 
#         SB_LFC_TEMP             SB_MIXR             ML_CAPE        ML_03km_CAPE 
#               17.70               14.88             1559.41              416.88 
#         ML_HGL_CAPE              ML_CIN          ML_LCL_HGT          ML_LFC_HGT 
#             1187.96                0.00              975.00              975.00 
#           ML_EL_HGT               ML_LI             ML_WMAX          ML_EL_TEMP 
#             7800.00               -7.15               55.85              -25.80 
#         ML_LCL_TEMP         ML_LFC_TEMP             ML_MIXR             LR_01km 
#               15.25               15.25               13.02              -10.00 
#             LR_03km             LR_24km             LR_36km        LR_500700hPa 
#               -9.05               -5.77               -4.29               -4.29 
#        LR_500800hPa            FRZG_HGT    FRZG_wetbulb_HGT HGT_max_thetae_03km 
#               -6.67             2500.00             2300.00                0.00 
# HGT_min_thetae_04km        Delta_thetae               DCAPE  Cold_Pool_Strength 
#             3700.00               28.46              595.13               12.77 
#        Wind_Index            PRCP_WATER  Moisture_Flux_02km             RH_02km 
#               34.12               27.10               28.49                0.72 
#             RH_25km              RH_HGL             BS_01km             BS_02km 
#                0.58                0.46                3.83                8.78 
#             BS_03km             BS_06km             BS_08km             BS_36km 
#               12.66               18.01               17.41                9.37 
#             BS_18km           BS_EFF_MU           BS_EFF_SB           BS_EFF_ML 
#               20.28               14.14               14.14               13.82 
#       BS_SFC_to_HGL    BS_MU_LFC_to_HGL    BS_SB_LFC_to_HGL    BS_ML_LFC_to_HGL 
#               15.51               14.07               14.07               13.69 
#             MW_01km             MW_02km             MW_06km             MW_13km 
#                2.36                2.81                5.14                6.88 
#         SRH_100m_RM         SRH_500m_RM          SRH_1km_RM          SRH_3km_RM 
#                2.87               14.37               29.47              136.42 
#         SRH_100m_LM         SRH_500m_LM          SRH_1km_LM          SRH_3km_LM 
#                0.30                1.51                3.10              -30.53 
#             K_Index     Showalter_Index   TotalTotals_Index         SWEAT_Index 
#               24.35                3.90               44.35              106.42 
#                 STP             STP_new                 SCP             SCP_new 
#                0.26                0.14                5.39                4.23 
#                SHIP                 DCP        MU_WMAXSHEAR        SB_WMAXSHEAR 
#                0.61                0.73             1193.11             1193.11 
#        ML_WMAXSHEAR    MU_EFF_WMAXSHEAR    SB_EFF_WMAXSHEAR    ML_EFF_WMAXSHEAR 
#             1005.54              936.94              936.94              771.71
```

Developers
-------------

**thundeR** package has been developed by atmospheric scientists, each having an equal contribution (listed in alphabetical order):
- Bartosz Czernecki (Adam Mickiewicz University in Poznań, Poland)
- Piotr Szuster (University of Technology in Cracow, Poland)
- Mateusz Taszarek (CIMMS/NSSL in Norman, Oklahoma, United States)


Contributions
-------------

[Feel free to submit issues and enhancement requests.](https://github.com/bczernecki/thunder/issues)
