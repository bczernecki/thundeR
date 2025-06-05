# thundeR 1.1.5

* Fixes for SRH-related calculations that were not working properly in thundeR 1.1.4
* Added tests to check whether calculations of indices matches expected results
* Removed dependency to the RadioSonde package due to "heavy"" OS-related dependencies on some Linux distributions
* Minor bug fixes


# thundeR 1.1.4

* Fixing minor bugs in C++ code that caused the package to show warnings uninitialized variable with GCC 14.2 on Windows 
* Adding RadioSonde package as dependency
* Minor bug fixes and changes in CI/CD pipeline

# thundeR 1.1.3

* Adding new tests for thermodynamic indices failing with metPy 1.5, but working smoothly with thundeR and sharppy

# thundeR 1.1.2

* Allowing to keep NA values in `get_sounding()` which might be not harmful if missing data occur on higher altitudes for dew point temperatures

# thundeR 1.1.1

* Improvements in C++ algorithm responsible for more accurate computations of parcel trajectories and related thermodynamic indices
  * Minor bug fixes in CIN calculation procedure
  * Minor bug fixes in mixed-layer calculation procedure
  * Minor bug fixes in DCAPE calculation procedure
* Fixing code leaks in C++ module
* New parameters added. Currently the list of computed parameters contains over 200 indices
* Documentation and pkgdown site updates
* Removing dependency to the `climate` package by including internal functions responsible for downloading soundings from Wyoming database
* Improving graphical layout:
  * Adding parcel trajectories shading for positive area
  * shading SRH area for left- and right movers on sounding hodographs
* Implementation of different interpolating mechanisms to `sounding_compute()`
* UX improvements and interactive editing of rawinsonde data browsed in shiny app (http://rawinsonde.com)

# thundeR 0.3.0

* Fixing code leaks in C++ module (big thanks to Jakub Nowosad!)
* Adding new convective parameters
* Adding set of parameters for interpolation accuracy that allow to control computation performance
* Further improvements in graphical layout for `sounding_plot()`, `skewt_plot()` and `sounding_hodograph`

# thundeR 0.0.12

* First public release on CRAN
