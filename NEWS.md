# thundeR 1.1.2

* Allowing to keep NA values in `get_sounding()` which might be not harmful if missing data occur on higher altitudes

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
