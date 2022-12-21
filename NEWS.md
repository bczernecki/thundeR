# thundeR 1.1.0

* Improvements in C++ algorithm responsible for more accurate computations of parcel trajectories and related thermodynamic indices
* Removing dependency to the `climate` package
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

* Public release on CRAN.
