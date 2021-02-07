"wobf" <-
function(temp)
{
#
# Copyright 2001,2002 Tim Hoar
#
# This file is part of the RadioSonde library for R and related languages.
#
# RadioSonde is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# RadioSonde is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RadioSonde; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

	#-----------------------------------------------------------------------
	#
	# this function calculates the difference of the wet-bulb potential
	# temperatures for saturated and dry air given the temperature.
	#
	#-----------------------------------------------------------------------
	# include 'lib_dev:[gudoc]edfvaxbox.for/list'
	# baker, schlatter  17-may-1982	  original version.
	#      let wbpts = wet-bulb potential temperature for saturated
	# air at temperature t (celsius). let wbptd = wet-bulb potential
	# temperature for completely dry air at the same temperature t.
	# the wobus function wobf (in degrees celsius) is defined by
	#                    wobf(t) = wbpts-wbptd.
	# although wbpts and wbptd are functions of both pressure and
	# temperature, their difference is a function of temperature only.
	#      to understand why, consider a parcel of dry air at tempera-
	# ture t and pressure p. the thermodynamic state of the parcel is
	# represented by a point on a pseudoadiabatic chart. the wet-bulb
	# potential temperature curve (moist adiabat) passing through this
	# point is wbpts. now t is the equivalent temperature for another
	# parcel saturated at some lower temperature tw, but at the same
	# pressure p.  to find tw, ascend along the dry adiabat through
	# (t,p). at a great height, the dry adiabat and some moist
	# adiabat will nearly coincide. descend along this moist adiabat
	# back to p. the parcel temperature is now tw. the wet-bulb
	# potential temperature curve (moist adiabat) through (tw,p) is wbptd.
	# the difference (wbpts-wbptd) is proportional to the heat imparted
	# to a parcel saturated at temperature tw if all its water vapor
	# were condensed. since the amount of water vapor a parcel can
	# hold depends upon temperature alone, (wbptd-wbpts) must depend
	# on temperature alone.
	#      the wobus function is useful for evaluating several thermo-
	# dynamic quantities.  by definition:
	#		    wobf(t) = wbpts-wbptd.               (1)
	# if t is at 1000 mb, then t is a potential temperature pt and
	# wbpts = pt. thus
	#		    wobf(pt) = pt-wbptd.                 (2)
	# if t is at the condensation level, then t is the condensation
	# temperature tc and wbpts is the wet-bulb potential temperature
	# wbpt. thus
	#		    wobf(tc) = wbpt-wbptd.               (3)
	# if wbptd is eliminated from (2) and (3), there results
	#		    wbpt = pt-wobf(pt)+wobf(tc).
	# if wbptd is eliminated from (1) and (2), there results
	#		    wbpts = pt-wobf(pt)+wobf(t).
	#      if t is an equivalent potential temperature ept (implying
	# that the air at 1000 mb is completely dry), then wbpts = ept
	# and wbptd = wbpt. thus
	#		    wobf(ept) = ept-wbpt.
	# this form is the basis for a polynomial approximation to wobf.
	# in table 78 on pp.319-322 of the smithsonian meteorological
	# tables by roland list (6th revised edition), one finds wet-bulb
	# potential temperatures and the corresponding equivalent potential
	# temperatures listed together. herman wobus, a mathematician for-
	# merly at the navy weather research facility, norfolk, virginia, 
	# and now retired, computed the coefficients for the polynomial
	# approximation from numbers in this table.
	#
	#                                  notes by t.w. schlatter
	#                                  noaa/erl/profs program office
	#                                  august 1981
	x <- temp - 20.
	if(x <= 0.) {
		pol <- 1. + x * (-0.0088416604999999992 + x * (
			0.00014714143000000001 + x * (-9.6719890000000006e-07 +
			x * (-3.2607217000000002e-08 + x * (
			-3.8598072999999999e-10)))))
		wbts <- 15.130000000000001/pol^4
	}
	else {
		pol <- 1. + x * (0.0036182989000000001 + x * (-1.3603273e-05 +
			x * (4.9618921999999997e-07 + x * (
			-6.1059364999999998e-09 + x * (3.9401550999999998e-11 +
			x * (-1.2588129e-13 + x * (1.668828e-16)))))))
		wbts <- 29.93/pol^4 + 0.95999999999999996 * x - 
			14.800000000000001
	}
	return(wbts)
}
