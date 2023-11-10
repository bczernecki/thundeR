
context("sounding_default and sounding_compute should be equal")

test_that("sounding_default should return the same values as sounding_compute for thermodynamics", {
        pressure <- c(1000, 855, 700, 500, 300, 100, 10)
        altitude <- c(0, 1500, 2500, 6000, 8500, 12000, 25000)
        temp <- c(25, 10, 0, -15, -30, -50, -92)
        dpt <- c(20, 5, -5, -30, -55, -80, -99)
        wd <- c(0, 90, 135, 180, 270, 350, 0)
        ws <- c(5, 10, 20, 30, 40, 5, 0)
        options(scipen = 999) # change formatting
        res_default = sounding_default(pressure, altitude, temp, dpt, wd, ws,
                                       export_profile = 0,
                                       accuracy = 2,
                                       interpolate_step = 1000,
                                       meanlayer_bottom_top = c(0, 500),
                                       storm_motion = c(999, 999))
        res_compute = sounding_compute(pressure, altitude, temp, dpt, wd, ws,
                                       accuracy = 2, 
                                       interpolate_step = 1000,
                                       meanlayer_bottom_top = c(0, 500),
                                       storm_motion = c(999, 999))
        
        indices_cape = grep(x = names(res_compute), "CAPE")
        indices_cin = grep(x = names(res_compute), "CIN")
        
        expect_equal(as.numeric(res_compute[indices_cape]), res_default[indices_cape])
        expect_equal(as.numeric(res_compute[indices_cin]), res_default[indices_cin])
})



test_that("sounding_default should have CAPE 2k+ values for profile having CAPE=1 in MetPy", {
  p = c(1009.57453125, 1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 
        350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50)
  t = c(25.92, 25.07, 23.13, 21.55, 20.54, 19.39, 18.19, 16.96, 15.84, 14.75, 13.55, 12.3, 9.39, 6.73, 3.08, 0.01,
        -4.14, -8.98, -14.8, -22.84, -32.18, -42.92, -48.6, -53.6, -57.96, -67.31, -72.61, -76.33, -75.49, -64.64) 
  alt = c(0, 84, 307, 534, 766, 1003, 1246, 1495, 1750, 2012, 2281, 2557, 3135, 3748, 4402, 5103, 5860, 6683, 7584, 
          8579, 9688, 10945, 11646, 12411, 13261, 14212, 15297, 16586, 18636, 20653)
  q = c(24.13, 21.86, 21.28, 21.13, 20.07, 18.6, 17.6, 16.7, 15.79, 14.75, 13.55, 12.21, 8.52, 3.77, -4.53, -16.51, 
        -21.2, -23.21, -32.89, -40.98, -50.96, -59.87, -64.33, -66.62, -65.41, -75.21, -82.58, -85.87, -87.36, -86.24)
  ws = runif(length(q))
  wd = runif(length(q))

  res_default = sounding_default(p, alt, t, q, wd, ws,
                                 export_profile = 0,
                                 accuracy = 2,
                                 interpolate_step = 10,
                                 meanlayer_bottom_top = c(0, 500),
                                 storm_motion = c(999, 999))
  res_compute = sounding_compute(p, alt, t, q, wd, ws,
                                 accuracy = 2, 
                                 interpolate_step = 1000,
                                 meanlayer_bottom_top = c(0, 500),
                                 storm_motion = c(999, 999))
  
  indices_cape = grep(x = names(res_compute), "SB_CAPE$")
  expect_equal(as.numeric(res_compute[indices_cape]), res_default[indices_cape])
  expect_gt(res_default[indices_cape], 2000)
  
})