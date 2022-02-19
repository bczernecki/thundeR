
test_that("sounding_save test png", {

  data("sounding_vienna")
  attach(sounding_vienna)
  fname = paste0(tempfile(), ".png")
  sounding_save(filename = fname, 
                pressure, altitude, temp, dpt, wd, ws, parcel = "MU", 
                title = "Vienna - 23 August 2011, 12:00 UTC")
  
  expect_true(file.size(fname) > 0)
  
})


# TODO:
# temporarily this test is turned off due to CI/CD for MacOS that
# does not support SVG graphics
# test_that("sounding_save test svg", {
#   
#   data("sounding_vienna")
#   attach(sounding_vienna)
#   fname = paste0(tempfile(), ".svg")
#   sounding_save(filename = fname, 
#                 pressure, altitude, temp, dpt, wd, ws, parcel = "MU", 
#                 title = "Vienna - 23 August 2011, 12:00 UTC")
#   
#   expect_true(file.size(fname) > 0)
#   
# })
# 
