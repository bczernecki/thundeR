test_that ("sounding_wind test", {
  
  data("sounding_vienna")
  attach(sounding_vienna)
  fname = paste0(tempfile(), ".png")
  png(fname)
  sounding_wind(pressure = pressure, ws = ws, yaxs = TRUE)
  dev.off()
  
  expect_true(file.size(fname) > 0)
  
})
