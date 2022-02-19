test_that ("Check example for sounding_plot", {
  
  data("sounding_vienna")
  sounding_vienna <- na.omit(sounding_vienna)
  
  
  expect_message(sounding_plot(sounding_vienna$pressure, sounding_vienna$altitude,
                              sounding_vienna$temp, sounding_vienna$dpt,
                              sounding_vienna$wd, sounding_vienna$ws,
                              parcel = "MU", title = "Vienna - 23 August 2011, 12:00 UTC"
  ))
  
  expect_message(sounding_plot(sounding_vienna$pressure, sounding_vienna$altitude,
                               sounding_vienna$temp, sounding_vienna$dpt,
                               sounding_vienna$wd, sounding_vienna$ws,
                               parcel = "MU", 
                               title = "Vienna - 23 August 2011, 12:00 UTC", hazards = TRUE
                              
  ))
  
})
