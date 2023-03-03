test_that("Downloading sounding profile", {
  
  profile = get_sounding(wmo_id = 12120, 
                         yy = 2010,
                         mm = 8, 
                         dd = 20, 
                         hh = 12)
  
  expect_equal(nrow(profile), 72)
  expect_true(is.data.frame(profile))
})

test_that("Sounding profile error", {
  
  expect_error(get_sounding(wmo_id = 99999, 
                         yy = 2010,
                         mm = 8, 
                         dd = 20, 
                         hh = 12)
               )
  
})