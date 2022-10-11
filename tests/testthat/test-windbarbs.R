
context("Check whether there's message if wind speed is currently not supported")

test_that ("windbarbs message", {
  plot(1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", frame = FALSE)
  expect_message(windbarbs(cx = 1, cy = 1, direction = 120, speed = 199, cex = 5), 
                 "^Currently wind barbs are supported only up to 190 knots. Your wind speed is 199\\n")
})

test_that("windbarbs as in examples", {
  
  expect_silent(for (i in 1:38){
      sc = 5
       plot(0:2, xaxt = 'n', yaxt = 'n', type = "n", xlab = "", ylab = "")
       text(1.4,1, i*sc, cex = 1.5)
       windbarbs(cx = 2, cy = 1, direction = 60, speed = i*sc, cex = 3)
  })
})
