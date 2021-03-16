
context("Check whether there's message if wind speed is currently not supported")

test_that ("windbarbs message", {
  plot(1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", frame = FALSE)
  expect_message(windbarbs(cx = 1, cy = 1, direction = 120, speed = 199, cex = 5), 
                 "^Currently wind barbs are supported only up to 190 knots. Your wind speed is 199\\n")
})
