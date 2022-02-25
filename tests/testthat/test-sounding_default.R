
context("Computation efficiency")

test_that("sounding_default should compute below 0.5s", {
        pressure <- c(1000, 855, 700, 500, 300, 100, 10)
        altitude <- c(0, 1500, 2500, 6000, 8500, 12000, 25000)
        temp <- c(25, 10, 0, -15, -30, -50, -92)
        dpt <- c(20, 5, -5, -30, -55, -80, -99)
        wd <- c(0, 90, 135, 180, 270, 350, 0)
        ws <- c(5, 10, 20, 30, 40, 5, 0)
        options(scipen = 999) # change formatting
        t1 <- system.time(sounding_default(pressure, altitude, temp, dpt, wd, ws, export_profile = 0, accuracy = 2))
        t1 <- as.numeric(t1[3])
        expect_lt(t1, 2.0)  # changed for valgrind compilator, which might be superslow in debugging mode
        
})
