# derivnum(): the Revelle factor from the derivative of CO2 (or fCO2) with respect to DIC at constant TAlk agrees
# with the analytic buffsun(), and the derivative of pCO2 differs by the fugacity-coefficient term.
test_that("derivnum Rf from CO2 and fCO2 agrees with buffsun, and pCO2 differs by the fugacity term", {
  a <- list(flag = 15, var1 = 2300e-6, var2 = 2000e-6, S = 35, T = 18,
            k1k2 = "l", kf = "pf", ks = "d", b = "u74", pHscale = "T")
  ref <- do.call(buffsun, c(a, warn = "n"))$Rf
  cs  <- do.call(carb, c(a, warn = "n"))
  dn  <- do.call(derivnum, c(list(varid = "var2"), a, list(warn = "n")))
  rf_co2  <- dn$CO2  * a$var2 / cs$CO2
  rf_fco2 <- dn$fCO2 * a$var2 / cs$fCO2
  rf_pco2 <- dn$pCO2 * a$var2 / cs$pCO2
  expect_lt(abs(rf_co2  / ref - 1), 1e-8)    # 0.01 ppm; 0.0013 ppm at the time of writing
  expect_lt(abs(rf_fco2 / ref - 1), 1e-8)
  expect_equal(rf_co2, rf_fco2, tolerance = 1e-9)
  f  <- cs$fCO2 * 1e-6                       # fCO2 in atm
  tk <- 273.15 + a$T
  term <- 4 * (1 - f) * (57.7 - 0.118 * tk) * f / (82.05736 * tk)
  expect_equal(rf_pco2 / rf_co2 - 1, term, tolerance = 1e-2)   # about 1.16e-6, to 1 percent
})
