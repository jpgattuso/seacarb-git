## ---------------------------------------------------------------------------
## buffsun(), buffzwg(), buffesm(), and buffer(): Rf must not depend on the input pair.
##
## What this protects against
##   The Revelle factor is a property of the state, not of the pair of variables used
##   to describe it. The four routines call carb() or carbfull() with the flag given,
##   so any inconsistency between flags in calculate_carb() reaches Rf. Until seacarb
##   3.4.2, flag 8 (pH and ALK) returned a DIC that was wrong by 1.5 ppm at pH 8.15 on
##   the total scale (see test-carb-flags.R), and Rf from flag 8 was off by 1.7e-6
##   relative on this grid. Nothing in the suite noticed, because the other buffer tests
##   use flag 15.
##
## How the test works
##   For each state of GRID (helper-buffergen.R) and each pH scale, Rf is computed with
##   flag 15 (ALK and DIC) and with the other 19 flags, each given the matching pair of
##   reference values from carbfull(). Flags 21 to 25 (pCO2 input) are only used where
##   P = 0. A tiny NH4t (1e-20 mol/kg) is added so that every routine uses carbfull()
##   (diprotic silicate), as in test-buff-scales.R. Rf from each flag is compared with
##   Rf from flag 15 for the same routine, so no other implementation is involved.
##
## Tolerance
##   1e-7 relative. With the fix the largest difference on this grid is 3.2e-9, from flags
##   21 to 25 (pCO2 input), where the pCO2 given under the default gas = "potential" is
##   not returned exactly (observed, not investigated; these routines have no gas
##   argument). All other flags agree to 1e-13. The unfixed flag 8 deviates by 1.7e-6.
##   Tolerances are not to be widened to make a failing test pass.
## ---------------------------------------------------------------------------

BUFF_PAIRS <- list(
  `1`  = c("pH", "CO2"),   `2`  = c("CO2", "HCO3"), `3`  = c("CO2", "CO3"),
  `4`  = c("CO2", "ALK"),  `5`  = c("CO2", "DIC"),  `6`  = c("pH", "HCO3"),
  `7`  = c("pH", "CO3"),   `8`  = c("pH", "ALK"),   `9`  = c("pH", "DIC"),
  `10` = c("HCO3", "CO3"), `11` = c("HCO3", "ALK"), `12` = c("HCO3", "DIC"),
  `13` = c("CO3", "ALK"),  `14` = c("CO3", "DIC"),  `15` = c("ALK", "DIC"),
  `21` = c("pCO2", "pH"),  `22` = c("pCO2", "HCO3"), `23` = c("pCO2", "CO3"),
  `24` = c("pCO2", "ALK"), `25` = c("pCO2", "DIC"))

## Rf from the four routines for one state, one pH scale, and one input pair
rf_by_flag <- function(flag, var1, var2, g, sc) {
  a <- list(flag = flag, var1 = var1, var2 = var2, S = g$S, T = g$T, P = g$P,
            Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt,
            k1k2 = OPT$k1k2, kf = OPT$kf, ks = OPT$ks, pHscale = sc, b = OPT$b)
  c(buffsun = as.numeric(do.call(buffsun, c(a, warn = "n"))$Rf),
    buffzwg = as.numeric(do.call(buffzwg, c(a, warn = "n"))$RF),
    buffesm = as.numeric(do.call(buffesm, c(a, warn = "n"))$R),
    buffer  = as.numeric(do.call(buffer,  c(a, warn = "n"))$BetaD))
}

## Flags for which any routine's Rf differs from its flag-15 value by 'tol' or more
buff_flags_off <- function(i, sc, tol = 1e-7) {
  g <- GRID[i, ]
  ref <- suppressWarnings(carbfull(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T,
    P = g$P, Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt,
    k1k2 = OPT$k1k2, kf = OPT$kf, ks = OPT$ks, pHscale = sc, b = OPT$b))
  r0 <- suppressWarnings(rf_by_flag(15, g$ALK, g$DIC, g, sc))
  bad <- character()
  for (f in names(BUFF_PAIRS)) {
    fl <- as.numeric(f)
    if (fl == 15 || (fl >= 21 && g$P > 0)) next
    v <- vapply(BUFF_PAIRS[[f]], function(n)
      switch(n, ALK = g$ALK, DIC = g$DIC, as.numeric(ref[[n]])), numeric(1))
    r <- suppressWarnings(rf_by_flag(fl, v[1], v[2], g, sc))
    dev <- abs(r / r0 - 1)
    if (any(!is.finite(dev)) || max(dev) >= tol) {
      k <- which.max(replace(dev, !is.finite(dev), Inf))
      bad <- c(bad, sprintf("flag %s: %s %.2e", f, names(dev)[k], dev[[k]]))
    }
  }
  bad
}

for (sc in c("T", "F", "SWS")) {
  test_that(sprintf("Rf from buffsun, buffzwg, buffesm, and buffer is the same for all 20 flags on the %s scale", sc), {
    for (i in seq_len(nrow(GRID))) {
      bad <- buff_flags_off(i, sc)
      expect_true(length(bad) == 0,
                  info = sprintf("GRID row %d: %s", i, paste(bad, collapse = "; ")))
    }
  })
}
