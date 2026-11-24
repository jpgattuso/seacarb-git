## ---------------------------------------------------------------------------
## carb() and carbfull(): every input pair (flag) must give the same state.
##
## What this protects against
##   Until seacarb 3.4.2, the flags that obtain the alkalinity from pH (1, 2, 3, 5, 6,
##   7, 9, 10, 12, 14, 21, 22, 23, 25) and flag 8 (which obtains DIC from pH and ALK)
##   computed the total alkalinity (or, for flag 8, the DIC) with the free proton and
##   bisulfate terms on the wrong side of the pH scale conversion: [H+]F was taken as
##   [H+]T * (1 + ST/Ks) instead of [H+]T / (1 + ST/Ks), and likewise on the seawater
##   scale. The HF term also used the wrong proton variable on F and SWS. The result
##   was an alkalinity off by 1.5 ppm at pH 8.15 on the total scale, 5.6 ppm at pH
##   7.57, and 16 ppm at pH 7.1 (SWS: about 10 percent more). Nothing in the suite
##   noticed, because every other test uses flag 15 (ALK and DIC), which solves for
##   pH with the correct alkalinity.
##
## How the test works
##   A reference state is computed with flag 15 (ALK and DIC given). Each of the 20
##   flags is then called with the matching pair of reference values, and every
##   output that exists for all flags is compared with the reference. The inputs
##   are not rounded, which is not the case for example 4 in carb.Rd (the only
##   earlier all-flag check, which has no assertions and uses the total scale only).
##
## Tolerance
##   1e-7 relative. The unfixed code deviates by 1.5e-6 (total scale, pH 8.15) to
##   3e-6 (SWS), so 1e-7 fails loudly. The fixed code deviates by at most 6.5e-9 on
##   this grid. That residual comes from flags 21 to 25 (pCO2 input), where pCO2 and
##   fCO2 are converted with a fugacity factor evaluated at slightly different
##   arguments in the two directions (observed, not investigated). Tolerances are not
##   to be widened to make a failing test pass.
##
## Notes
##   * GRID, from helper-buffergen.R, has 9 states with and without nutrients, with
##     pressure, and cold and warm water. Flags 21 to 25 are only used where P = 0.
##   * For flags 21 to 25, gas = "standard" is used. At P = 0 it round trips exactly.
##     The other two options differ at P = 0 (about 4e-9 relative for "potential",
##     about 1e-6 for "insitu"), although the help page says that all three should
##     agree at surface pressure. That is a separate matter and is not tested here.
##   * k1k2, kf, ks, and b are always passed explicitly (see CLAUDE.md).
## ---------------------------------------------------------------------------

FLAG_PAIRS <- list(
  `1`  = c("pH", "CO2"),   `2`  = c("CO2", "HCO3"), `3`  = c("CO2", "CO3"),
  `4`  = c("CO2", "ALK"),  `5`  = c("CO2", "DIC"),  `6`  = c("pH", "HCO3"),
  `7`  = c("pH", "CO3"),   `8`  = c("pH", "ALK"),   `9`  = c("pH", "DIC"),
  `10` = c("HCO3", "CO3"), `11` = c("HCO3", "ALK"), `12` = c("HCO3", "DIC"),
  `13` = c("CO3", "ALK"),  `14` = c("CO3", "DIC"),  `15` = c("ALK", "DIC"),
  `21` = c("pCO2", "pH"),  `22` = c("pCO2", "HCO3"), `23` = c("pCO2", "CO3"),
  `24` = c("pCO2", "ALK"), `25` = c("pCO2", "DIC"))

FLAG_OUT <- c("pH", "CO2", "fCO2", "HCO3", "CO3", "DIC", "ALK",
              "OmegaAragonite", "OmegaCalcite")

## Flags whose result differs from the flag-15 reference by more than 'tol' (relative)
## for one state and one pH scale. Returns a character vector, empty when all agree.
flags_off <- function(fun, i, sc, tol = 1e-7) {
  g <- GRID[i, ]
  base <- list(S = g$S, T = g$T, P = g$P, Pt = g$Pt, Sit = g$Sit,
               k1k2 = OPT$k1k2, kf = OPT$kf, ks = OPT$ks, b = OPT$b,
               pHscale = sc)
  if (identical(fun, carbfull)) base <- c(base, list(NH4t = g$NH4, HSt = g$HSt))
  else base <- c(base, list(warn = "n"))          # carbfull() has no warn argument
  ref <- suppressWarnings(do.call(fun, c(list(flag = 15, var1 = g$ALK, var2 = g$DIC), base)))
  bad <- character()
  for (f in names(FLAG_PAIRS)) {
    fl <- as.numeric(f)
    if (fl >= 21 && g$P > 0) next
    v <- vapply(FLAG_PAIRS[[f]], function(n)
      switch(n, ALK = g$ALK, DIC = g$DIC, as.numeric(ref[[n]])), numeric(1))
    a <- if (fl >= 21) c(base, list(gas = "standard")) else base
    r <- suppressWarnings(do.call(fun, c(list(flag = fl, var1 = v[1], var2 = v[2]), a)))
    dev <- vapply(FLAG_OUT, function(q) abs(as.numeric(r[[q]]) / as.numeric(ref[[q]]) - 1), numeric(1))
    if (any(!is.finite(dev)) || max(dev) >= tol) {
      q <- names(dev)[which.max(replace(dev, !is.finite(dev), Inf))]
      bad <- c(bad, sprintf("flag %s (%s): %s %.2e", f, sc, q, dev[[q]]))
    }
  }
  bad
}

for (fn in c("carb", "carbfull")) {
  for (sc in c("T", "F", "SWS")) {
    test_that(sprintf("%s: all 20 flags agree with flag 15 on the %s scale", fn, sc), {
      fun <- get(fn)
      for (i in seq_len(nrow(GRID))) {
        bad <- flags_off(fun, i, sc)
        expect_true(length(bad) == 0,
                    info = sprintf("GRID row %d: %s", i, paste(bad, collapse = "; ")))
      }
    })
  }
}
