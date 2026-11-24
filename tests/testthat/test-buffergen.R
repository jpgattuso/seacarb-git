## ---------------------------------------------------------------------------
## Regression tests for buffergen() (Hagens and Middelburg, 2016).
##
## What these tests protect against
##   1. A wrong sign in the derivative terms of the fluoride and sulfate systems
##      (the 3.4.1 code had both reversed), which made RF wrong by 5 to 20 ppm.
##   2. Counting HSO4 (or HF) as a separate alkalinity term on a pH scale where it
##      is already part of the definition of [H+].
##   3. Returning plausible looking numbers for derivatives that are not defined on
##      a given pH scale. Those outputs must be NA.
##
## How the tests are independent of buffergen()
##   Test 4 does not call any buffergen formula. It rebuilds the total alkalinity
##   of each pH scale from the species concentrations that carbfull() returns, and
##   obtains every derivative by NUMERICAL partial derivatives of that alkalinity
##   (implicit function theorem). Only the equilibrium constants are taken from
##   carbfull()'s own species (K = [H+][base]/[acid]), so a bug in the derivative
##   algebra of buffergen() cannot hide behind a copy of the same algebra.
##
## Alkalinity on each pH scale (what carbfull() solves; same as SolveSAPHE):
##   free  (F):  ... + [OH-] - [H+]F - [HSO4-] - [HF]
##   total (T):  ... + [OH-] - [H+]T - [HF]      ([HSO4-] is inside [H+]T)
##   SWS:        ... + [OH-] - [H+]SWS           ([HSO4-] and [HF] are inside [H+]SWS)
##
## Settings: k1k2 and kf are always passed explicitly. Their defaults ("x") switch
## formulation silently outside certain ranges (see CLAUDE.md).
## Tolerances are not to be widened to make a failing test pass without numerical
## justification recorded here (see the fluoride/sulfate test below for a documented
## exception).
## ---------------------------------------------------------------------------

## The helpers (GRID, bg(), indep(), d1(), relerr()) are in helper-buffergen.R.

## ---------------------------------------------------------------------------
test_that("RF agrees with buffsun() on the total scale, same acid-base systems", {
  ## buffsun() and buffergen() share no derivative code. With Sit = 0 both solve
  ## the same alkalinity. With silicate, carb() is monoprotic and carbfull() is
  ## diprotic, so a tiny NH4t (1e-20 mol/kg) is added to force buffsun() onto
  ## carbfull() too. That changes nothing else.
  for (i in seq_len(nrow(GRID))) {
    g <- GRID[i, ]
    rf <- bg(i, "T")$RF
    ref <- buffsun(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T, P = g$P,
      Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt, k1k2 = OPT$k1k2,
      kf = OPT$kf, ks = OPT$ks, pHscale = "T", b = OPT$b)$Rf
    expect_equal(as.numeric(rf), as.numeric(ref), tolerance = 1e-9)
  }
})

test_that("RF at the reference state is unchanged", {
  ## S=35, T=18, ALK=2300, DIC=2000 umol/kg, no nutrients. The 3.4.1 code gave
  ## 9.6914383246 here, which is 6.0e-6 too high.
  g <- suppressWarnings(buffergen(flag = 15, var1 = 2300e-6, var2 = 2000e-6, S = 35, T = 18,
    k1k2 = "l", kf = "pf", ks = "d", pHscale = "T", b = "u74"))
  expect_equal(as.numeric(g$RF), 9.6913805691, tolerance = 1e-9)
})

test_that("RF does not depend on the pH scale", {
  ## RF is dimensionless: d ln(s) / d ln(DIC) at constant alkalinity.
  ## Uses kf = "pf". With kf = "dg", carbfull() itself returns a slightly different
  ## physical state on the SWS scale (a finite difference of carbfull() agrees with
  ## buffergen() on every scale to 2e-9, but differs between SWS and T by up to 1.6e-6
  ## in RF at some states). That is a property of the constants, not of buffergen().
  for (i in seq_len(nrow(GRID))) {
    r <- sapply(c("F", "T", "SWS"), function(sc) as.numeric(bg(i, sc, kf = "pf")$RF))
    expect_lt(diff(range(r)) / r[["T"]], 1e-8)
  }
})

test_that("a vector of pH scales gives the same result as separate calls", {
  sc <- c("F", "T", "SWS")
  v <- suppressWarnings(buffergen(flag = 15, var1 = rep(2300e-6, 3), var2 = rep(2000e-6, 3),
    S = 35, T = 18, k1k2 = "l", kf = "dg", ks = "d", pHscale = sc, b = "u74"))
  one <- sapply(sc, function(s) as.numeric(suppressWarnings(buffergen(flag = 15, var1 = 2300e-6,
    var2 = 2000e-6, S = 35, T = 18, k1k2 = "l", kf = "dg", ks = "d", pHscale = s, b = "u74"))$RF))
  expect_equal(as.numeric(v$RF), as.numeric(one), tolerance = 1e-12)
})

test_that("fluoride and sulfate derivatives match the alkalinity definition where defined", {
  ## Defined: free scale, both systems; total scale, fluoride system only.
  for (i in seq_len(nrow(GRID))) {
    for (sc in c("F", "T")) {
      ind <- indep(i, sc)
      sp <- if (sc == "F") c(FLUORIDE, SULFATE) else FLUORIDE
      for (n in sp) for (q in OUTS8) {
        ## Tolerance widened from 1e-7 to 1e-6 after SolveSAPHE 2.1.1 (2026-11): that
        ## release fixed an R continuation-line parsing bug in anw.R (dropped terms in
        ## the fixed-CO2/HCO3/CO3 analytic derivative branches), which nudges the
        ## converged H by ~2e-10 relative (within the solver's own convergence
        ## tolerance). buffergen()'s own analytic dX.dALK for sulfate at GRID row 3
        ## (S=35, T=-1.5 degC) is essentially unchanged by that shift (~3e-10 relative),
        ## but the true value of dX.dALK at that cold/polar state is itself sensitive
        ## enough to H that the same shift in H moves it by ~1.1e-7 relative. Confirmed
        ## via a finite-difference step-size sweep (1e-2 to 1e-5) that this is a real,
        ## step-independent sensitivity, not a discretization artifact, and that no
        ## other grid row-scale-species-quantity is affected. After correspondence with
        ## J.-M. Epitalon and G. Munhoven, following SolveSAPHE bug fix from
        ## Otmane Bentaleb (circa Nov 2026).
        expect_equal(as.numeric(ind$g[[q]][n]), ind$truth[[n]][[q]], tolerance = 1e-6)
      }
    }
  }
})

test_that("derivatives with respect to a system that is part of the pH scale are NA", {
  ## Not separate alkalinity terms: sulfate on T and SWS, fluoride on SWS.
  for (i in seq_len(nrow(GRID))) {
    for (sc in c("F", "T", "SWS")) {
      g <- bg(i, sc)
      na.sys <- c(if (sc != "F") SULFATE, if (sc == "SWS") FLUORIDE)
      ok.sys <- setdiff(c(FLUORIDE, SULFATE), na.sys)
      for (q in OUTS12) {
        for (n in na.sys) expect_true(is.na(g[[q]][n]))
        for (n in ok.sys) expect_false(is.na(g[[q]][n]))
      }
      ## RF and beta.H never depend on this
      expect_false(is.na(g$RF))
      expect_false(is.na(g$beta.H))
    }
  }
})

test_that("dH.dALK is the same for every species where it is defined", {
  ## The proton response to an alkalinity change at fixed totals cannot depend on
  ## which species of which system is used to describe it. With the 3.4.1 code the
  ## fluoride and sulfate rows differed from the others by about 5e-6 relative.
  for (i in seq_len(nrow(GRID))) {
    g <- bg(i, "F")
    v <- g$dH.dALK
    v <- v[!is.na(v)]
    expect_lt(max(abs(v / v[["CO2"]] - 1)), 1e-8)
  }
})
