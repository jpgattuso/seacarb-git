## ---------------------------------------------------------------------------
## carb() and carbfull(): the physical state must not depend on the pH scale.
##
## What this protects against
##   The pH scale (total, free, or seawater) only changes how [H+] is expressed, so
##   with the same ALK and DIC the concentrations of CO2, HCO3, CO3, and everything
##   derived from them (fCO2, pCO2, the saturation states) cannot depend on it. Every
##   equilibrium constant is converted to the chosen scale, so one constant converted
##   with a different Kf, or not converted at all, breaks this. Past examples are the
##   free-to-total conversion applied to Perez and Fraga (1987) Kf, fixed in 3.4.1, and
##   Kn and Khs converted with the default kf = "x" instead of the kf requested, fixed
##   in 3.4.2 (an inconsistency of up to about 0.9 percent in those two constants, see
##   the comments in Kn.R).
##
## How the test works
##   Each state of GRID (helper-buffergen.R) is solved with flag 15 (ALK and DIC) on the
##   total, free, and seawater scales. Outputs that do not depend on the scale are
##   compared with the total-scale result. carbfull() is given the nutrients of GRID,
##   including ammonium and sulfide, so that Kn and Khs are used.
##
## Tolerance
##   1e-9 relative. Without any scale error the spread is 2.2e-14 over GRID, which
##   includes pressures up to 3000 dbar, so 1e-9 is four orders of slack for round-off.
##   A constant converted with the wrong Kf is far larger (see above). Tolerances are not
##   to be widened to make a failing test pass.
##
## Notes
##   * Only k1k2 = "l" (Lueker et al., 2000) is tested. Millero (2010) and Waters et al.
##     (2014) provide K1 and K2 as separate fits for each pH scale, so their results
##     differ between scales by construction (about 3 percent and 0.1 percent in CO2).
##   * kf = "dg" and "pf" are both tested, and always passed explicitly (see CLAUDE.md).
## ---------------------------------------------------------------------------

SCALE_OUT <- c("CO2", "fCO2", "pCO2", "HCO3", "CO3", "OmegaAragonite", "OmegaCalcite")

## Largest relative difference from the total-scale result, over the scale-independent
## outputs and over the free and seawater scales, for one state of GRID.
scale_spread <- function(fun, i, kf) {
  g <- GRID[i, ]
  one <- function(sc) {
    a <- list(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T, P = g$P,
              Pt = g$Pt, Sit = g$Sit, k1k2 = "l", kf = kf, ks = OPT$ks, b = OPT$b,
              pHscale = sc)
    a <- if (identical(fun, carbfull)) c(a, list(NH4t = g$NH4, HSt = g$HSt))
         else c(a, list(warn = "n"))
    suppressWarnings(do.call(fun, a))
  }
  r <- lapply(c("T", "F", "SWS"), one)
  max(vapply(SCALE_OUT, function(q)
    max(abs(vapply(r, function(x) as.numeric(x[[q]]), numeric(1)) /
            as.numeric(r[[1]][[q]]) - 1)), numeric(1)))
}

for (fn in c("carb", "carbfull")) {
  for (kf in c("dg", "pf")) {
    test_that(sprintf("%s with kf = \"%s\": the state is the same on T, F, and SWS", fn, kf), {
      fun <- get(fn)
      for (i in seq_len(nrow(GRID))) {
        expect_lt(scale_spread(fun, i, kf), 1e-9,
                  label = sprintf("spread across pH scales, GRID row %d", i))
      }
    })
  }
}
