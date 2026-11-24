# seacarb 3.4.3 (2026-11-24)

* Fix: `carb()`, `carbfull()`, and `carbb()` computed the alkalinity from pH (flags 1, 2, 3, 5, 6, 7, 9, 10, 12, 14, 21, 22, 23, and 25) and the DIC from pH and alkalinity (flag 8) with the free proton and bisulfate on the wrong side of the pH scale conversion. The alkalinity was off by 1.5 ppm at pH 8.15 and up to 16 ppm at pH 7.1 (about 10 percent more on the seawater scale). The DIC from flag 8 was off by a similar amount in the opposite direction. The hydrogen fluoride term also used the total-scale proton on the free and seawater scales (up to about 0.1 ppm). Flags 4, 11, 13, 15, and 24, which take ALK as input, were not affected.
* Update: `derivnum()` and `derivnum.Rd` now describe the perturbations correctly. Each perturbation is a relative delta times a fixed reference value, so it is the same absolute amount whatever the input is. Results are unchanged.
* Update: `derivnum.Rd` and the header of `derivnum()` now explain that the derivative of pCO2 is not exactly proportional to the derivatives of fCO2 and CO2. The Revelle factor computed from the pCO2 derivative is larger by about 1 ppm (up to about 4 ppm at high CO2) than the one computed from the derivative of CO2 or fCO2, which agrees with `buffsun()` to 0.01 ppm or less.
* New tests: `test-derivnum.R`, `test-carb-flags.R`, `test-carb-scales.R`, and `test-buff-flags.R`.
* Update: the tolerance in `test-buffergen.R` for the fluoride and sulfate derivative test was widened from 1e-7 to 1e-6, needed after SolveSAPHE 2.1.1 fixed a parsing bug. No `buffergen()` code was changed.

# seacarb 3.4.2 (2026-11-06)

* Fix: `buffergen()` had reversed signs in the derivative terms of the fluoride and sulfate systems, and counted sulfate twice on the total pH scale. The Revelle factor was too high or too low by 5 to 20 ppm. Outputs that are not defined on a pH scale are now `NA`.
* Fix: `buffsun()`, `buffzwg()`, `buffesm()`, and `buffer()` always evaluated the fluoride and sulfate terms of the noncarbonate alkalinity for the total pH scale. On the free and seawater scales, Rf was off by up to about 13 ppm. Results on the total scale (the default) are unchanged.
* Fix: the ammonium and sulfide constants `Kn` and `Khs` were converted to the chosen pH scale with a different Kf than all other constants whenever `kf` was set explicitly to "dg". With nonzero NH4t or HSt, pCO2 from `carbfull()` differed between pH scales by up to about 4e-6 relative, and Rf by up to about 1.6 ppm on the seawater scale only. The default `kf = "x"` is unchanged.
* Update: the help pages of 21 functions now state that Kf from both Perez and Fraga (1987) and Dickson and Riley (1979) is on the free scale. The code was already correct since 3.4.1.

# seacarb 3.4.1 (2026-07-30)

* Fix: `Kf()` now treats the equation of Perez and Fraga (1987), `kf = "pf"`, as being on the free hydrogen ion scale, and the unwarranted total-to-free conversion was removed. The value of Kf changes by up to about 34 percent (largest at high temperature) for `kf = "pf"` on the total or seawater scale, but the effect on carbonate-system results is much smaller. For example, the Revelle factor typically changes only in its fourth significant figure. The Dickson and Riley (1979) option (`kf = "dg"`) is unaffected.
* Update: the unused `Ks_p0` argument was removed from `Kf()`.

# seacarb 3.4.0 (2026-07-17)

Major revision of the buffer-factor calculations, with a unified equation for the Revelle factor Rf and the Bolin sensitivity Be (the partial derivative dDIC/d[CO2*]) in the form of Sundquist et al. (1979).

* New: `buffsun()` computes Rf and Be from the equation of Sundquist et al. (1979).
* New: `buffzwg()` computes Be and Rf following Zeebe and Wolf-Gladrow (2001), with the typo in their Eq. 1.5.88 corrected.
* New: `buffderiv()` computes the four analytic partial derivatives of the Bolin sensitivity with respect to T, S, DIC, and ALK. It is validated against finite differences through `carb()` and against an independent automatic-differentiation implementation.
* Fix: the noncarbonate alkalinity used by the buffer functions is now complete (borate, water, silicate, phosphate, ammonium, hydrogen sulfide, and fluoride).
* Fix: silicate is treated consistently with the solver in the buffer functions (monoprotic in `carb()`, diprotic in `carbfull()`).
* Fix: in `buffer()`, PhiB and BetaB were in error by up to 0.16 percent, and PhiC and BetaC by up to 1.5 percent, in waters rich in phosphate and silicate.
* Update: `Kn()` and `Khs()` accept the scale-conversion arguments already accepted by `K1()`, `K2()`, and the other constants. Defaults are unchanged, so existing results are identical.
* Documentation: the `buffer()` help now states that PhiH and PiH are per mole of strong acid added (and so negative), that the Beta and Pi factors use fCO2 and not pCO2, and that the Pi factors are per umol/kg. The code was already right.
* New: a regression test suite under `tests/testthat`, marked `skip_on_cran()`.

# seacarb 3.3.4 (2026-04-16)

* James Orr takes over as maintainer from Jean-Pierre Gattuso.

# seacarb 3.3.3 (2024-02-14)

* New: K1 and K2 from Cai and Wang (1998).
* New: `pHnbs2sws()`, `pHsws2nbs()`, and `fH()`.

# seacarb 3.3.2 (2023-05-21)

* Fix: bugs in `carb()`, `calculate_carb()`, and `derivnum()` identified by S. Petton.

# seacarb 3.3.0 (2021-10-15)

* Update: implementation of the Solver Suite for Alkalinity-PH Equations (SolveSAPHE) of Munhoven (2013, 2021), which calculates pH from total alkalinity and DIC.

Older entries are in the `ChangeLog` file.
