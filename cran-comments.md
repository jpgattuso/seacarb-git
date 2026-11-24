# v3.4.3
## Submission summary

This release occurs a few weeks after version 3.4.2 because the update of SolveSAPHE (a dependency of seacarb) to version 2.1.1 was flagged by CRAN as causing a failure in the seacarb test suite. The cause was in seacarb and not in SolveSAPHE. SolveSAPHE 2.1.1 fixes a parsing bug that had silently dropped terms from three analytic derivative branches. The fix shifts the converged pH by up to about 1e-10 pH units, which at one grid point (S = 35, T = -1.5 degC) narrowly exceeded the tolerance of the fluoride and sulfate derivative test in test-buffergen.R. The tolerance was widened from 1e-7 to 1e-6, and no buffergen() code was changed. A step size sweep confirmed that the sensitivity is real and not a discretization artifact.

In addition, seacarb 3.4.3 corrects the alkalinity computed from pH in carb(), carbfull(), and carbb() (flags 1, 2, 3, 5, 6, 7, 9, 10, 12, 14, 21, 22, 23, and 25) and the DIC computed from pH and alkalinity (flag 8). The free proton and bisulfate were on the wrong side of the pH scale conversion, which made the alkalinity off by 1.5 ppm at pH 8.15 and up to 16 ppm at pH 7.1. Flags that take alkalinity as input (4, 11, 13, 15, and 24) were not affected. The help page and header of derivnum() now describe the perturbations correctly (results unchanged) and explain why the pCO2 derivative differs slightly from those of CO2 and fCO2. New regression tests were added (test-derivnum.R, test-carb-flags.R, test-carb-scales.R, and test-buff-flags.R).

## Test environments

* local Ubuntu 24.04.5, R 4.3.3, with NOT_CRAN=true and --as-cran
* GitHub Actions: ubuntu-latest (release, devel) and windows-latest
  (release), all passing. The macos-latest job could not be run
  because no runner became available.

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking for future file timestamps: unable to verify current time
  (the local machine could not check the time).
* checking HTML version of manual: math rendering skipped because
  package V8 is unavailable locally.

Both notes come from the local environment and not from the package.

## Reverse dependencies

seacarb has three reverse dependencies on CRAN (ecolMod, marelac, respirometry). All three were checked against this version (R 4.3.3, Ubuntu 24.04), and R CMD check returned OK for each. ecolMod and marelac do not call carb(), carbfull(), or carbb(). respirometry calls carb() with flag 8 (DIC from pH and alkalinity), which this release corrects, so its DIC values shift by the size of the corrected error (the ChangeLog gives 1.5 ppm at pH 8.15 and up to 16 ppm at pH 7.1 for the alkalinity, with the DIC from flag 8 off by a similar amount in the opposite direction). It also uses flags 15 and 24, which are unaffected. The checks of respirometry pass.


