# seacarb

R package for seawater carbonate chemistry. Maintainer: James Orr (LSCE-IPSL).
This code is used in published science. Wrong numbers are worse than no numbers.

## The one rule that matters

**Never change a carbonate-chemistry formula on the basis of reading it.**

Every equilibrium constant, every buffer factor, every derivative in this package was
validated numerically. If you modify any of them, you must show the number, not the
argument. A formula that "looks wrong" or "looks simplifiable" is not evidence.

Before proposing any change to a formula, run the check and paste the output:

```bash
cd notebooks && Rscript validate_buffderiv.R && Rscript crossvalidate_buff.R
```

`notebooks/reference_values.csv` holds 10 reference points. Any change that moves
them is wrong until proved otherwise.

If you cannot validate a change numerically, do not make it. Say so instead.

## Do not "clean up" these files

- **`tests/testthat/helper-ad.R`** contains `Kall()`, which duplicates seacarb's own
  equilibrium constants. **This duplication is deliberate and is the entire point.** It is
  an independent second implementation, so that automatic differentiation through it can
  detect drift in the hand-derived derivatives in `R/dlnK.R`. Refactoring `Kall()` to call
  `K1()`, `K2()`, `Kb()` etc. would destroy the test while making it still pass. Do not
  do it. Do not suggest it.

- **`R/dlnK.R`** is hand-derived analytic algebra (dlnK/dT and dlnK/dS for every constant).
  It is verbose on purpose. Do not "simplify" it. Every line was checked against a finite
  difference of seacarb's own K functions.

- **`man/*.Rd`** are hand-written, not roxygen. **Never run `devtools::document()`** in
  this repo: it would overwrite them.

## Domain traps that have already caused bugs here

These are real, they are silent, and they have each cost a debugging session.

1. **`k1k2 = "x"` silently switches formulation.** `R/K1.R` does
   `is_outrange <- T>35 | T<2 | S<19 | S>43`, then switches to Waters et al. (2014).
   Below 2 degC you get a different formulation with no warning. Polar CMIP6 cells live
   exactly there. Always pass `k1k2 = "l"` explicitly in tests and examples.

2. **`kf = "x"` silently switches formulation.** `R/Kf.R` does
   `is_outrange <- T>33 | T<10 | S<10 | S>40`, Perez & Fraga in range, Dickson & Riley
   outside. This changes kSWS2total and hence Kw, Ksi, K1p, K2p, K3p by 0.5 to 0.9 percent.
   Always pass `kf = "dg"` explicitly.

3. **`carb()` and `carbfull()` do not solve the same alkalinity.**
   `carb()` (via `calculate_carb`, `fullresult=FALSE`) treats silicate as **monoprotic**
   (K1si only) and cannot take NH4t or HSt. `carbfull()` is **diprotic** and takes both.
   Any buffer factor must use the alkalinity of the solver actually called. Getting this
   wrong is a real bug that shipped for years.

4. **Total alkalinity includes `-[HF]`.** SolveSAPHE carries it explicitly. Any expression
   for `-dAlk/dh` must include `FT*Kfa/(Kfa+h)^2` with `Kfa = Kf_free * (1 + ST/Ks)`, i.e.
   Kf on the **same pH scale as h**. Using free-scale Kf with a total-scale h is wrong by
   28 percent.

5. **`buffer()`'s factors are per mole of the substance ADDED**, with strong acid meaning
   dALK = -1, dDIC = 0. So `PhiH` is correctly negative. The Beta and Pi factors use
   **fCO2, not pCO2** (0.34 percent apart at 20 degC), and Pi is per **umol/kg**. The code
   is right; the docs were wrong until 3.3.5. Do not "fix" the code to match old docs.

## How to verify any change

```r
devtools::load_all()
devtools::test()          # expect ~98 assertions, 0 failures
```
```bash
R CMD build .
R CMD check --as-cran seacarb_*.tar.gz    # expect 0 errors, 0 warnings
```

A pre-existing NOTE about an invalid IAEA URL in `man/oa.Rd` is not your problem.

## Layout

- `R/` package code. `dlnK.R` and `calculate_carb.R` are internal (not in NAMESPACE).
- `man/` hand-written Rd files.
- `tests/testthat/` the regression suite. Shipped to CRAN but `skip_on_cran()`,
  so it runs locally and in GitHub Actions, not on CRAN's machines.
- `notebooks/` NOT part of the package (see `.Rbuildignore`). Validation scripts,
  notebooks, and `reference_values.csv`.
- `.Rbuildignore` every line is a Perl regex. **Comments are not allowed** and will break
  `R CMD build`.

## Working style here

Prefer showing a number over making an argument. If asked to explain a discrepancy,
measure it first. If a test fails, do not adjust the tolerance to make it pass; find out
why. State uncertainty explicitly rather than asserting with false confidence.
