# tests/testthat

## helper-ad.R: the duplication is deliberate

`Kall()` in `helper-ad.R` re-implements seacarb's equilibrium constants. It looks like
code duplication. It is not a mistake, and it must not be removed.

**Why it exists.** `R/buffderiv.R` takes K *values* from seacarb and uses *hand-derived*
dK/dT and dK/dS from `R/dlnK.R`. If a coefficient in Lueker (2000) or Millero (1995) is
ever corrected in `R/K1.R`, `R/Kb.R`, `R/Ksi.R` ..., then `dlnK.R` keeps returning the OLD
derivative and **nothing else in this package would notice**. `Kall()` is written straight
from the published formulations, so automatic differentiation through it follows any such
change automatically, and `test-buffderiv-ad.R` fails immediately.

**If you refactor `Kall()` to call `K1()`, `K2()`, `Kb()` and so on, the test will still
pass, and it will have stopped testing anything.** That is the worst possible outcome: a
green tick that means nothing.

This has been verified to work. Injecting a 0.1 percent error into `dlnK.R`'s K1
temperature derivative produces 6 failures in `test-buffderiv-ad.R`, while every other
test in the suite still passes.

## Do not relax tolerances

If a test fails, find out why. Do not widen `tolerance=` to make it green. The tolerances
here are not arbitrary:

- `1e-9` on the analytic-vs-AD comparison: the two agree to ~1e-11 in practice, so 1e-9 is
  already two orders of slack.
- `1e-9` on the alkalinity residual at `npolish = 0`: loose enough for `carb()`'s own
  ~1e-10 convergence, tight enough that a **missing fluoride term** (which lands at 5e-7)
  fails loudly. That is exactly the bug this test was written to catch.

## The residual test must use npolish = 0

`buffderiv(..., npolish = 0)` uses `carb()`'s raw h. With `npolish > 0` the residual is
driven to zero by construction and the test would be vacuous.
