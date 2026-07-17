## ---------------------------------------------------------------------------
## Regression tests for buffderiv(), the analytic partial derivatives of the
## Bolin sensitivity Be.
##
## These are cross-implementation tests, not golden-value tests.  They compare
## buffderiv()'s HAND-DERIVED derivatives against an INDEPENDENT forward-mode
## automatic-differentiation implementation (see helper-ad.R).  The two share no
## code: buffderiv takes K values from seacarb and hand-codes dK/dT and dK/dS in
## dlnK.R, while the AD route codes the K formulas and differentiates them
## automatically.  Agreement therefore validates BOTH.
##
## Why this matters: if a coefficient in seacarb's K1.R, Kb.R, Ksi.R ... is ever
## corrected, dlnK.R goes SILENTLY STALE.  No other test in the suite would
## catch it.  Test 2 below would fail immediately.
##
## skip_on_cran(): these are slow-ish and their purpose is developer-side drift
## detection, not user-facing correctness on CRAN's machines.
## ---------------------------------------------------------------------------

OPT <- list(k1k2 = "l", kf = "dg", ks = "d", pHscale = "T", b = "u74", warn = "n")

## ---------------------------------------------------------------------------
test_that("Kall() reproduces seacarb's own equilibrium constants", {
  skip_on_cran()

  ## This is the load-bearing test for the AD route.  Kall() is a second copy of
  ## seacarb's formulations; if it drifts from R/K1.R etc., every AD-based check
  ## below becomes meaningless.  Check the copy BEFORE trusting it.
  for (ts in list(c(-1, 34), c(5, 33.5), c(15, 30), c(25, 35), c(31, 37))) {
    T <- ts[1]; S <- ts[2]

    Ks_P0 <- Ks(S = S, T = T, P = 0, ks = "d", warn = "n")
    Kff   <- Kf(S = S, T = T, P = 0, pHscale = "F", kf = "dg", Ks_P0, Ks_P0, warn = "n")
    kc    <- kconv(S = S, T = T, P = 0, kf = "dg", Ks = Ks_P0, Kff = Kff, warn = "n")
    t2s   <- kc$ktotal2SWS
    s2c   <- kc$kSWS2total

    a <- Kall(T, S)
    b <- list(
      K1  = as.numeric(K1 (S=S,T=T,P=0,pHscale="T",k1k2="l",s2c,t2s,warn="n")),
      K2  = as.numeric(K2 (S=S,T=T,P=0,pHscale="T",k1k2="l",s2c,t2s,warn="n")),
      Kb  = as.numeric(Kb (S=S,T=T,P=0,pHscale="T",s2c,t2s,warn="n")),
      Kw  = as.numeric(Kw (S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
      Ksi = as.numeric(Ksi(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
      K1p = as.numeric(K1p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
      K2p = as.numeric(K2p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
      K3p = as.numeric(K3p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
      Ks  = as.numeric(Ks_P0),
      Kf  = as.numeric(Kff),
      BOR = as.numeric(bor(S = S, b = "u74")))

    for (k in names(b))
      expect_equal(a[[k]], b[[k]], tolerance = 1e-10,
                   info = sprintf("%s at T=%g S=%g. If this fails, seacarb's formulation
changed and helper-ad.R's Kall() must be updated to match.", k, T, S))
  }
})

## ---------------------------------------------------------------------------
test_that("buffderiv()'s analytic derivatives match automatic differentiation", {
  skip_on_cran()

  ## THE drift detector.  buffderiv's dlnK.R is hand-derived; AD is not.
  cases <- data.frame(
    T   = c(-1,   10,   20,   25,   29,    5),
    S   = c(34,   34.8, 35,   35,   36.5,  33.5),
    CT  = c(2150, 2050, 2000, 2000, 1950,  2200) * 1e-6,
    AT  = c(2300, 2290, 2300, 2300, 2380,  2320) * 1e-6,
    Pt  = c(1.5,  0.5,  0.2,  0.1,  0.05,  1.8)  * 1e-6,
    Sit = c(60,   10,   3,    2,    1,     70)   * 1e-6)

  for (i in seq_len(nrow(cases))) {
    with(cases[i, ], {
      ana <- do.call(buffderiv,
               c(list(15, AT, CT, S = S, T = T, Pt = Pt, Sit = Sit), OPT))
      h0  <- 10^(-do.call(carb,
               c(list(15, AT, CT, S = S, T = T, P = 0, Pt = Pt, Sit = Sit), OPT))$pH)
      ad  <- be_ad(T, S, CT, AT, Pt, Sit, h0)

      expect_equal(ana$Be,      ad$Be,      tolerance = 1e-9)
      expect_equal(ana$dBe_dT,  ad$dBe_dT,  tolerance = 1e-9)
      expect_equal(ana$dBe_dS,  ad$dBe_dS,  tolerance = 1e-9)
      expect_equal(ana$dBe_dCT, ad$dBe_dCT, tolerance = 1e-9)
      expect_equal(ana$dBe_dAT, ad$dBe_dAT, tolerance = 1e-9)
    })
  }
})

## ---------------------------------------------------------------------------
test_that("buffderiv() differentiates the alkalinity that carb() actually solves", {
  skip_on_cran()

  ## The diagnostic that catches a WRONG ALKALINITY (e.g. a missing -[HF] term, or
  ## silicate treated as diprotic when carb() is monoprotic).  It MUST be evaluated
  ## at carb()'s own unpolished h: with npolish > 0 the residual is driven to zero
  ## by construction and tells you nothing.
  set.seed(42); n <- 100
  S   <- runif(n, 28, 38)
  T   <- runif(n, -1.8, 32)
  AT  <- runif(n, 2000, 2450) * 1e-6
  CT  <- AT * runif(n, 0.83, 0.99)
  Pt  <- runif(n, 0, 2.5) * 1e-6
  Sit <- runif(n, 0, 120) * 1e-6

  raw <- do.call(buffderiv,
           c(list(15, AT, CT, S = S, T = T, Pt = Pt, Sit = Sit), OPT, npolish = 0))

  ## 1e-9 is loose enough for carb()'s own ~1e-10 convergence, tight enough that a
  ## missing fluoride term (which lands at 5e-7) fails loudly.
  expect_lt(max(abs(raw$alk_residual) / AT), 1e-9)
})

## ---------------------------------------------------------------------------
test_that("buffsun, buffzwg, buffesm and buffer agree on Be and Rf", {
  skip_on_cran()

  ## All four evaluate the same Sundquist et al. (1979) template with the same
  ## non-carbonate term w.  They must agree to floating point, since they take h
  ## from the same carb()/carbfull() call.  Exercises BOTH solver paths.
  cases <- list(
    list(T = 20, S = 35, Pt = 0,     Sit = 0,     NH4t = 0,    HSt = 0),
    list(T = -1, S = 34, Pt = 1.5e-6, Sit = 60e-6, NH4t = 0,    HSt = 0),
    list(T = 18, S = 35, Pt = 2.5e-6, Sit = 120e-6, NH4t = 3e-6, HSt = 2e-6))

  for (cs in cases) {
    a <- c(list(flag = 15, var1 = 2300e-6, var2 = 2000e-6, P = 0), cs, OPT)
    s <- do.call(buffsun, a)
    z <- do.call(buffzwg, a)
    e <- do.call(buffesm, a)
    u <- do.call(buffer,  a)

    expect_equal(s$Be, z$Be,    tolerance = 1e-9,  ignore_attr = TRUE)
    expect_equal(s$Rf, z$RF,    tolerance = 1e-12,  ignore_attr = TRUE)
    expect_equal(s$Rf, e$R,     tolerance = 1e-12,  ignore_attr = TRUE)
    expect_equal(s$Rf, u$BetaD, tolerance = 1e-12,  ignore_attr = TRUE)
  }
})
