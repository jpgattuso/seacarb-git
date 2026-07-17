## ---------------------------------------------------------------------------
## Forward-mode automatic differentiation by operator overloading.
## A dual number carries a value v and a vector d of partial derivatives with
## respect to 5 seed directions: (T, S, CT, AT, h).
## ---------------------------------------------------------------------------
ND <- 5L
# v is a numeric vector of length n; d is an n x ND matrix of partials.
dual    <- function(v, d) structure(list(v = v, d = d), class = "dual")
as.dual <- function(x, n) if (inherits(x, "dual")) x else
                          dual(rep(x, length.out = n), matrix(0, n, ND))
seed    <- function(v, i) { d <- matrix(0, length(v), ND); d[, i] <- 1; dual(v, d) }

Ops.dual <- function(e1, e2) {
  if (missing(e2)) {                                   # unary + and -
    if (.Generic == "-") return(dual(-e1$v, -e1$d))
    if (.Generic == "+") return(e1)
    stop("unsupported unary op: ", .Generic)
  }
  n <- max(if (inherits(e1,"dual")) NROW(e1$d) else 1L,
           if (inherits(e2,"dual")) NROW(e2$d) else 1L)
  a <- as.dual(e1, n); b <- as.dual(e2, n)
  switch(.Generic,
    "+" = dual(a$v + b$v, a$d + b$d),
    "-" = dual(a$v - b$v, a$d - b$d),
    "*" = dual(a$v * b$v, a$d * b$v + a$v * b$d),
    "/" = dual(a$v / b$v, (a$d * b$v - a$v * b$d) / b$v^2),
    "^" = { v <- a$v^b$v
            dual(v, v * (b$d * log(a$v) + b$v * a$d / a$v)) },
    stop("unsupported op: ", .Generic))
}

Math.dual <- function(x, ...) {
  switch(.Generic,
    "exp"  = dual(exp(x$v),  exp(x$v) * x$d),
    "log"  = dual(log(x$v),  x$d / x$v),
    "sqrt" = dual(sqrt(x$v), x$d / (2 * sqrt(x$v))),
    stop("unsupported math fn: ", .Generic))
}
val <- function(x) if (inherits(x, "dual")) x$v else x

## ---------------------------------------------------------------------------
## Kall(): a DELIBERATE, INDEPENDENT re-implementation of seacarb's constant set
## (k1k2="l", ks="d", kf="dg", b="u74", total pH scale).
##
## This duplication is the whole point.  buffderiv() takes K values FROM seacarb
## and hand-codes dK/dT and dK/dS in dlnK.R.  If a coefficient in Lueker (2000)
## or Millero (1995) is ever corrected in seacarb's R/*.R, dlnK.R will keep
## returning the OLD derivative and nothing else in the test suite would notice.
## Kall() is written straight from the published formulations, so AD through it
## follows any such change automatically, and test-buffderiv-ad.R then fails.
##
## Plain arithmetic only: no as.numeric(), no attr(), no branching on values,
## so dual numbers propagate.
## ---------------------------------------------------------------------------
Kall <- function(T, S) {
  TK <- T + 273.15 ; lnTK <- log(TK) ; sqS <- sqrt(S) ; ln10 <- log(10)
  Io <- 19.924*S/(1000 - 1.005*S) ; sqIo <- sqrt(Io)
  Cl <- S/1.80655 ; ST <- 0.14*Cl/96.062 ; FT <- 6.7e-5*Cl/18.9984
  lgS <- log(1 - 0.001005*S)
  K1 <- 10^(-3633.86/TK + 61.2172 - 9.67770*lnTK + 0.011555*S - 0.0001152*S*S)
  K2 <- 10^(-471.78/TK - 25.9290 + 3.16967*lnTK + 0.01781*S - 0.0001122*S*S)
  tb1 <- -8966.90 - 2890.53*sqS - 77.942*S + 1.728*S^1.5 - 0.0996*S*S
  Kb <- exp(tb1/TK + 148.0248 + 137.1942*sqS + 1.62142*S +
            (-24.4344 - 25.085*sqS - 0.2474*S)*lnTK + 0.053105*sqS*TK)
  Ks <- exp(-4276.1/TK + 141.328 - 23.093*lnTK +
            (-13856/TK + 324.57 - 47.986*lnTK)*sqIo +
            (35474/TK - 771.54 + 114.723*lnTK)*Io +
            (-2698/TK)*Io*sqIo + (1776/TK)*Io*Io + lgS)
  Kf <- exp(1590.2/TK - 12.641 + 1.525*sqIo + lgS)
  aph <- 1 + ST/Ks
  Sig <- aph/(aph + FT/Kf)                       # kSWS2total
  Kw <- exp(-13847.26/TK + 148.9802 - 23.6521*lnTK +
            (118.67/TK - 5.977 + 1.0495*lnTK)*sqS - 0.01615*S)*Sig
  Bsi <- -8904.2 - 458.79*sqIo + 188.74*Io - 12.1652*Io*Io
  Ksi <- exp(117.40 + 3.5913*sqIo - 1.5998*Io + 0.07871*Io*Io + Bsi/TK -
             19.334*lnTK + lgS)*Sig
  K1p <- exp(-4576.752/TK + 115.54 - 18.453*lnTK + (-106.736/TK+0.69171)*sqS +
             (-0.65643/TK-0.01844)*S)*Sig
  K2p <- exp(-8814.715/TK + 172.1033 - 27.927*lnTK + (-160.34/TK+1.3566)*sqS +
             (0.37335/TK-0.05778)*S)*Sig
  K3p <- exp(-3070.75/TK - 18.126 + (17.27039/TK+2.81197)*sqS +
             (-44.99486/TK-0.09984)*S)*Sig
  list(K1=K1,K2=K2,Kb=Kb,Kw=Kw,Ksi=Ksi,K1p=K1p,K2p=K2p,K3p=K3p,
       Ks=Ks,Kf=Kf,ST=ST,FT=FT,aph=aph,Ksa=Ks*aph,Kfa=Kf*aph,
       BOR=0.1284*S*1e-3/10.811)
}

## alkalinity residual F(h;T,S,CT,AT), exactly what carb() inverts
Fres <- function(h, T, S, CT, AT, Pt, Sit) {
  K <- Kall(T, S)
  Dc <- h^2 + K$K1*h + K$K1*K$K2
  Ac <- CT*(K$K1*h + 2*K$K1*K$K2)/Dc
  DP <- h^3 + K$K1p*h^2 + K$K1p*K$K2p*h + K$K1p*K$K2p*K$K3p
  NP <- K$K1p*K$K2p*h + 2*K$K1p*K$K2p*K$K3p - h^3
  Ac + K$BOR*K$Kb/(K$Kb+h) + K$Kw/h + Sit*K$Ksi/(K$Ksi+h) + Pt*NP/DP -
      h/K$aph - K$ST*h/(K$Ksa+h) - K$FT*h/(K$Kfa+h) - AT
}

## Be(h;T,S,CT) -- the explicit algebraic function buffsun() evaluates
Bfun <- function(h, T, S, CT, Pt, Sit) {
  K <- Kall(T, S)
  Dc <- h^2 + K$K1*h + K$K1*K$K2
  s  <- CT*h^2/Dc ; b <- CT*K$K1*h/Dc ; cc <- CT*K$K1*K$K2/Dc
  DP <- h^3 + K$K1p*h^2 + K$K1p*K$K2p*h + K$K1p*K$K2p*K$K3p
  NP <- K$K1p*K$K2p*h + 2*K$K1p*K$K2p*K$K3p - h^3
  DPh<- 3*h^2 + 2*K$K1p*h + K$K1p*K$K2p ; NPh <- K$K1p*K$K2p - 3*h^2
  w  <- 1/K$aph + K$ST*K$Ksa/(K$Ksa+h)^2 + K$FT*K$Kfa/(K$Kfa+h)^2 +
        K$Kw/h^2 + K$BOR*K$Kb/(K$Kb+h)^2 + Sit*K$Ksi/(K$Ksi+h)^2 +
        Pt*(NP*DPh - NPh*DP)/DP^2
  X <- h*w
  1 + cc/s + (X*b - 4*cc^2)/(s*(b + 4*cc + X))
}

## ---------------------------------------------------------------------------
## The four derivatives, entirely by AD + the implicit function theorem.
## NO hand-derived dK/dT, dK/dS, no G_h, no F_h.  AD supplies all of them.
## ---------------------------------------------------------------------------
be_ad <- function(T0, S0, CT0, AT0, Pt0, Sit0, h0) {
  nn <- length(T0)
  # 1. refine h so that F(h)=0 to machine precision (plain doubles, no AD needed)
  for (i in 1:4) {
    e <- h0*1e-7
    fp <- Fres(h0+e,T0,S0,CT0,AT0,Pt0,Sit0); fm <- Fres(h0-e,T0,S0,CT0,AT0,Pt0,Sit0)
    h0 <- h0 - Fres(h0,T0,S0,CT0,AT0,Pt0,Sit0)/((fp-fm)/(2*e))
  }
  # 2. ONE AD pass, seeding all five directions (T,S,CT,AT,h)
  Td<-seed(T0,1); Sd<-seed(S0,2); Cd<-seed(CT0,3); Ad<-seed(AT0,4); hd<-seed(h0,5)
  Fd <- Fres(hd, Td, Sd, Cd, Ad, Pt0, Sit0)
  Gd <- Bfun(hd, Td, Sd, Cd, Pt0, Sit0)
  # 3. implicit function theorem:  dBe/dY = G_Y - G_h * F_Y / F_h
  F_h <- Fd$d[,5] ; G_h <- Gd$d[,5]
  dBe <- Gd$d[,1:4,drop=FALSE] - (G_h/F_h) * Fd$d[,1:4,drop=FALSE]
  list(Be = Gd$v, dBe_dT = dBe[,1], dBe_dS = dBe[,2],
       dBe_dCT = dBe[,3], dBe_dAT = dBe[,4], h = h0)
}
