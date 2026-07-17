source("setup.R")   # library(seacarb) + the buffer functions from ../R
# Validation of buffderiv() against seacarb.
#
# buffderiv() now takes every K VALUE from seacarb's own functions, so there is
# nothing left to check there: it is seacarb's constants by construction.
# What must still be tested is (1) the hand-coded dK/dT and dK/dS in dlnK.R,
# (2) that our Anc is the alkalinity carb() actually solves, and (3) the four
# derivatives themselves.
## Richardson-extrapolated central difference, O(eps^4)
rich <- function(f, x, eps) {
  d1 <- (f(x+eps)   - f(x-eps))   / (2*eps)
  d2 <- (f(x+2*eps) - f(x-2*eps)) / (4*eps)
  (4*d1 - d2)/3
}

cat("seacarb", as.character(packageVersion("seacarb")), "/", R.version.string, "\n\n")

## seacarb's K's, fetched exactly as buffderiv() (hence calculate_carb()) fetches them
scK <- function(S, T) {
  Ks_P0  <- Ks(S=S, T=T, P=0, ks="d", warn="n")
  Kff_P0 <- Kf(S=S, T=T, P=0, pHscale="F", kf="dg", Ks_P0, Ks_P0, warn="n")
  t2s    <- kconv(S=S,T=T,P=0,kf="dg",Ks=Ks_P0,Kff=Kff_P0,warn="n")$ktotal2SWS
  s2c    <- kconv(S=S,T=T,P=0,kf="dg",Ks=Ks_P0,Kff=Kff_P0,warn="n")$kSWS2total
  list(K1 =as.numeric(K1 (S=S,T=T,P=0,pHscale="T",k1k2="l",s2c,t2s,warn="n")),
       K2 =as.numeric(K2 (S=S,T=T,P=0,pHscale="T",k1k2="l",s2c,t2s,warn="n")),
       Kb =as.numeric(Kb (S=S,T=T,P=0,pHscale="T",s2c,t2s,warn="n")),
       Kw =as.numeric(Kw (S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
       Ksi=as.numeric(Ksi(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
       K1p=as.numeric(K1p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
       K2p=as.numeric(K2p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
       K3p=as.numeric(K3p(S=S,T=T,P=0,pHscale="T",s2c,warn="n")),
       Ks =as.numeric(Ks_P0),
       Kf =as.numeric(Kff_P0),
       BOR=as.numeric(bor(S=S,b="u74")))
}

## dK/dT, dK/dS reassembled the way buffderiv() reassembles them
dK <- function(S, T) {
  K <- scK(S,T)
  L <- dlnK(S, T, K$Ks, K$Kf)
  list(K1_T =K$K1 *L$K1_T,             K1_S =K$K1 *L$K1_S,
       K2_T =K$K2 *L$K2_T,             K2_S =K$K2 *L$K2_S,
       Kb_T =K$Kb *L$Kb_T,             Kb_S =K$Kb *L$Kb_S,
       Ks_T =K$Ks *L$Ks_T,             Ks_S =K$Ks *L$Ks_S,
       Kf_T =K$Kf *L$Kf_T,             Kf_S =K$Kf *L$Kf_S,
       Kw_T =K$Kw *(L$Kw_T +L$Sig_T),  Kw_S =K$Kw *(L$Kw_S +L$Sig_S),
       Ksi_T=K$Ksi*(L$Ksi_T+L$Sig_T),  Ksi_S=K$Ksi*(L$Ksi_S+L$Sig_S),
       K1p_T=K$K1p*(L$K1p_T+L$Sig_T),  K1p_S=K$K1p*(L$K1p_S+L$Sig_S),
       K2p_T=K$K2p*(L$K2p_T+L$Sig_T),  K2p_S=K$K2p*(L$K2p_S+L$Sig_S),
       K3p_T=K$K3p*(L$K3p_T+L$Sig_T),  K3p_S=K$K3p*(L$K3p_S+L$Sig_S),
       BOR_T=0,                        BOR_S=K$BOR/S)
}

nm <- c("K1","K2","Kb","Ks","Kf","Kw","Ksi","K1p","K2p","K3p","BOR")
TS <- list(c(-1,34), c(5,33.5), c(15,30), c(25,35), c(31,37))

cat("CHECK 1: dK/dT and dK/dS from dlnK.R  vs  Richardson FD of seacarb's own K()\n\n")
cat(sprintf("%6s%6s%5s","T","S","wrt")); for (k in nm) cat(sprintf("%10s",k)); cat("\n")
w1 <- 0
for (ts in TS) {
  d <- dK(ts[2], ts[1])
  for (v in c("T","S")) {
    cat(sprintf("%6.1f%6.1f%5s", ts[1], ts[2], v))
    for (k in nm) {
      num <- if (v=="T") rich(function(x) scK(ts[2],x)[[k]], ts[1], 1e-3)
             else        rich(function(x) scK(x,ts[1])[[k]], ts[2], 1e-3)
      ana <- d[[paste0(k,"_",v)]]
      rd  <- if (num==0 && ana==0) 0 else 100*(ana-num)/num
      w1  <- max(w1, abs(rd)); cat(sprintf("%10.1e", rd))
    }
    cat("\n")
  }
}
cat(sprintf("\nWORST: %.3e %%   (finite-difference noise floor)\n\n", w1))

cases <- data.frame(
  T   = c(-1.0,  2.0, 10.0, 20.0, 25.0, 29.0, 28.0, 15.0, 25.0,  5.0),
  S   = c(34.0, 34.5, 34.8, 35.0, 35.0, 36.5, 33.0, 30.0, 35.0, 33.5),
  CT  = c(2150, 2100, 2050, 2000, 2000, 1950, 1900, 1950, 2100, 2200)*1e-6,
  AT  = c(2300, 2280, 2290, 2300, 2300, 2380, 2150, 2100, 2300, 2320)*1e-6,
  Pt  = c(1.5, 1.2, 0.5, 0.2, 0.1, 0.05, 0.3, 0.8, 0.1, 1.8)*1e-6,
  Sit = c( 60,  40,  10,   3,   2,   1,   5,  15,   2,  70)*1e-6,
  lab = c("polar","subpolar","temperate","notebook_base","subtropical",
          "warm_salty_gyre","warm_fresh","brackish_coastal",
          "high_CO2_2100","southern_ocean"), stringsAsFactors=FALSE)

bd <- function(S,T,CT,AT,Pt,Sit)
  buffderiv(flag=15, var1=AT, var2=CT, S=S, T=T, Pt=Pt, Sit=Sit,
            k1k2="l", kf="dg", ks="d", pHscale="T", b="u74", warn="n")

cat("CHECK 2: alkalinity residual, and Be vs buffsun() at Pt=Sit=0\n\n")
w2 <- 0; wBe <- 0
for (i in 1:nrow(cases)) {
  a  <- bd(cases$S[i],cases$T[i],cases$CT[i],cases$AT[i],cases$Pt[i],cases$Sit[i])
  a0 <- bd(cases$S[i],cases$T[i],cases$CT[i],cases$AT[i],0,0)
  b0 <- buffsun(flag=15, var1=cases$AT[i], var2=cases$CT[i], S=cases$S[i],
                T=cases$T[i], P=0, Pt=0, Sit=0, k1k2="l", kf="dg", ks="d",
                pHscale="T", b="u74", warn="n")
  rr <- abs(a$alk_residual)/cases$AT[i]
  w2 <- max(w2, rr); wBe <- max(wBe, abs(100*(a0$Be-b0$Be)/b0$Be))
  cat(sprintf("  %-17s |resid|/AT=%9.2e   Be=%10.6f\n", cases$lab[i], rr, a$Be))
}
cat(sprintf("\nWORST |alk residual|/AT : %.3e  -> %s\n", w2,
    ifelse(w2 < 1e-12, "Anc IS carb()'s alkalinity", "*** MISMATCH ***")))
cat(sprintf("WORST Be vs buffsun(), Pt=Sit=0 : %.3e %%  (buffsun inherits carb()'s h tolerance)\n\n", wBe))

cat("CHECK 3: analytic dBe/dY  vs  Richardson FD through carb()\n")
cat("         Pt and Sit nonzero: exercises phosphate + monoprotic silicate\n")
w3 <- 0
for (i in 1:nrow(cases)) {
  T0<-cases$T[i]; S0<-cases$S[i]; CT0<-cases$CT[i]; AT0<-cases$AT[i]
  P0<-cases$Pt[i]; Si0<-cases$Sit[i]
  a  <- bd(S0,T0,CT0,AT0,P0,Si0)
  Bf <- function(S,T,CT,AT) bd(S,T,CT,AT,P0,Si0)$Be
  num <- c(dBe_dT  = rich(function(x) Bf(S0,x,CT0,AT0), T0, 1e-3),
           dBe_dS  = rich(function(x) Bf(x,T0,CT0,AT0), S0, 1e-3),
           dBe_dCT = rich(function(x) Bf(S0,T0,x,AT0), CT0, CT0*1e-5),
           dBe_dAT = rich(function(x) Bf(S0,T0,CT0,x), AT0, AT0*1e-5))
  cat(sprintf("\n--- %s: T=%.1f S=%.1f CT=%.0f AT=%.0f Pt=%.2f Sit=%.0f  Be=%.6f\n",
      cases$lab[i], T0,S0,CT0*1e6,AT0*1e6,P0*1e6,Si0*1e6, a$Be))
  for (k in c("dBe_dT","dBe_dS","dBe_dCT","dBe_dAT")) {
    an <- a[[k]]; nu <- num[[k]]; rd <- 100*(an-nu)/nu
    w3 <- max(w3, abs(rd))
    cat(sprintf("    %-8s analytic=%19.10e  FD=%19.10e  rel.diff=%10.2e%%\n", k, an, nu, rd))
  }
}
cat(sprintf("\nWORST |rel.diff| dBe/dY vs FD : %.3e %%   (bar 1e-4 %%)  -> %s\n",
    w3, ifelse(w3 < 1e-4, "PASS", "FAIL")))
