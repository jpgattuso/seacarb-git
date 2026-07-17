# This file is part of seacarb.
#
# buffderiv: analytic partial derivatives of the Bolin sensitivity
#            Be = dCT/d[CO2*]  at constant AT, T, S
#
#   dBe/dT   [per degC]      dBe/dCT  [per (mol/kg)]
#   dBe/dS   [per psu]       dBe/dAT  [per (mol/kg)]
#
# Surface only (P = 0).  Acid-base systems: carbonate, borate, water, phosphate
# and silicate.  Silicate is MONOPROTIC (K1si only), exactly as in carb(), which
# calls calculate_carb() with fullresult=FALSE and therefore never forms K2si.
#
# ------------------------------------------------------------------------------
# METHOD
# ------------------------------------------------------------------------------
# h = [H+] is fixed implicitly by the alkalinity constraint
#
#   F(h; T,S,CT,AT) = Ac(h,CT) + Anc(h) - AT = 0
#   Ac  = CT (K1 h + 2 K1 K2)/Dc,          Dc = h^2 + K1 h + K1 K2
#   Anc = BOR Kb/(Kb+h) + Kw/h - h + Sit Ksi/(Ksi+h) + Pt NP/DP
#
# and the buffer term of buffsun() is exactly minus the slope of Anc:
#
#   w = -dAnc/dh = 1 + Kw/h^2 + BOR Kb/(Kb+h)^2 + Sit Ksi/(Ksi+h)^2 + wP
#   X = h w
#
# Be is then the EXPLICIT algebraic function G(h,CT,T,S) evaluated by buffsun().
# For Y in {T,S,CT,AT} the implicit function theorem gives
#
#   dBe/dY = G_Y + G_h ( -F_Y / F_h )
#
# Everything is closed form.  The only nonlinear solve is the single carb() call
# that supplies h.
#
# ------------------------------------------------------------------------------
# CONSTANTS
# ------------------------------------------------------------------------------
# All K VALUES come from seacarb's own K1(), K2(), Kb(), Kw(), Ksi(), K1p(),
# K2p(), K3p(), Ks(), Kf(), kconv() and bor(), fetched exactly as
# calculate_carb() fetches them.  Nothing is transcribed or duplicated here.
#
# What seacarb does NOT supply is dK/dT and dK/dS.  Those are in dlnK.R, as
# logarithmic derivatives, and are reassembled below via
#
#   total-scale-native (K1, K2, Kb):   dK/dY = K (dlnK/dY)
#   SWS-native (Kw, Ksi, K1p..K3p):    dK/dY = K (dlnKsws/dY + dlnSig/dY)
#
# Because dlnK.R hard-codes particular formulations, buffderiv() STOPS unless
# the matching options are passed.  Two of those guards are not pedantry:
#
#   k1k2="x"  -> K1.R does:  is_outrange <- T>35 | T<2 | S<19 | S>43
#                            k1k2[is_x & is_outrange] <- "w14"
#                i.e. it silently switches to Waters et al. (2014) below 2 degC,
#                which is exactly where the polar CMIP6 surface cells live.
#
#   kf="x"    -> silently switches between Perez & Fraga and Dickson & Riley by
#                T,S range.  That changes kSWS2total and hence Kw, Ksi, K1p, K2p
#                and K3p by 0.5 to 0.9 %.
#
# Neither switch emits a warning.  Pass k1k2="l" and kf="dg" explicitly.
# ------------------------------------------------------------------------------

buffderiv <- function(flag, var1, var2, S=35, T=25, Patm=1, Pt=0, Sit=0,
                      NH4t=0, HSt=0,
                      k1k2="l", kf="dg", ks="d", pHscale="T", b="u74", warn="y",
                      npolish=3)
{
    if (any(k1k2 != "l"))
        stop('buffderiv: dlnK.R is hard-wired to Lueker et al. (2000). Pass k1k2="l", NOT "x".')
    if (any(kf != "dg"))
        stop('buffderiv: dlnK.R is hard-wired to Dickson & Riley (1979). Pass kf="dg", NOT "x".')
    if (any(ks != "d"))
        stop('buffderiv: dlnK.R is hard-wired to Dickson (1990b). Pass ks="d".')
    if (any(b != "u74"))
        stop('buffderiv: dlnK.R is hard-wired to Uppstrom (1974). Pass b="u74".')
    if (any(pHscale != "T"))
        stop('buffderiv: total pH scale only. Pass pHscale="T".')

    n <- max(length(flag), length(var1), length(var2), length(S), length(T),
             length(Patm), length(Pt), length(Sit), length(NH4t), length(HSt))
    if (length(flag)!=n) flag <- rep(flag[1], n)
    if (length(var1)!=n) var1 <- rep(var1[1], n)
    if (length(var2)!=n) var2 <- rep(var2[1], n)
    if (length(S)   !=n) S    <- rep(S[1],    n)
    if (length(T)   !=n) T    <- rep(T[1],    n)
    if (length(Patm)!=n) Patm <- rep(Patm[1], n)
    if (length(Pt)  !=n) Pt   <- rep(Pt[1],   n)
    if (length(Sit) !=n) Sit  <- rep(Sit[1],  n)
    if (length(NH4t)!=n) NH4t <- rep(NH4t[1], n)
    if (length(HSt) !=n) HSt  <- rep(HSt[1],  n)

    #---------------------------------------------------------------------
    #-------------- speciation from seacarb ------------------------------
    #---------------------------------------------------------------------
    # carb() hard-codes NH4t=HSt=0 AND treats silicate as monoprotic
    # (calculate_carb, fullresult=FALSE). Only carbfull() carries ammonia, sulfide
    # and the second silicate dissociation. Pick the solver, and let K2si follow it.
    use_carbfull <- any(NH4t != 0) || any(HSt != 0)
    if (use_carbfull) {
        Carb <- carbfull(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=0,
                         Pt=Pt, Sit=Sit, NH4t=NH4t, HSt=HSt, k1k2=k1k2, kf=kf,
                         ks=ks, pHscale=pHscale, b=b)
    } else {
        Carb <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=0,
                     Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale,
                     b=b, warn=warn)
    }
    h  <- 10^(-Carb$pH)
    CT <- Carb$DIC
    ALK <- Carb$ALK
    # CO2, HCO3 and CO3 are NOT taken from carb(); they are recomputed below from
    # the polished h, so that speciation, Be and the derivatives all follow from
    # the SAME h to machine precision.

    #---------------------------------------------------------------------
    #-------------- compute K's and Boron, as calculate_carb() does -------
    #---------------------------------------------------------------------
    BOR <- as.numeric(bor(S=S, b=b))

    # Ks (free pH scale).  P = 0, so Ks_P0 and Ks coincide.
    Ks_P0 <- Ks(S=S, T=T, P=0, ks=ks, warn=warn)

    # Kf on free pH scale
    Kff_P0 <- Kf(S=S, T=T, P=0, pHscale="F", kf=kf, Ks_P0, Ks_P0, warn=warn)

    # Conversion factor from total to SWS pH scale at zero pressure
    ktotal2SWS_P0 <- kconv(S=S, T=T, P=0, kf=kf, Ks=Ks_P0, Kff=Kff_P0)$ktotal2SWS

    # Conversion factor from SWS to the chosen (total) pH scale
    conv        <- kconv(S=S, T=T, P=0, kf=kf, Ks=Ks_P0, Kff=Kff_P0, warn=warn)
    kSWS2chosen <- conv$kSWS2total

    K1  <- as.numeric(K1 (S=S, T=T, P=0, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn))
    K2  <- as.numeric(K2 (S=S, T=T, P=0, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn))
    Kw  <- as.numeric(Kw (S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, warn=warn))
    Kb  <- as.numeric(Kb (S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn))
    K1p <- as.numeric(K1p(S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, warn=warn))
    K2p <- as.numeric(K2p(S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, warn=warn))
    K3p <- as.numeric(K3p(S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, warn=warn))
    Ksi <- as.numeric(Ksi(S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, warn=warn))
    Ksv <- as.numeric(Ks_P0)
    Kfv <- as.numeric(Kff_P0)

    # Silicate: MONOPROTIC on the carb() path, DIPROTIC on the carbfull() path.
    # The diprotic wSi reduces exactly to Sit*Ksi/(Ksi+h)^2 as K2si -> 0, so zeroing
    # K2si is all that is needed to switch.
    if (use_carbfull) {
        K2si <- as.numeric(K2si(S=S, T=T, P=0, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0))
        Kn   <- as.numeric(Kn (S=S, T=T, P=0, pHscale=pHscale, warn=warn))
        Khs  <- as.numeric(Khs(S=S, T=T, P=0, pHscale=pHscale, warn=warn))
    } else {
        K2si <- rep(0, n) ; Kn <- rep(1, n) ; Khs <- rep(1, n)   # inert; NH4t=HSt=0
    }

    # Sulfate and fluoride. These are NOT bookkeeping for the pH scale: they are
    # species in Dickson's total alkalinity, and SolveSAPHE (which carb() inverts)
    # carries them explicitly. From solve_pH_from_AT.R, for pHscale="T":
    #     aphscale = 1 + ST/Ks ;  api1_so4 = Ks*aphscale ;  api1_flu = Kf*aphscale
    # and from equation_at.R / anw.R the alkalinity contains
    #     Kw/h - h/aphscale - ST*h/(api1_so4 + h) - FT*h/(api1_flu + h)
    # The first two collapse to -h to within ~1e-16, but -[HF] does NOT: it is
    # worth up to 5e-7 of AT, and omitting it displaces h by ~6e-7 relative.
    Cl  <- S/1.80655
    ST  <- 0.14*Cl/96.062         # total sulfate  (Dickson et al. 2007, Table 2)
    FT  <- 6.7e-5*Cl/18.9984      # total fluoride (Dickson et al. 2007, Table 2)
    aph <- 1 + ST/Ksv             # SolveSAPHE aphscale, total pH scale
    Ksa <- Ksv*aph                # == Ks + ST
    Kfa <- Kfv*aph

    #---------------------------------------------------------------------
    #-------------- dK/dT and dK/dS --------------------------------------
    #---------------------------------------------------------------------
    L <- dlnK(S, T, Ksv, Kfv)

    dK1_dT  <- K1 *L$K1_T ;  dK1_dS  <- K1 *L$K1_S            # total-scale native
    dK2_dT  <- K2 *L$K2_T ;  dK2_dS  <- K2 *L$K2_S
    dKb_dT  <- Kb *L$Kb_T ;  dKb_dS  <- Kb *L$Kb_S

    dKw_dT  <- Kw *(L$Kw_T  + L$Sig_T) ; dKw_dS  <- Kw *(L$Kw_S  + L$Sig_S)   # SWS native
    dKsi_dT <- Ksi*(L$Ksi_T + L$Sig_T) ; dKsi_dS <- Ksi*(L$Ksi_S + L$Sig_S)
    dK1p_dT <- K1p*(L$K1p_T + L$Sig_T) ; dK1p_dS <- K1p*(L$K1p_S + L$Sig_S)
    dK2p_dT <- K2p*(L$K2p_T + L$Sig_T) ; dK2p_dS <- K2p*(L$K2p_S + L$Sig_S)
    dK3p_dT <- K3p*(L$K3p_T + L$Sig_T) ; dK3p_dS <- K3p*(L$K3p_S + L$Sig_S)

    dBOR_dT <- 0                                              # bor() is linear in S
    dBOR_dS <- BOR/S

    dK2si_dT <- K2si*L$K2si_T ; dK2si_dS <- K2si*L$K2si_S     # total-scale native
    dKhs_dT  <- Khs *L$Khs_T  ; dKhs_dS  <- Khs *L$Khs_S      # total-scale native
    # Kn is SWS native, but seacarb converts it with kf="x" (see dlnK.R), NOT with our
    # kf="dg". We must use the SAME factor seacarb used, or dKn/dY is wrong by up to 8%.
    dKn_dT   <- Kn*(L$Kn_T + L$Sigx_T) ; dKn_dS <- Kn*(L$Kn_S + L$Sigx_S)

    dKs_dT <- Ksv*L$Ks_T ; dKs_dS <- Ksv*L$Ks_S               # free scale
    dKf_dT <- Kfv*L$Kf_T ; dKf_dS <- Kfv*L$Kf_S               # free scale
    dST_dT <- 0 ; dST_dS <- 0.14/96.062/1.80655
    dFT_dT <- 0 ; dFT_dS <- 6.7e-5/18.9984/1.80655

    #---------------------------------------------------------------------
    #-------------- Newton polish of h -----------------------------------
    #---------------------------------------------------------------------
    # carb() (SolveSAPHE) converges h to about 1e-8 relative.  That is ample for h
    # itself, but Be is steep in h and its derivatives steeper still: an FD check
    # of dBe/dT with eps = 1e-3 has amplification Be/(|dBe/dT| 2 eps) ~ 1e5, so a
    # 1e-8 error in h masquerades as a 1e-3 % error in the derivative.  Two or
    # three Newton steps on F(h) = 0, using F_h which we need anyway, cost nothing
    # and take h to machine precision.
    for (it in seq_len(npolish)) {
        Dc   <- h^2 + K1*h + K1*K2
        Ac   <- CT*(K1*h + 2*K1*K2)/Dc
        DP   <- h^3 + K1p*h^2 + K1p*K2p*h + K1p*K2p*K3p
        NP   <- K1p*K2p*h + 2*K1p*K2p*K3p - h^3
        DP_h <- 3*h^2 + 2*K1p*h + K1p*K2p
        NP_h <- K1p*K2p - 3*h^2
        DSi  <- h^2 + Ksi*h + Ksi*K2si
        Anc  <- BOR*Kb/(Kb+h) + Kw/h + Sit*(Ksi*h + 2*Ksi*K2si)/DSi + Pt*NP/DP +
                NH4t*Kn/(Kn+h) + HSt*Khs/(Khs+h) -
                h/aph - ST*h/(Ksa+h) - FT*h/(Kfa+h)
        w    <- 1/aph + ST*Ksa/(Ksa+h)^2 + FT*Kfa/(Kfa+h)^2 +
                Kw/h^2 + BOR*Kb/(Kb+h)^2 +
                Sit*Ksi*(h^2 + 4*K2si*h + Ksi*K2si)/DSi^2 +
                NH4t*Kn/(Kn+h)^2 + HSt*Khs/(Khs+h)^2 +
                Pt*(NP*DP_h - NP_h*DP)/DP^2
        F_h  <- -CT*K1*(h^2 + 4*K2*h + K1*K2)/Dc^2 - w
        h    <- h - (Ac + Anc - ALK)/F_h
    }

    #---------------------------------------------------------------------
    #-------------- speciation, w and X at the polished h ----------------
    #---------------------------------------------------------------------
    Dc   <- h^2 + K1*h + K1*K2
    CO2  <- CT*h^2/Dc
    HCO3 <- CT*K1*h/Dc
    CO3  <- CT*K1*K2/Dc

    # phosphate: Ap = Pt NP/DP ;  wP = -dAp/dh
    DP    <- h^3 + K1p*h^2 + K1p*K2p*h + K1p*K2p*K3p
    NP    <- K1p*K2p*h + 2*K1p*K2p*K3p - h^3
    DP_h  <- 3*h^2 + 2*K1p*h + K1p*K2p
    NP_h  <- K1p*K2p - 3*h^2
    DP_hh <- 6*h + 2*K1p
    NP_hh <- -6*h
    Ap    <- Pt*NP/DP
    wP    <- Pt*(NP*DP_h - NP_h*DP)/DP^2
    wP_h  <- Pt*((NP*DP_hh - NP_hh*DP)*DP - 2*(NP*DP_h - NP_h*DP)*DP_h)/DP^3

    NP_1 <- K2p*h + 2*K2p*K3p ; DP_1 <- h^2 + K2p*h + K2p*K3p ; NP_h1 <- K2p ; DP_h1 <- 2*h + K2p
    NP_2 <- K1p*h + 2*K1p*K3p ; DP_2 <- K1p*h + K1p*K3p       ; NP_h2 <- K1p ; DP_h2 <- K1p
    NP_3 <- 2*K1p*K2p         ; DP_3 <- K1p*K2p               ; NP_h3 <- 0   ; DP_h3 <- 0

    dAp_dK1p <- Pt*(NP_1*DP - NP*DP_1)/DP^2
    dAp_dK2p <- Pt*(NP_2*DP - NP*DP_2)/DP^2
    dAp_dK3p <- Pt*(NP_3*DP - NP*DP_3)/DP^2
    dwP_dK1p <- Pt*((NP_1*DP_h + NP*DP_h1 - NP_h1*DP - NP_h*DP_1)*DP - 2*(NP*DP_h - NP_h*DP)*DP_1)/DP^3
    dwP_dK2p <- Pt*((NP_2*DP_h + NP*DP_h2 - NP_h2*DP - NP_h*DP_2)*DP - 2*(NP*DP_h - NP_h*DP)*DP_2)/DP^3
    dwP_dK3p <- Pt*((NP_3*DP_h + NP*DP_h3 - NP_h3*DP - NP_h*DP_3)*DP - 2*(NP*DP_h - NP_h*DP)*DP_3)/DP^3

    # proton + sulfate + fluoride block, exactly as SolveSAPHE writes it
    DS <- Ksa + h
    DF <- Ksv*(Kfa + h)                    # == Kf*Ks + Kf*ST + Ks*h
    Ah   <- -h/aph - ST*h/DS - FT*h/(Kfa+h)
    wA   <- 1/aph + ST*Ksa/DS^2 + FT*Kfa/(Kfa+h)^2
    wA_h <- -2*ST*Ksa/DS^3 - 2*FT*Kfa/(Kfa+h)^3

    # silicate, diprotic (K2si = 0 on the carb() path gives the monoprotic form exactly)
    DSi     <- h^2 + Ksi*h + Ksi*K2si
    ASi     <- Sit*(Ksi*h + 2*Ksi*K2si)/DSi
    wSi     <- Sit*Ksi*(h^2 + 4*K2si*h + Ksi*K2si)/DSi^2
    wSi_h   <- 2*Ksi*Sit*(-(2*K2si + h)*(Ksi + 2*h)^2 + (2*K2si + Ksi + 3*h)*DSi)/DSi^3
    dASi_dKsi  <- Sit*(2*K2si + h)*h^2/DSi^2
    dASi_dK2si <- Ksi*Sit*h*(Ksi + 2*h)/DSi^2
    dwSi_dKsi  <- Sit*(-2*Ksi*(K2si + h)*(2*K2si + h)*(Ksi + 2*h) - DSi^2 +
                   DSi*(Ksi*(K2si + h) + Ksi*(2*K2si + h) + (2*K2si + h)*(Ksi + 2*h)))/DSi^3
    dwSi_dK2si <- Ksi*Sit*(-2*Ksi*(2*K2si + h)*(Ksi + 2*h) + (3*Ksi + 4*h)*DSi)/DSi^3

    # ammonia and sulfide: same functional form as borate
    AN <- NH4t*Kn/(Kn+h)  ; wN <- NH4t*Kn/(Kn+h)^2  ; wN_h <- -2*NH4t*Kn/(Kn+h)^3
    AS <- HSt*Khs/(Khs+h) ; wS <- HSt*Khs/(Khs+h)^2 ; wS_h <- -2*HSt*Khs/(Khs+h)^3
    dAN_dKn  <- NH4t*h/(Kn+h)^2  ; dwN_dKn  <- NH4t*(h - Kn)/(Kn+h)^3
    dAS_dKhs <- HSt *h/(Khs+h)^2 ; dwS_dKhs <- HSt *(h - Khs)/(Khs+h)^3

    Anc <- BOR*Kb/(Kb+h) + Kw/h + ASi + Ap + AN + AS + Ah
    w   <- wA + Kw/h^2 + BOR*Kb/(Kb+h)^2 + wSi + wP + wN + wS
    w_h <- wA_h - 2*Kw/h^3 - 2*BOR*Kb/(Kb+h)^3 + wSi_h + wP_h + wN_h + wS_h
    X   <- h*w

    # partials of Ah and wA wrt Ks, Kf, ST, FT (sympy-derived, transcribed verbatim)
    dAh_dKs <- -FT*Kfv*ST*h/DF^2 + ST*h/DS^2 - ST*h/Ksa^2
    dAh_dKf <-  FT*Ksv*h*Ksa/DF^2
    dAh_dST <-  FT*Kfv*Ksv*h/DF^2 - Ksv*h/DS^2 + Ksv*h/Ksa^2 - h^2/DS^2
    dAh_dFT <- -Ksv*h/DF
    dwA_dKs <-  FT*Kfv^2*Ksv*ST/DF^3 + FT*Kfv^2*ST^2/DF^3 - FT*Kfv*Ksv*ST*h/DF^3 -
                Ksv*ST/DS^3 - ST^2/DS^3 + ST*h/DS^3 + ST/Ksa^2
    dwA_dKf <-  FT*Ksv*Ksa*(Ksv*h - Kfv*Ksa)/DF^3
    dwA_dST <- -FT*Kfv^2*Ksv^2/DF^3 - FT*Kfv^2*Ksv*ST/DF^3 + FT*Kfv*Ksv^2*h/DF^3 -
                Ksv*ST/DS^3 - Ksv*h/DS^3 - Ksv/Ksa^2 - ST^2/DS^3 - h^2/DS^3 + 1/DS
    dwA_dFT <-  Kfv*Ksv*Ksa/DF^2

    dAh_dT <- dAh_dKs*dKs_dT + dAh_dKf*dKf_dT + dAh_dST*dST_dT + dAh_dFT*dFT_dT
    dAh_dS <- dAh_dKs*dKs_dS + dAh_dKf*dKf_dS + dAh_dST*dST_dS + dAh_dFT*dFT_dS
    dwA_dT <- dwA_dKs*dKs_dT + dwA_dKf*dKf_dT + dwA_dST*dST_dT + dwA_dFT*dFT_dT
    dwA_dS <- dwA_dKs*dKs_dS + dwA_dKf*dKf_dS + dwA_dST*dST_dS + dwA_dFT*dFT_dS

    # dX/dKj = h dw/dKj ;  dAnc/dKj
    dX_dKw   <- 1/h                                ; dAnc_dKw  <- 1/h
    dX_dKb   <- h*BOR*(h - Kb)/(Kb+h)^3            ; dAnc_dKb  <- BOR*h/(Kb+h)^2
    dX_dBOR  <- h*Kb/(Kb+h)^2                      ; dAnc_dBOR <- Kb/(Kb+h)
    dX_dKsi  <- h*dwSi_dKsi                        ; dAnc_dKsi <- dASi_dKsi
    dX_dK2si <- h*dwSi_dK2si                       ; dAnc_dK2si<- dASi_dK2si
    dX_dKn   <- h*dwN_dKn                          ; dAnc_dKn  <- dAN_dKn
    dX_dKhs  <- h*dwS_dKhs                         ; dAnc_dKhs <- dAS_dKhs
    dX_dK1p  <- h*dwP_dK1p
    dX_dK2p  <- h*dwP_dK2p
    dX_dK3p  <- h*dwP_dK3p

    #---------------------------------------------------------------------
    #-------------- Be and the partials of G ------------------------------
    #---------------------------------------------------------------------
    r1    <- K1/h
    r2    <- K1*K2/h^2
    alpha <- r1 +   r2
    beta  <- r1 + 2*r2
    gam   <- r1 + 4*r2
    dlt   <- r1 + 8*r2
    Q     <- CO2*gam + X                       # == HCO3 + 4 CO3 + X

    Be <- 1 + CO3/CO2 + (X*HCO3 - 4*CO3^2)/(CO2*(HCO3 + 4*CO3 + X))

    G_CO2  <- -beta^2*X/Q^2
    G_beta <- -2*CO2*beta/Q
    G_gam  <-  CO2^2*beta^2/Q^2
    G_X    <-  CO2*beta^2/Q^2
    # G_alpha = 1

    dCO2_dh <- CO2^2*beta/(h*CT)
    dX_dh   <- w + h*w_h
    G_h     <- G_CO2*dCO2_dh - beta/h + G_beta*(-gam/h) + G_gam*(-dlt/h) + G_X*dX_dh
    G_CT    <- G_CO2*CO2/CT

    dCO2_dK1 <- -CO2*(h + K2)/Dc
    dCO2_dK2 <- -CO2*K1/Dc
    GK1 <- G_CO2*dCO2_dK1 + alpha/K1 + G_beta*beta/K1        + G_gam*gam/K1
    GK2 <- G_CO2*dCO2_dK2 + r2/K2    + G_beta*(2*r2/K2)      + G_gam*(4*r2/K2)

    G_T <- GK1*dK1_dT + GK2*dK2_dT + G_X*(dX_dKw *dKw_dT  + dX_dKb *dKb_dT  +
           dX_dBOR*dBOR_dT + dX_dKsi*dKsi_dT + dX_dK1p*dK1p_dT +
           dX_dK2p*dK2p_dT + dX_dK3p*dK3p_dT + h*dwA_dT +
           dX_dK2si*dK2si_dT + dX_dKn*dKn_dT + dX_dKhs*dKhs_dT)
    G_S <- GK1*dK1_dS + GK2*dK2_dS + G_X*(dX_dKw *dKw_dS  + dX_dKb *dKb_dS  +
           dX_dBOR*dBOR_dS + dX_dKsi*dKsi_dS + dX_dK1p*dK1p_dS +
           dX_dK2p*dK2p_dS + dX_dK3p*dK3p_dS + h*dwA_dS +
           dX_dK2si*dK2si_dS + dX_dKn*dKn_dS + dX_dKhs*dKhs_dS)

    #---------------------------------------------------------------------
    #-------------- partials of the alkalinity constraint F ---------------
    #---------------------------------------------------------------------
    Nc   <- K1*h + 2*K1*K2
    Ac   <- CT*Nc/Dc
    F_h  <- -CT*K1*(h^2 + 4*K2*h + K1*K2)/Dc^2 - w

    dAc_dK1 <- CT*((h + 2*K2)*Dc - Nc*(h + K2))/Dc^2
    dAc_dK2 <- CT*K1*(2*Dc - Nc)/Dc^2

    F_T <- dAc_dK1*dK1_dT + dAc_dK2*dK2_dT +
           dAnc_dKw *dKw_dT  + dAnc_dKb *dKb_dT + dAnc_dBOR*dBOR_dT +
           dAnc_dKsi*dKsi_dT +
           dAp_dK1p*dK1p_dT + dAp_dK2p*dK2p_dT + dAp_dK3p*dK3p_dT + dAh_dT +
           dAnc_dK2si*dK2si_dT + dAnc_dKn*dKn_dT + dAnc_dKhs*dKhs_dT
    F_S <- dAc_dK1*dK1_dS + dAc_dK2*dK2_dS +
           dAnc_dKw *dKw_dS  + dAnc_dKb *dKb_dS + dAnc_dBOR*dBOR_dS +
           dAnc_dKsi*dKsi_dS +
           dAp_dK1p*dK1p_dS + dAp_dK2p*dK2p_dS + dAp_dK3p*dK3p_dS + dAh_dS +
           dAnc_dK2si*dK2si_dS + dAnc_dKn*dKn_dS + dAnc_dKhs*dKhs_dS
    F_CT <- Ac/CT

    #---------------------------------------------------------------------
    #-------------- assemble ----------------------------------------------
    #---------------------------------------------------------------------
    dh_dT  <- -F_T /F_h
    dh_dS  <- -F_S /F_h
    dh_dCT <- -F_CT/F_h
    dh_dAT <-  1   /F_h

    dBe_dT  <- G_T  + G_h*dh_dT
    dBe_dS  <- G_S  + G_h*dh_dS
    dBe_dCT <- G_CT + G_h*dh_dCT
    dBe_dAT <-        G_h*dh_dAT

    # Diagnostic. Ac(h) + Anc(h) - ALK must vanish to machine precision. If it did
    # not, our Anc would not be the alkalinity that carb() itself solves, and every
    # derivative above would be the derivative of the wrong function. A structural
    # error, e.g. diprotic instead of monoprotic silicate at Sit = 60 umol/kg,
    # shows up here at the 1e-6 level.
    alk_residual <- Ac + Anc - ALK

    data.frame(Be, dBe_dT, dBe_dS, dBe_dCT, dBe_dAT,
               dh_dT, dh_dS, dh_dCT, dh_dAT, alk_residual)
}
