# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# buffsun: computes the Bolin sensitivity Be = d(DIC)/d[CO2*] and the classical Revelle
#          factor Rf, at constant total alkalinity (AT), temperature (T), and salinity (S),
#          using the condensed closed-form template shown to be the common algebraic result
#          of four independent derivations: Sundquist, Plummer & Wigley (1979); Zeebe &
#          Wolf-Gladrow (2001, corrected); Egleston et al. (2010, corrected); and
#          Frankignoulle (1994, as correctly implemented in seacarb's buffer()).
#
# --------------------------------------------------------------------------------------------
# BACKGROUND
# --------------------------------------------------------------------------------------------
#
# Writing s = [CO2*], b = [HCO3-], c = [CO3--], and X for the non-carbonate buffering term,
# all four of the above reduce algebraically to:
#
#     Be = 1 + c/s + (X*b - 4*c^2) / (s*(b + 4*c + X))
#     Rf = DIC / (s * Be)
#
# with X = h*w, where w is built additively, one acid-base system at a time (identical
# construction to seacarb's buffzwg() and buffesm()):
#
#     w = w_bw + w_Si + w_P + w_N + w_S
#
#     w_bw = 1 + Kw/h^2 + BT*Kb/(Kb+h)^2                  (borate + water)
#     w_Si = Sit*Ksi*(h^2+4*K2si*h+Ksi*K2si) / (h^2+Ksi*h+Ksi*K2si)^2
#                                                          (silicic acid, full diprotic system,
#                                                           matching calculate_carb.R's own
#                                                           Asi = [SiO(OH)3-] + 2[SiO2(OH)2--])
#     w_P  = phosphoric acid (three dissociation steps; see buffzwg.R for the NP/DP form)
#     w_N  = NH4t*Kn/(Kn+h)^2                              (ammonia)
#     w_S  = HSt*Khs/(Khs+h)^2                             (sulfide)
#
# Unlike buffzwg() and buffesm(), this routine needs no explicit K1 or K2: since b, c, and s
# (i.e. HCO3, CO3, and CO2) are taken directly from carb()/carbfull()'s own speciation output,
# the equation is written purely in terms of the three dissolved inorganic carbon species and
# X. This is the point of the condensed form: it makes explicit that Be depends on the
# carbonate system only through s, b, and c, with all non-carbonate chemistry isolated in X.
#
# Setting X = 0 collapses this exactly to the constant-Ac (not constant-AT) special case,
# Be = 1 + b*c/(s*(b+4*c)) -- the form obtained if carbonate alkalinity, rather than total
# alkalinity, is (incorrectly) assumed to remain fixed as CO2 is added to the ocean.
#
# --------------------------------------------------------------------------------------------

buffsun <-
  function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, NH4t=0, HSt=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", warn="y", eos="eos80", long=1.e20, lat=1.e20){

    n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(NH4t), length(HSt), length(k1k2), length(kf), length(pHscale), length(ks), length(b))
    if(length(flag)!=n){flag <- rep(flag[1],n)}
    if(length(var1)!=n){var1 <- rep(var1[1],n)}
    if(length(var2)!=n){var2 <- rep(var2[1],n)}
    if(length(S)!=n){S <- rep(S[1],n)}
    if(length(T)!=n){T <- rep(T[1],n)}
    if(length(Patm)!=n){Patm <- rep(Patm[1],n)}
    if(length(P)!=n){P <- rep(P[1],n)}
    if(length(Pt)!=n){Pt <- rep(Pt[1],n)}
    if(length(Sit)!=n){Sit <- rep(Sit[1],n)}
    if(length(NH4t)!=n){NH4t <- rep(NH4t[1],n)}
    if(length(HSt)!=n){HSt <- rep(HSt[1],n)}
    if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
    if(length(kf)!=n){kf <- rep(kf[1],n)}
    if(length(ks)!=n){ks <- rep(ks[1],n)}
    if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
    if(length(b)!=n){b <- rep(b[1],n)}

    # if the concentrations of total silicate, total phosphate, total ammonia, and
    # total sulfide are NA, they are set at 0
    Sit[is.na(Sit)]  <- 0
    Pt[is.na(Pt)]    <- 0
    NH4t[is.na(NH4t)] <- 0
    HSt[is.na(HSt)]   <- 0

    # Only two options for eos
    if (eos != "teos10" && eos != "eos80")
        stop ("invalid parameter eos: ", eos)

    # if use of EOS-10 standard
    if (eos == "teos10")
    {
        # Must convert temperature and salinity from TEOS-10 to EOS-80
        # convert temperature: from Conservative (CT) to in-situ temperature
        # and salinity from Absolute to Practical (SP)
        eos <- teos2eos_geo (S, T, P, long, lat)
        InsT <- eos$T
        SP <- eos$SP
    }
    else
    {
        InsT <- T
        SP <- S
    }

    # carb() hard-codes NH4t=0, HSt=0 internally and does not expose them as arguments;
    # only carbfull() does. So if the caller has supplied any nonzero NH4t or HSt, we must
    # use carbfull() to get pH/CO2/HCO3/CO3/DIC that are actually consistent with those
    # inputs. Otherwise carb() is used, for speed and backward compatibility (verified to
    # agree with carbfull(NH4t=0, HSt=0) to solver tolerance, ~1e-9 relative).
    use_carbfull <- any(NH4t != 0) || any(HSt != 0)
    if (use_carbfull) {
        Carb <- carbfull(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, NH4t=NH4t, HSt=HSt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    } else {
        Carb <- carb(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    }
    P    <- Carb$P
    pH   <- Carb$pH
    h    <- 10^(-pH)
    s    <- Carb$CO2    # s := [CO2*]
    b_   <- Carb$HCO3   # b := [HCO3-]  (renamed locally to avoid clashing with the 'b' boron-formulation argument)
    c    <- Carb$CO3    # c := [CO3--]
    DIC  <- Carb$DIC

    #-------Constants----------------
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK = InsT + tk;           # TK [K]; InsT[C]

    Cl = SP / 1.80655;            # Cl = chlorinity; SP = practical salinity (psu)
    ST = 0.14 * Cl/96.062        # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    FT = 6.7e-5 * Cl/18.9984     # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
    FLUO = 6.7e-5 * Cl/18.9984   # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
    bor = bor(S=SP , b=b)         # (mol/kg) total boron

    #---------------------------------------------------------------------
    #--------------------- compute K's ----------------------------------
    #---------------------------------------------------------------------

    # Ks (free pH scale) at zero pressure and given pressure
    Ks_P0 <- Ks(S=SP, T=InsT, P=0, ks=ks, warn=warn)
    Ks    <- Ks(S=SP, T=InsT, P=P, ks=ks, warn=warn)

    # Kf on free pH scale
    Kff_P0 <- Kf(S=SP, T=InsT, P=0, pHscale="F", kf=kf, Ks_P0, Ks)
    Kff <- Kf(S=SP, T=InsT, P=P, pHscale="F", kf=kf, Ks_P0, Ks)
    # Kf on given pH scale
    Kf <- Kf(S=SP, T=InsT, P=P, pHscale=pHscale, kf=kf, Ks_P0, Ks)

    # Conversion factor from total to SWS pH scale at zero pressure
    ktotal2SWS_P0 <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks_P0,Kff=Kff_P0)$ktotal2SWS

    # Conversion factor from SWS to chosen pH scale
    conv <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks,Kff=Kff)
    kSWS2chosen <- rep(1.,n)
    kSWS2chosen [pHscale == "T"] <- conv$kSWS2total [pHscale == "T"]
    kSWS2chosen [pHscale == "F"] <- conv$kSWS2free [pHscale == "F"]

    # Note: unlike buffzwg() and buffesm(), K1 and K2 are NOT needed here -- the condensed
    # equation is written directly in terms of s, b (=HCO3), and c (=CO3), which already
    # encode the carbonate equilibrium via carb()/carbfull()'s own speciation.
    Kw  <- Kw(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Kb  <- Kb(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    K1p <- K1p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K2p <- K2p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K3p <- K3p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Ksi <- Ksi(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K2si <- K2si(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0)
    # --- silicate: monoprotic or diprotic MUST follow the solver -----------------
    # calculate_carb.R branches on fullresult:
    #   fullresult=FALSE  (carb(),     used when NH4t=HSt=0) -> siooh3 only, MONOPROTIC
    #   fullresult=TRUE   (carbfull(), used otherwise)        -> siooh3 + 2*sio2oh2, DIPROTIC
    # The diprotic wSi below reduces EXACTLY (bit for bit) to the monoprotic form
    # Sit*Ksi/(Ksi+h)^2 in the limit K2si -> 0, so zeroing K2si on the carb() path is
    # all that is needed. Without this, wSi is differentiating an alkalinity that
    # carb() does not use (error ~6e-8 of AT at Sit = 60 umol/kg).
    if (!use_carbfull) K2si <- 0*K2si
    Kn  <- Kn(S=SP, T=InsT, P=P, pHscale=pHscale, warn=warn)
    Khs <- Khs(S=SP, T=InsT, P=P, pHscale=pHscale, warn=warn)

    #---------------------------------------------------------------------
    #--------------------- Bolin sensitivity Be and Revelle factor Rf -----
    #---------------------------------------------------------------------

    # Non-carbonate buffering term w, built additively, one acid-base system at a time
    # (identical construction to buffzwg.R and buffesm.R).
    wbw <- 1 + Kw/h^2 + bor*Kb/(Kb+h)^2

    DSi <- h^2 + Ksi*h + Ksi*K2si
    wSi <- Sit*Ksi*(h^2 + 4*K2si*h + Ksi*K2si) / DSi^2

    NP <- K1p*K2p*h + 2*K1p*K2p*K3p - h^3
    DP <- h^3 + K1p*h^2 + K1p*K2p*h + K1p*K2p*K3p
    wP <- Pt * ( (3*h^2 + 2*K1p*h + K1p*K2p)*NP - (K1p*K2p - 3*h^2)*DP ) / DP^2

    wN <- NH4t*Kn/(Kn+h)^2
    wS <- HSt*Khs/(Khs+h)^2


    # --- fluoride: -[HF] is a species in Dickson's total alkalinity ---------------
    # SolveSAPHE (which carb()/carbfull() invert) carries it explicitly:
    #   aphscale = 1 + ST/Ks ;  api1_flu = Kf_free*aphscale ;  A_HF = -FT*h/(api1_flu+h)
    # so w = -dAlk/dh must include +FT*Kfa/(Kfa+h)^2.  Kf MUST be put on the same pH
    # scale as h (a 28% factor); using the free-scale Kf with a total-scale h is wrong.
    aph <- 1 + ST/Ks                 # SolveSAPHE aphscale (total pH scale)
    Kfa <- Kff*aph                   # Kf on the chosen (total) pH scale
    wF  <- FT*Kfa/(Kfa + h)^2

    w <- wbw + wSi + wP + wN + wS + wF

    # X is the non-carbonate term in the condensed template (X = h*w).
    X <- h * w

    # The condensed common-template equation (see BACKGROUND above), written directly in
    # terms of the three dissolved inorganic carbon species and X:
    Be <- 1 + c/s + (X*b_ - 4*c^2) / (s*(b_ + 4*c + X))
    Rf <- DIC / (s * Be)

    col <- c("Be", "Rf")
    res <- data.frame(Be, Rf)
    names(res) <- col

    return(res)
  }
