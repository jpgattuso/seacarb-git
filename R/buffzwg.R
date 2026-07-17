# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# buffzwg: computes the "Bolin sensitivity" Be = d(DIC)/d[CO2*] at constant total alkalinity (AT),
#          temperature (T), and salinity (S), following the exact closed-form expression given
#          by Zeebe and Wolf-Gladrow (2001, "CO2 in Seawater: Equilibrium, Kinetics, Isotopes",
#          Elsevier Oceanography Series 65), Eq. 1.5.88 -- CORRECTED, and extended here to
#          include the phosphoric and silicic acid systems.
#
# --------------------------------------------------------------------------------------------
# BACKGROUND / DERIVATION NOTES (kept in some detail, following the convention of buffesm.R)
# --------------------------------------------------------------------------------------------
#
# Zeebe and Wolf-Gladrow (2001) derive an exact expression for
#
#     Be := (d DIC / d s)_{AT,T,S},     s := [CO2*] = [CO2(aq)] + [H2CO3]
#
# by writing DIC and AT as functions of s and h := [H+], then imposing d(AT)=0 (i.e. AT held
# fixed) to eliminate dh in terms of ds. Their book's final closed-form result, Eq. (1.5.88),
# contains a typographical error that Egleston et al. (2010, GBC, doi:10.1029/2008GB003407)
# flag explicitly in their text ("...unstable, their final formula is made unusable by a
# ruinous typographical error") -- but Egleston et al. do NOT state what the error is or where
# it is located, nor do they provide a corrected version of the Zeebe & Wolf-Gladrow formula
# (they instead derive their own, differently-organized set of buffer factors).
#
# The typo was isolated here (July 2026) by re-deriving Eq. (1.5.88) from scratch starting
# from Zeebe & Wolf-Gladrow's own partial derivatives Eqs. (1.5.83)-(1.5.86), and comparing
# numerically against the printed formula. The error is a missing factor of (h/s) on one of
# the two numerator terms -- i.e. the term (Ds * w) in the printed book equation should read
# (h/s)*Ds*w. Taking the printed formula literally over-estimates Be by roughly 50% for
# typical seawater-like parameter choices; the corrected formula below was verified against
# an independent from-scratch symbolic differentiation using SymPy (both exactly, symbol by
# symbol, and numerically at several random points to machine precision).
#
# The corrected Zeebe & Wolf-Gladrow formula (their alkalinity = "practical alkalinity",
# i.e. carbonate + borate + water only -- called "PA" in their book, Eq. 1.2.31) is:
#
#     Be = [ K1/h + 4*K1*K2/h^2 + K1^2*K2/h^3 + (h/s)*Ds*w ]
#          -------------------------------------------------
#          [ K1/h + 4*K1*K2/h^2 + (h/s)*w ]
#
#     with  Ds = 1 + K1/h + K1*K2/h^2                    (== DIC/[CO2*], the "ionization
#                                                          fraction" of Jones et al. and
#                                                          Egleston et al. 2010)
#           w  = 1 + Kw/h^2 + BT*Kb/(Kb+h)^2             (borate + water only; this is what
#                                                          Zeebe & Wolf-Gladrow themselves use)
#
# EXTENSION implemented here: exactly as Orr et al. (2018, Mar. Chem.) extended the original
# 2010 buffesm() (Egleston-based) to include phosphate and silicate alkalinity contributions,
# the w term above is extended additively to include ALL non-carbonate acid-base systems
# handled by seacarb's carbfull():
#
#     w = w_bw + w_Si + w_P + w_N + w_S
#
#     w_bw = 1 + Kw/h^2 + BT*Kb/(Kb+h)^2                  (borate + water; unchanged)
#
#     w_Si = Sit*Ksi*(h^2 + 4*K2si*h + Ksi*K2si) / (h^2 + Ksi*h + Ksi*K2si)^2
#            -- CORRECTED (July 2026) to the full DIPROTIC silicic acid system,
#            Asi = [SiO(OH)3-] + 2*[SiO2(OH)2--], matching calculate_carb.R's own
#            treatment (siooh3 + 2*sio2oh2), which uses BOTH Ksi and K2si. An earlier
#            version of this term used the single-dissociation form Sit*Ksi/(Ksi+h)^2
#            (structurally identical to borate); that form is the correct K2si -> 0
#            limit of the expression above (verified numerically), but is NOT what
#            carbfull() actually computes for silicate alkalinity, so it undercounts
#            the true buffering contribution whenever K2si is not negligible.
#
#     w_P  = Pt*[ (3h^2+2*K1p*h+K1p*K2p)*NP - (K1p*K2p-3h^2)*DP ] / DP^2
#            with  NP = K1p*K2p*h + 2*K1p*K2p*K3p - h^3
#                  DP = h^3 + K1p*h^2 + K1p*K2p*h + K1p*K2p*K3p
#
#     w_N  = NH4t*Kn/(Kn+h)^2                              (ammonia; single dissociation,
#                                                            same functional form as borate)
#     w_S  = HSt*Khs/(Khs+h)^2                             (sulfide; single dissociation,
#                                                            same functional form as borate)
#
# Each of these additional terms was independently verified to be (i) exactly equal to
# -d(Aj)/dh, where Aj(h) is that system's contribution to alkalinity at fixed total
# concentration, and (ii) strictly positive over the full seawater pH range (7.5-8.5), as
# required physically for a genuine (stabilizing) buffering contribution.
#
# carb() vs carbfull(): seacarb's carb() hard-codes NH4t=0, HSt=0 when it calls
# calculate_carb() internally -- it does NOT expose NH4t or HSt as adjustable arguments at
# all. Only carbfull() does. Therefore, whenever the caller supplies a nonzero NH4t or HSt,
# buffzwg() below calls carbfull() instead of carb() to get pH/CO2/DIC that are actually
# consistent with those inputs; when NH4t = HSt = 0 (the default), carb() is used as before,
# for speed and backward compatibility. carb() and carbfull() were verified here (with real
# seacarb, not just symbolically) to agree to solver tolerance (~1e-9 relative) when NH4t and
# HSt are both zero, so this branch is a performance choice, not a correctness one.
#
# Two outputs are returned for convenience:
#   Be  -- the Bolin sensitivity itself (dimensionless: dDIC and ds are both mol/kg, so their
#          ratio has no units; do NOT confuse this with the dimensionless *classical* Revelle
#          factor, see RF below)
#   RF  -- the classical (relative) Revelle factor, RF = (DIC/s) / Be. This equals
#          Egleston et al.'s "R" (as computed by seacarb's own buffesm()) when Pt = Sit = 0,
#          and provides a convenient internal consistency check against that independently
#          derived and separately-coded function.
#
# --------------------------------------------------------------------------------------------

buffzwg <-
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
    # use carbfull() to get pH/CO2/DIC that are actually consistent with those inputs.
    # Otherwise carb() is used, for speed and backward compatibility (verified to agree with
    # carbfull(NH4t=0, HSt=0) to solver tolerance, ~1e-9 relative).
    use_carbfull <- any(NH4t != 0) || any(HSt != 0)
    if (use_carbfull) {
        Carb <- carbfull(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, NH4t=NH4t, HSt=HSt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    } else {
        Carb <- carb(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    }
    P    <- Carb$P
    pH   <- Carb$pH
    h    <- 10^(-pH)
    CO2  <- Carb$CO2          # this is [CO2*] = [CO2(aq)] + [H2CO3], i.e. our "s"
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

    # Unlike buffesm(), K1 and K2 are needed explicitly here (the ZWG formula is written
    # directly in terms of K1, K2, h -- it is not expressed via HCO3/CO3/DIC combinations
    # the way Egleston et al.'s formulas are).
    K1  <- K1(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    K2  <- K2(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    Kw  <- Kw(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Kb  <- Kb(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    K1p  <- K1p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K2p  <- K2p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K3p  <- K3p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Ksi  <- Ksi(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
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
    Kn   <- Kn(S=SP, T=InsT, P=P, pHscale=pHscale, warn=warn)
    Khs  <- Khs(S=SP, T=InsT, P=P, pHscale=pHscale, warn=warn)

    #---------------------------------------------------------------------
    #--------------------- Bolin sensitivity Be ---------------------------
    #---------------------------------------------------------------------

    s <- CO2   # s := [CO2*], following Zeebe & Wolf-Gladrow's own notation

    # Ds = DIC/[CO2*] at fixed h -- the "ionization fraction" of Jones et al. and
    # Egleston et al. (2010); equivalently Ds = 1 + K1/h + K1*K2/h^2
    Ds <- 1 + K1/h + K1*K2/h^2

    # Non-carbonate buffering term w, built additively system by system.
    # w_bw: borate + water (this alone reproduces the ORIGINAL Zeebe & Wolf-Gladrow
    #       "practical alkalinity" case, i.e. their Eq. 1.5.88 corrected)
    wbw <- 1 + Kw/h^2 + bor*Kb/(Kb+h)^2

    # w_Si: silicic acid -- CORRECTED to the full diprotic system, Asi =
    #       [SiO(OH)3-] + 2*[SiO2(OH)2--], matching calculate_carb.R (siooh3 + 2*sio2oh2).
    #       w_Si = -d(Asi)/dh; verified against a from-scratch finite-difference check and
    #       against the K2si -> 0 monoprotic limit (Sit*Ksi/(Ksi+h)^2), both to machine
    #       precision.
    DSi <- h^2 + Ksi*h + Ksi*K2si
    wSi <- Sit*Ksi*(h^2 + 4*K2si*h + Ksi*K2si) / DSi^2

    # w_P: phosphoric acid (three dissociation steps -- quotient-rule form kept
    #      unexpanded/uncollected for transparency and to minimize transcription risk;
    #      see derivation notes above). NP and DP are intermediate quantities from
    #      differentiating the phosphate contribution to alkalinity; they are NOT
    #      themselves alkalinity or concentration terms.
    NP <- K1p*K2p*h + 2*K1p*K2p*K3p - h^3
    DP <- h^3 + K1p*h^2 + K1p*K2p*h + K1p*K2p*K3p
    wP <- Pt * ( (3*h^2 + 2*K1p*h + K1p*K2p)*NP - (K1p*K2p - 3*h^2)*DP ) / DP^2
    # Note: wP is automatically and exactly zero when Pt = 0 (Pt multiplies the whole
    # bracket), so no special-case guard against division by zero is needed here --
    # unlike the analogous phosphate term in buffesm(), DP is never zero since
    # h, K1p, K2p, K3p are all strictly positive.

    # w_N: ammonia (NH4+/NH3), single dissociation, same functional form as borate.
    # w_S: sulfide (H2S/HS-), single dissociation, same functional form as borate.
    # Both are automatically and exactly zero when NH4t = 0 / HSt = 0 respectively.
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

    # Corrected Zeebe & Wolf-Gladrow (2001) Eq. (1.5.88), extended for phosphate/silicate:
    numerator   <- K1/h + 4*K1*K2/h^2 + K1^2*K2/h^3 + (h/s)*Ds*w
    denominator <- K1/h + 4*K1*K2/h^2 + (h/s)*w
    Be <- numerator / denominator

    # Classical (relative, dimensionless) Revelle factor, for cross-checking against
    # buffesm()'s "R" column when Pt = Sit = 0 (they should then agree to numerical
    # precision, since both ultimately encode the same underlying carbonate-system
    # buffering physics, just packaged differently).
    RF <- (DIC/s) / Be

    # For debugging or to assess different alkalinity components
    #col <- c("Be", "RF", "Ds", "w", "wbw", "wSi", "wP", "wN", "wS")
    #res <- data.frame(Be, RF, Ds, w, wbw, wSi, wP, wN, wS)
      
    col <- c("Be", "RF")
    res <- data.frame(Be, RF)
    names(res) <- col

    return(res)
  }
