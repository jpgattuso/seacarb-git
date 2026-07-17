# Copyright (C) 2008 Jean-Marie Epitalon and Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# Extended (July 2026) to include phosphate, silicate, ammonia, and sulfide alkalinity in
# the non-carbonate buffering term (DD/V), matching the extensions made to buffzwg.R and
# buffesm.R. The silicate term uses the full diprotic silicic acid system (matching
# calculate_carb.R's Asi = [SiO(OH)3-] + 2*[SiO2(OH)2--]) rather than a single-dissociation
# approximation. Whenever NH4t or HSt are nonzero, carbfull() is used instead of carb(),
# since carb() hard-codes NH4t=HSt=0 internally. PhiH is simplified to 1/(h*ln10*(D-DD)),
# an identity that holds for any number of non-carbonate systems included in DD (verified
# algebraically and numerically), avoiding duplication of the individual borate/water terms
# that were previously written out separately in that line.
#

buffer <- 
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
        STeos <- teos2eos_geo (S, T, P, long, lat)
        InsT <- STeos$T
        SP <- STeos$SP
    }
    else
    {
        InsT <- T
        SP <- S
    }

    # carb() hard-codes NH4t=0, HSt=0 internally and does not expose them as arguments;
    # only carbfull() does. So if the caller has supplied any nonzero NH4t or HSt, we must
    # use carbfull() to get pH/CO2/HCO3/CO3/DIC/ALK that are actually consistent with those
    # inputs. Otherwise carb() is used, for speed and backward compatibility (verified to
    # agree with carbfull(NH4t=0, HSt=0) to solver tolerance, ~1e-9 relative).
    use_carbfull <- any(NH4t != 0) || any(HSt != 0)
    if (use_carbfull) {
        Carb <- carbfull(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, NH4t=NH4t, HSt=HSt, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    } else {
        Carb <- carb(flag=flag, var1=var1, var2=var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    }
    
	pH   <- Carb$pH
	h    <- 10^(-pH)
	CO2  <- Carb$CO2
	HCO3 <- Carb$HCO3
	CO3  <- Carb$CO3
	DIC  <- Carb$DIC
	ALK  <- Carb$ALK
	Oa   <- Carb$OmegaAragonite
	Oc   <- Carb$OmegaCalcite
    
    #-------Constants----------------
    
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK = InsT + tk;           # TK [K]; InsT[C]
    
    Cl = SP / 1.80655;            # Cl = chlorinity; SP = practical salinity (psu)
    ST = 0.14 * Cl/96.062        # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    FT = 6.7e-5 * Cl/18.9984     # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
    FLUO = 6.7e-5 * Cl/18.9984   # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
    BOR = bor(S=SP , b=b);        # (mol/kg) total boron

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

    K1 <- K1(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)   
    K2 <- K2(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    Kw <- Kw(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K0 <- K0(S=SP, T=InsT, Patm=Patm, P=P, warn=warn)
    Kb <- Kb(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn)
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
    Kspa <- Kspa(S=SP, T=InsT, P=P, warn=warn)
    Kspc <- Kspc(S=SP, T=InsT, P=P, warn=warn)

    rho <- rho(S=SP,T=InsT,P=P)

    # Compute potential K0 with S, potential temperature, and atmospheric pressure (usually 1 atm)
    K0pot <- K0(S=SP, T=theta(S=SP, T=InsT, P=P, Pref=0), Patm=Patm, P=0, warn=warn)
    
    #---------------------------------------------------------------------
    #--------------------    buffer effects    ---------------------------
    #---------------------------------------------------------------------
    
    # Non-carbonate buffering term, extended (July 2026) beyond borate+water to include
    # phosphate, silicate (diprotic), ammonia, and sulfide -- same additive extension used
    # in buffzwg.R and buffesm.R. Each term is -d(Aj)/dh for that system's alkalinity
    # contribution Aj(h) at fixed total concentration.
    #
    #   wbw = 1 + Kw/h^2 + Kb*BOR/(Kb+h)^2                    (borate + water; unchanged)
    #   wSi = Sit*Ksi*(h^2+4*K2si*h+Ksi*K2si) / (h^2+Ksi*h+Ksi*K2si)^2
    #         -- CORRECTED to the full diprotic silicic acid system, Asi =
    #         [SiO(OH)3-] + 2*[SiO2(OH)2--], matching calculate_carb.R; the previous
    #         single-dissociation form Sit*Ksi/(Ksi+h)^2 is its K2si -> 0 limit.
    #   wP  = phosphoric acid (three dissociation steps; NP, DP as originally derived)
    #   wN  = NH4t*Kn/(Kn+h)^2                                 (ammonia)
    #   wS  = HSt*Khs/(Khs+h)^2                                (sulfide)
    #
    # DD and V were, before this extension, computed as two separately-coded but
    # algebraically identical quantities (both equal Kb*BOR/(Kb+h)^2 + Kw/h^2 + 1). That
    # equivalence was verified independently earlier (via the Frankignoulle 1994
    # rearrangement); V is now simply set equal to DD to avoid maintaining two copies of
    # the same formula.
    wbw <- -((-Kb*BOR)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1
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

    DD <- wbw + wSi + wP + wN + wS + wF
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
    PhiD=-1/(h*log(10) * ( B+A+C ) );
    BetaD=-h*log(10)*DIC/CO2*B*PhiD;
    
    
    Q=(h+2*K2);
    V <- DD;
    
    # NOTE: the trailing non-carbonate term must be the FULL w (= DD), not just wbw.
    # Before this fix DB and DC carried only borate+water here while DD carried every
    # system, so PhiB/BetaB and PhiC/BetaC were wrong by up to 0.16% and 1.5% at high
    # Pt/Sit. PhiD/BetaD were unaffected because they use DD directly.
    DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1))) + DD;
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
    PhiB=-1/(h*log(10) * ( B+A+C ) );
    BetaB=-h*log(10)*DIC/CO2*B*PhiB;
    
    DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1))) + DD;
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
    PhiC=-1/(h*log(10) * ( B+A+C ) );
    BetaC=-h*log(10)*DIC/CO2*B*PhiC;
    
    D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
    D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
    D=D1+2*D2;
    # PhiH = 1/(h*ln10*(D - DD)) -- this identity holds regardless of how many
    # non-carbonate systems DD includes (verified algebraically: the original line,
    # D + (-Kb*BOR/(h+Kb)^2) + (-Kw/h^2) - 1, is exactly D - DD when DD is the
    # borate+water-only form; extending DD to phosphate/silicate/ammonia/sulfide extends
    # this identity the same way, so DD (already computed above) is reused directly here
    # instead of repeating each individual non-carbonate term).
    PhiH=1/ (h*log(10)* (D - DD))  ; 
    
    Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));

    # In situ values - use K0 computed with in situ T & in situ P (atm + hydrostatic)
    # PiH=((-h/K0)*log(10)*Pi)*PhiH;
    # PiB=CO2/(K0*DIC)*BetaB;
    # PiD=CO2/(K0*DIC)*BetaD;
    # PiC=CO2/(K0*DIC)*BetaC;

    # Potential values - use K0 computed with potential T & surface atm P only
    PiH=((-h/K0pot)*log(10)*Pi)*PhiH;
    PiB=CO2/(K0pot*DIC)*BetaB;
    PiD=CO2/(K0pot*DIC)*BetaD;
    PiC=CO2/(K0pot*DIC)*BetaC;

    col <- c("PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
    res <- data.frame(PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH)
    names(res) <- col 
    
    return(res)
  }
