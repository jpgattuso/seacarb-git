# Copyright (C) 2008 Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye and Jean-Pierre Gattuso  
# with a most valuable contribution of Bernard Gentili
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
#
# This function is the common implementation of  three public funtions : carb(), carbb() and carbfull()
# As such, it accepts as input the parameters of all three functions.
#
# The last parameter, "envir", is an optional environment. When given, it is the environment of the caller funtion.
# This environment may contain redefinition of some public seacarb functions such as K0(), K1(), K2(), ....

calculate_carb<-
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, NH4t=0, HSt=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", badd=0,
         warn="y", eos="eos80", long=1.e20, lat=1.e20, fullresult=FALSE, envir=NULL){
    n <- max(length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), length(ks), length(b))
    if(length(flag)!=n){flag <- rep(flag[1],n)}
    if(length(var1)!=n){var1 <- rep(var1[1],n)}
    if(length(var2)!=n){var2 <- rep(var2[1],n)}
    if(length(S)!=n){S <- rep(S[1],n)}
    if(length(T)!=n){T <- rep(T[1],n)}
    if(length(Patm)!=n){Patm <- rep(Patm[1],n)}
    if(length(P)!=n){P <- rep(P[1],n)}
    if(length(Pt)!=n){Pt <- rep(Pt[1],n)}
    if(length(Sit)!=n){Sit <- rep(Sit[1],n)}
    if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
    if(length(kf)!=n){kf <- rep(kf[1],n)}
    if(length(ks)!=n){ks <- rep(ks[1],n)}
    if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
    if(length(b)!=n){b <- rep(b[1],n)}
    if(length(gas)!=n){gas <- rep(gas[1],n)}
    if (fullresult) {
        if(length(NH4t)!=n){NH4t <- rep(NH4t[1],n)}
        if(length(HSt)!=n){HSt <- rep(HSt[1],n)}
    } else {
        NH4t <- rep(0,n)
        HSt <- rep(0,n)
    }
    
    # if the concentrations of total silicate, total phosphate, total nitrite,
    # total ammonium, and total hydrogen sulfide are NA, they are set at 0
    Sit[is.na(Sit)] <- 0
    Pt[is.na(Pt)] <- 0
    NH4t[is.na(NH4t)] <- 0
    HSt[is.na(HSt)] <- 0
    
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
 
    #-------Constants----------------   
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])  
    TK = InsT + tk;           # TK [K]; T[C]
    
    #---- issues de equic----
    Cl = SP / 1.80655;             # Cl = chlorinity; S = salinity (per mille)
    ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    FLUO = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
 
    # get functions that compute K0, K1, K2, Kb, Kw, Kspa, Kspc and Boron concentration
    # --------------------------------------------------------------------------------------------------------------------
    
    if (! is.null(envir)) 
    {
        # These functions may, in certain cases, have a replacement function (when computing derivatives, namely)
        # We get the original one or, when it exists, the replacement one by accessing them through the "environment" given as input parameter
        
        K0_fct <-get("K0", envir=envir)
        K1_fct <-get("K1", envir=envir)
        K2_fct <-get("K2", envir=envir)
        Kw_fct <-get("Kw", envir=envir) 
        Kb_fct <-get("Kb", envir=envir)
        Kspa_fct <-get("Kspa", envir=envir)
        Kspc_fct <-get("Kspc", envir=envir)
        bor_fct <-get("bor", envir=envir)
    }
    else
    {
        # Get the original functions from Seacarb
        K0_fct <-K0
        K1_fct <-K1
        K2_fct <-K2
        Kw_fct <-Kw
        Kb_fct <-Kb
        Kspa_fct <-Kspa
        Kspc_fct <-Kspc
        bor_fct <-bor
    }
   
    
    #---------------------------------------------------------------------
    #-------------- compute K's and Boron ---------------------
    #---------------------------------------------------------------------
    
    BOR = bor_fct(S=S , b=b) + badd;         # (mol/kg) total boron + boron added
    
    # Ks (free pH scale) at zero pressure and given pressure
    Ks_P0 <- Ks(S=SP, T=InsT, P=0, ks=ks, warn=warn)
    Ks    <- Ks(S=SP, T=InsT, P=P, ks=ks, warn=warn)
    
    # Kf on free pH scale
    Kff_P0 <- Kf(S=SP, T=InsT, P=0, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
    Kff <- Kf(S=SP, T=InsT, P=P, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
    # Kf on given pH scale
    Kf <- Kf(S=SP, T=InsT, P=P, pHscale=pHscale, kf=kf, Ks_P0, Ks, warn=warn)
    
    # Conversion factor from total to SWS pH scale at zero pressure
    ktotal2SWS_P0 <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks_P0,Kff=Kff_P0)$ktotal2SWS

    # Conversion factor from SWS to chosen pH scale
    conv <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks,Kff=Kff, warn=warn)
    kSWS2chosen <- rep(1.,n)
    kSWS2chosen [pHscale == "T"] <- conv$kSWS2total [pHscale == "T"]
    kSWS2chosen [pHscale == "F"] <- conv$kSWS2free [pHscale == "F"]  

    K1 <- K1_fct(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)   
    K2 <- K2_fct(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    Kw <- Kw_fct(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Kb <- Kb_fct(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn)
    K1p <- K1p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K2p <- K2p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    K3p <- K3p(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Ksi <- Ksi(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn)
    Kspa <- Kspa_fct(S=SP, T=InsT, P=P, warn=warn)
    Kspc <- Kspc_fct(S=SP, T=InsT, P=P, warn=warn)
    if (fullresult) {
        Kn <- Kn(S=S, T=T, P=P, pHscale=pHscale)
        Khs <- Khs(S=S, T=T, P=P, pHscale=pHscale)
        K2si <- K2si(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0)
    } else {
        Kn <- rep(0,n)
        Khs <- rep(0,n)
        K2si <- rep(0,n)
    }
    
    rho <- rho(S=SP,T=InsT,P=P)

    # Compute "standard" K0 with S, in situ T, and atmospheric pressure
    K0 <- K0_fct(S=SP, T=InsT, Patm=Patm, P=0, warn=warn)                         
    # Compute potential K0 with S, potential temperature, and atmospheric pressure (usually 1 atm)
    K0pot <- K0_fct(S=SP, T=theta(S=SP, T=InsT, P=P, Pref=0), Patm=Patm, P=0, warn=warn)
    # Compute in situ K0 with S, in situ temperature, and total pressure pressure (atmospheric + hydrostatic)
    K0insitu <- K0_fct(S=SP, T=InsT, Patm=Patm, P=P, warn=warn)
 
    #------------------------------------------------------------------#
    #------------------------------------------------------------------#
    #                            VARIABLES                             #
    #------------------------------------------------------------------#
    #------------------------------------------------------------------#
    
    # flag = 1      pH-CO2 given
    # flag = 2      CO2-HCO3 given
    # flag = 3      CO2-CO3 given
    # flag = 4      CO2-ALK given
    # flag = 5      CO2-DIC given
    # flag = 6      pH and HCO3 given
    # flag = 7      pH and CO3 given
    # flag = 8      pH and ALK given
    # flag = 9      pH and DIC given
    # flag = 10     HCO3 and CO3 given
    # flag = 11     HCO3 and ALK given
    # flag = 12     HCO3 and DIC given
    # flag = 13     CO3 and ALK given
    # flag = 14     CO3 and DIC given
    # flag = 15     ALK and DIC given
    # flag = 21     pCO2-pH given
    # flag = 22     pCO2-HCO3 given
    # flag = 23     pCO2-CO3 given
    # flag = 24     pCO2-ALK given
    # flag = 25     pCO2-DIC given
    
    # Initialise output vectors
    H   <- rep(NA, n)
    PH   <- rep(NA, n)
    CO2  <- rep(NA, n)
    pCO2 <- rep(NA, n)
    fCO2 <- rep(NA, n)
    HCO3 <- rep(NA, n)
    CO3  <- rep(NA, n)
    DIC  <- rep(NA, n)
    ALK  <- rep(NA, n)

    # ------------ case 1.) PH and CO2 given ----
    # Indices of flag elements where flag = 1
    i_flag_1 <- which (flag == 1)    
    PH[i_flag_1]   <- var1[i_flag_1]
    CO2[i_flag_1]  <- var2[i_flag_1]
    h <- 10^(-PH[i_flag_1])
    fCO2[i_flag_1] <- CO2[i_flag_1] / K0[i_flag_1]
    HCO3[i_flag_1] <- (K1[i_flag_1] * CO2[i_flag_1]) / h
    CO3[i_flag_1]  <- (K2[i_flag_1] * HCO3[i_flag_1]) / h
    DIC[i_flag_1]  <- CO2[i_flag_1] + HCO3[i_flag_1] + CO3[i_flag_1]
    H[i_flag_1] <- h
    
    # ------------ case 2.) CO2 and HCO3 given ----
    # Indices of flag elements where flag = 2
    i_flag_2 <- which (flag == 2)     
    CO2[i_flag_2]  <- var1[i_flag_2]
    HCO3[i_flag_2] <- var2[i_flag_2]
    fCO2[i_flag_2] <- CO2[i_flag_2] / K0[i_flag_2]
    h <- K0[i_flag_2] * K1[i_flag_2] * fCO2[i_flag_2] / HCO3[i_flag_2]
    CO3[i_flag_2]  <- K0[i_flag_2] * K1[i_flag_2] * K2[i_flag_2] * fCO2[i_flag_2] / (h^2)
    DIC[i_flag_2]  <- CO2[i_flag_2] + HCO3[i_flag_2] + CO3[i_flag_2]
    PH[i_flag_2]   <- -log10(h)
    H[i_flag_2] <- h

    # ------------ case 3.) CO2 and CO3 given ----
    # Indices of flag elements where flag = 3
    i_flag_3 <- which (flag == 3)     
    CO2[i_flag_3]  <- var1[i_flag_3]
    CO3[i_flag_3]  <- var2[i_flag_3]
    fCO2[i_flag_3] <- CO2[i_flag_3] / K0[i_flag_3]
    k1co2 <- K1[i_flag_3] * CO2[i_flag_3]
    h <- sqrt((K2[i_flag_3]*k1co2) / CO3[i_flag_3])
    HCO3[i_flag_3] <- k1co2 / h
    DIC[i_flag_3]  <- CO2[i_flag_3] + HCO3[i_flag_3] + CO3[i_flag_3]
    PH[i_flag_3]   <- -log10(h)
    H[i_flag_3] <- h

    # ------------ case 4.) CO2 and ALK given ----
    # Indices of flag elements where flag = 4
    i_flag_4 <- which (flag == 4)     
    CO2[i_flag_4] <- var1[i_flag_4]
    ALK[i_flag_4] <- var2[i_flag_4]
    
    # Calculate [H+] from [CO2] and total alk
    h <- rep(NA, length(i_flag_4))
    j <- 1
    for(i in (i_flag_4))
    {
        # Create a list containing all the dissociation constants
        dissoc = list(K1_DIC=K1[i], K2_DIC=K2[i], K_BT=Kb[i], 
                        K1_PO4 = K1p[i], K2_PO4 = K2p[i], K3_PO4 = K3p[i], K_Sil = Ksi[i],
                        K2_Sil = K2si[i], K_NH4 = Kn[i], K_H2S = Khs[i],
                        K_HSO4 = Ks[i], K_HF = Kff[i], K_H2O = Kw[i])

        # use SolveSAPHE v2.0
        result <- solve_pH_from_AT(ALK[i], CO2[i], BOR[i], Pt[i], Sit[i], NH4t[i], HSt[i], 
                                    ST[i], FLUO[i], pHscale[i], "CO2", FALSE, dissoc)
        h[j] <- result[1]
        j <- j + 1
    }
 
    fCO2[i_flag_4] <- CO2[i_flag_4] / K0[i_flag_4]
    HCO3[i_flag_4] <- K1[i_flag_4]*CO2[i_flag_4]/h
    CO3[i_flag_4]  <- K2[i_flag_4]*HCO3[i_flag_4]/h
    PH[i_flag_4]   <- -log10(h)
    H[i_flag_4] <- h
    DIC[i_flag_4]  <- CO2[i_flag_4] + HCO3[i_flag_4] + CO3[i_flag_4]


    # ------------ case 5.) CO2 and DIC given ----
    # Indices of flag elements where flag = 5
    i_flag_5 <- which (flag == 5)     
    CO2[i_flag_5]  <- var1[i_flag_5]
    DIC[i_flag_5]  <- var2[i_flag_5]
    fCO2[i_flag_5] <- CO2[i_flag_5] / K0[i_flag_5]
    a <- K1[i_flag_5] * K2[i_flag_5] * CO2[i_flag_5]
    b <- K1[i_flag_5] * CO2[i_flag_5]
    c <- CO2[i_flag_5] - DIC[i_flag_5]
    D <- b*b - 4*a*c
    h <- (2*a) / (sqrt(D)-b)
    HCO3[i_flag_5] <- K0[i_flag_5] * K1[i_flag_5] * fCO2[i_flag_5] / h
    CO3[i_flag_5] <- DIC[i_flag_5] - CO2[i_flag_5] - HCO3[i_flag_5]
    PH[i_flag_5] <- -log10(h)
    H[i_flag_5] <- h

    # ------------ case 6.) PH and HCO3 given ----
    # Indices of flag elements where flag = 6
    i_flag_6 <- which (flag == 6)     
    PH[i_flag_6] <- var1[i_flag_6]
    HCO3[i_flag_6] <- var2[i_flag_6]
    h <- 10^(-PH[i_flag_6])
    CO2[i_flag_6] <- (HCO3[i_flag_6] * h)/K1[i_flag_6]
    CO3[i_flag_6] <- K2[i_flag_6] * HCO3[i_flag_6] /h
    DIC[i_flag_6] <- CO2[i_flag_6] + HCO3[i_flag_6] + CO3[i_flag_6]
    fCO2[i_flag_6] <- CO2[i_flag_6]/K0[i_flag_6]
    H[i_flag_6] <- h

    # ------------ case 7.) PH and CO3 given ---- 
    # Indices of flag elements where flag = 7
    i_flag_7 <- which (flag == 7)     
    PH[i_flag_7] <- var1[i_flag_7]
    CO3[i_flag_7] <- var2[i_flag_7]
    h <- 10^(-PH[i_flag_7])
    HCO3[i_flag_7] <- CO3[i_flag_7] * h/K2[i_flag_7]
    CO2[i_flag_7]  <- HCO3[i_flag_7] * h/K1[i_flag_7]
    fCO2[i_flag_7] <- CO2[i_flag_7]/K0[i_flag_7]
    DIC[i_flag_7]  <- CO2[i_flag_7] + HCO3[i_flag_7] + CO3[i_flag_7]
    H[i_flag_7] <- h

    # ------------ case 8.) PH and ALK given ----
    # Indices of flag elements where flag = 8
    i_flag_8 <- which(flag == 8)     
    PH[i_flag_8]  <- var1[i_flag_8]
    ALK[i_flag_8] <- var2[i_flag_8] 
    h <- 10^(-PH[i_flag_8])
    H[i_flag_8] <- h

    ## calculate Hfree anf Htot
    hfree <- rep(NA, length(i_flag_8))
    htot  <- rep(NA, length(i_flag_8))    
    sc <- pHscale[i_flag_8] 
    st <- ST[i_flag_8]
    ks  <- Ks[i_flag_8]
    fluo <- FLUO[i_flag_8]
    kff  <- Kff[i_flag_8]
    
    # Where pHscale=="F", pHscale = free scale
    i_sc_F <- which (sc == "F") 
    hfree[i_sc_F] <- h[i_sc_F]  
    htot[i_sc_F]  <- h[i_sc_F] / (1+st[i_sc_F]/ks[i_sc_F]) 
    # Where pHscale=="T", pHscale = total scale
    i_sc_T <- which (sc == "T") 
    hfree[i_sc_T] <- h[i_sc_T] * (1+st[i_sc_T]/ks[i_sc_T])
    htot[i_sc_T]  <- h[i_sc_T]
    # Where pHscale=="SWS", pHscale = SW scale
    i_sc_S <- which (sc == "SWS") 
    hfree[i_sc_S] <- h[i_sc_S] * (1 + st[i_sc_S]/ks[i_sc_S] + fluo[i_sc_S]/kff[i_sc_S])
    htot[i_sc_S]  <- hfree[i_sc_S] / (1+st[i_sc_S]/ks[i_sc_S])

    # Calculate some invariable components of total alkalinity
    boh4   <- BOR[i_flag_8] / (1+h/Kb[i_flag_8])
    oh     <- Kw[i_flag_8]/h
    temp   <- (h^3 + K1p[i_flag_8]*h^2 + K1p[i_flag_8]*K2p[i_flag_8]*h +
              K1p[i_flag_8]*K2p[i_flag_8]*K3p[i_flag_8]) 
    h3po4  <- Pt[i_flag_8]*h^3 / temp
    hpo4   <- Pt[i_flag_8]*K1p[i_flag_8]*K2p[i_flag_8]*h / temp
    po4    <- Pt[i_flag_8]*K1p[i_flag_8]*K2p[i_flag_8]*K3p[i_flag_8] / temp
    hso4   <- st/(1+ks/hfree)
    hf     <- fluo / (1+Kf[i_flag_8]/htot)

    if (fullresult) {
        # adapted to include second dissociation constant of silicate
        siooh3 <- Sit[i_flag_8]*h*Ksi[i_flag_8] / 
        (h*h + Ksi[i_flag_8]*h + Ksi[i_flag_8]*K2si[i_flag_8])
        sio2oh2 <- Sit[i_flag_8]*Ksi[i_flag_8]*K2si[i_flag_8] / 
        (h*h + Ksi[i_flag_8]*h + Ksi[i_flag_8]*K2si[i_flag_8])
        nh3 <- NH4t[i_flag_8]/(1+h/Kn[i_flag_8])
        hs <- HSt[i_flag_8]/(1+h/Khs[i_flag_8])

        # Sum of these components (partial alkalinity)
        alk_p <- boh4+oh+hpo4+2*po4+siooh3+2*sio2oh2+nh3+hs-hfree-hso4-hf-h3po4
    }
    else {
        siooh3 <- Sit[i_flag_8] / (1+h/Ksi[i_flag_8])
        # Sum of these components (partial alkalinity)
        alk_p <- boh4+oh+hpo4+2*po4+siooh3-hfree-hso4-hf-h3po4
    }
    
    # As [HCO3-] and  [CO3--] are both linearly dependent on DIC
    #    hco3  = co2*K1/h = DIC/(1+K1/h+K1*K2/(h^2)) * K1/h
    #   2*co3  = 2*hco3*K2/h = DIC/(1+K1/h+K1*K2/(h^2)) * 2*K1*K2/h^2
    # carbonate alk ([HCO3-] + 2*[CO3--]) is also linearly dependent on DIC
    # Calculate DIC from carbonate alk (total alk - partial alk)
    temp <- h*h + K1[i_flag_8]*h + K1[i_flag_8]*K2[i_flag_8]
    DIC[i_flag_8] <- (ALK[i_flag_8]-alk_p) * temp / 
         (h*K1[i_flag_8] + 2*K1[i_flag_8]*K2[i_flag_8])

    CO2[i_flag_8]  <- DIC[i_flag_8] / 
        (1 + K1[i_flag_8]/h + K1[i_flag_8]*K2[i_flag_8]/(h^2))
    HCO3[i_flag_8] <- CO2[i_flag_8]*K1[i_flag_8]/h
    CO3[i_flag_8]  <- HCO3[i_flag_8]*K2[i_flag_8]/h
    fCO2[i_flag_8] <- CO2[i_flag_8]/K0[i_flag_8]
    
    # ------------ case 9.) PH and DIC given ----
    # Indices of flag elements where flag = 9
    i_flag_9 <- which (flag == 9)     
    PH[i_flag_9]  <- var1[i_flag_9]
    DIC[i_flag_9] <- var2[i_flag_9]
    h <- 10^(-PH[i_flag_9])
    H[i_flag_9] <- h
    temp <- h*h + K1[i_flag_9]*h + K1[i_flag_9]*K2[i_flag_9]
    HCO3[i_flag_9] <- (DIC[i_flag_9]*K1[i_flag_9]*h) / temp
    CO3[i_flag_9]  <- (DIC[i_flag_9]*K1[i_flag_9]*K2[i_flag_9]) / temp
    CO2[i_flag_9]  <- h*HCO3[i_flag_9]/K1[i_flag_9]
    fCO2[i_flag_9] <- CO2[i_flag_9]/K0[i_flag_9]
    
    # ------------ case 10.) HCO3 and CO3 given ----
    # Indices of flag elements where flag = 10
    i_flag_10 <- which (flag == 10)     
    HCO3[i_flag_10] <- var1[i_flag_10]
    CO3[i_flag_10]  <- var2[i_flag_10]
    h <- K2[i_flag_10]*HCO3[i_flag_10]/CO3[i_flag_10]
    CO2[i_flag_10]  <- h*HCO3[i_flag_10]/K1[i_flag_10]
    DIC[i_flag_10]  <- CO2[i_flag_10] + HCO3[i_flag_10] + CO3[i_flag_10]
    fCO2[i_flag_10] <- CO2[i_flag_10]/K0[i_flag_10]
    PH[i_flag_10] <- -log10(h)
    H[i_flag_10] <- h

    # ------------ case 11.) HCO3 and ALK given ----
    # Indices of flag elements where flag = 11
    i_flag_11 <- which (flag == 11)     
    HCO3[i_flag_11] <- var1[i_flag_11]
    ALK[i_flag_11]  <- var2[i_flag_11]

    # Calculate [H+] from [HCO3] and total alk
    h <- rep(NA, length(i_flag_11))
    j <- 1 
    for(i in (i_flag_11))
    {
        # Create a list containing all the dissociation constants
        dissoc = list(K1_DIC=K1[i], K2_DIC=K2[i], K_BT=Kb[i], 
                        K1_PO4 = K1p[i], K2_PO4 = K2p[i], K3_PO4 = K3p[i], K_Sil = Ksi[i], 
                        K2_Sil = K2si[i], K_NH4 = Kn[i], K_H2S = Khs[i],
                        K_HSO4 = Ks[i], K_HF = Kff[i], K_H2O = Kw[i])

        # Use SolveSAPHE PACKAGE
        result <- solve_pH_from_AT(ALK[i], HCO3[i], BOR[i], Pt[i], Sit[i], NH4t[i], HSt[i],
                                    ST[i], FLUO[i], pHscale[i], "HCO3", FALSE, dissoc)
        h[j] <- result[1]
        j <- j + 1
    }
    
    CO2[i_flag_11] <- h*HCO3[i_flag_11]/K1[i_flag_11]
    CO3[i_flag_11] <- K2[i_flag_11]*HCO3[i_flag_11]/h
    DIC[i_flag_11] <- CO2[i_flag_11] + HCO3[i_flag_11] + CO3[i_flag_11]
    PH[i_flag_11] <- -log10(h)
    H[i_flag_11] <- h
    fCO2[i_flag_11] <- CO2[i_flag_11]/K0[i_flag_11]

    # ------------ case 12.) HCO3 and DIC given
    # Indices of flag elements where flag = 12
    i_flag_12 <- which (flag == 12)     
    HCO3[i_flag_12] <- var1[i_flag_12]
    DIC[i_flag_12]  <- var2[i_flag_12]
    a <- HCO3[i_flag_12]
    b <- K1[i_flag_12]*(HCO3[i_flag_12]-DIC[i_flag_12])
    c <- K1[i_flag_12]*K2[i_flag_12]*HCO3[i_flag_12]
    D <- b*b - 4*a*c
    h <- (-b-sqrt(D))/(2*a)
    CO2[i_flag_12]  <- h*HCO3[i_flag_12]/K1[i_flag_12]
    CO3[i_flag_12]  <- K2[i_flag_12]*HCO3[i_flag_12]/h
    fCO2[i_flag_12] <- CO2[i_flag_12]/K0[i_flag_12]
    PH[i_flag_12] <- -log10(h)
    H[i_flag_12] <- h

    # ------------ case 13.) CO3 and ALK given ----
    # Indices of flag elements where flag = 13
    i_flag_13 <- which (flag == 13)     
    CO3[i_flag_13] <- var1[i_flag_13]
    ALK[i_flag_13] <- var2[i_flag_13]

    h <- rep(NA, length(i_flag_13))
    j <- 1 
    for(i in (i_flag_13))
    {
        # Create a list containing all the dissociation constants
        dissoc = list(K1_DIC=K1[i], K2_DIC=K2[i], K_BT=Kb[i], 
                        K1_PO4 = K1p[i], K2_PO4 = K2p[i], K3_PO4 = K3p[i], K_Sil = Ksi[i], 
                        K2_Sil = K2si[i], K_NH4 = Kn[i], K_H2S = Khs[i],
                        K_HSO4 = Ks[i], K_HF = Kff[i], K_H2O = Kw[i])

        # Use solveSAPHE package
        result <- solve_pH_from_AT(ALK[i], CO3[i], BOR[i], Pt[i], Sit[i], NH4t[i], HSt[i],
                                    ST[i], FLUO[i], pHscale[i], "CO3", FALSE, dissoc)

        if (result[2] == -1) # only one solution for pH
        {                
            h[j] <- result[1]
        }
        else    # two solutions for pH
        {
            warning ("Given Alk and [CO3--], there were 2 solutions for pH. The one closest to pH7 was chosen !")
            # look for the solution closest to ph7
            h1 <- result[1]
            h2 <- result[2]
            if (abs(log10(h1) + 7) < abs(log10(h2) + 7)) {
                i_close7 <- 1
            } else {
                i_close7 <- 2
            }
            h[j] <- result[i_close7]
        }
        j <- j + 1
    }
    
    HCO3[i_flag_13] <- h*CO3[i_flag_13]/K2[i_flag_13]
    CO2[i_flag_13]  <- h*HCO3[i_flag_13]/K1[i_flag_13]
    fCO2[i_flag_13] <- CO2[i_flag_13]/K0[i_flag_13]
    DIC[i_flag_13]  <- HCO3[i_flag_13] + CO2[i_flag_13] + CO3[i_flag_13]
    PH[i_flag_13] <- -log10(h)
    H[i_flag_13] <- h

    # ------------ case 14.) CO3 and DIC given ----
    # Indices of flag elements where flag = 14
    i_flag_14 <- which (flag == 14)     
    CO3[i_flag_14] <- var1[i_flag_14]
    DIC[i_flag_14] <- var2[i_flag_14]

    a <- CO3[i_flag_14]
    b <- K1[i_flag_14]*CO3[i_flag_14]
    c <- K1[i_flag_14] * K2[i_flag_14] * (CO3[i_flag_14]-DIC[i_flag_14])
    D <- b*b - 4*a*c
    h <- (-b+sqrt(D))/(2*a)
    HCO3[i_flag_14] <- h*CO3[i_flag_14]/K2[i_flag_14]
    CO2[i_flag_14]  <- h*HCO3[i_flag_14]/K1[i_flag_14]
    fCO2[i_flag_14] <- CO2[i_flag_14]/K0[i_flag_14]
    PH[i_flag_14] <- -log10(h)
    H[i_flag_14] <- h
    
    # ------------ case 15.) ALK and DIC given ----
    # Indices of flag elements where flag = 15
    i_flag_15 <- which(flag == 15)     
    ALK[i_flag_15] <- var1[i_flag_15]
    DIC[i_flag_15] <- var2[i_flag_15]

    h <- rep(NA, length(i_flag_15))
    j <- 1 
    for(i in (i_flag_15))
    {
        # Create a list containing all the dissociation constants
        dissoc = list(K1_DIC=K1[i], K2_DIC=K2[i], K_BT=Kb[i], 
                        K1_PO4 = K1p[i], K2_PO4 = K2p[i], K3_PO4 = K3p[i], K_Sil = Ksi[i], 
                        K2_Sil = K2si[i], K_NH4 = Kn[i], K_H2S = Khs[i],
                        K_HSO4 = Ks[i], K_HF = Kff[i], K_H2O = Kw[i])

        # use SolveSAPHE v2.0
        result <- solve_pH_from_AT(ALK[i], DIC[i], BOR[i], Pt[i], Sit[i], NH4t[i], HSt[i],
                                    ST[i], FLUO[i], pHscale[i], "DIC", FALSE, dissoc)
        h[j] <- result[1]
        j <- j + 1
    }

    temp <- h*h + K1[i_flag_15]*h + K1[i_flag_15]*K2[i_flag_15]
    HCO3[i_flag_15] <- (DIC[i_flag_15]*K1[i_flag_15]*h) / temp 
    CO3[i_flag_15]  <- (DIC[i_flag_15]*K1[i_flag_15]*K2[i_flag_15]) / temp
    CO2[i_flag_15]  <- h*HCO3[i_flag_15]/K1[i_flag_15]
    fCO2[i_flag_15] <- CO2[i_flag_15]/K0[i_flag_15]
    PH[i_flag_15] <- -log10(h)
    H[i_flag_15] <- h
    
    #Initialize new 'potential' & "insitu' arrays (this are overwritten below with correct values)
    fCO2pot <- fCO2
    pCO2pot <- pCO2
    fCO2insitu <- fCO2
    pCO2insitu <- pCO2

    # ------------ Compute pCO2 for cases 1 to 15 (pCO2, pCO2insitu, & pCO2pot are NOT input variables) 
    # Indices of flag elements where 1 <= flag <= 15
    i_flag <- which (flag >= 1 & flag <= 15)
    # 1) Classic pCO2: compute from in situ T and surface P (Patm)
    tk <- TK[i_flag]
    B  <- -1636.75 + 12.0408*tk  - 0.0327957*(tk*tk)   + 0.0000316528*(tk*tk*tk);
    fugfac <- exp((Patm[i_flag])*(B + 2*((1-fCO2[i_flag])^2)*(57.7-0.118*tk))/(82.05736*tk))
    pCO2[i_flag]    <- fCO2[i_flag] / fugfac
    # 2) In-situ pCO2: compute from in situ T and total P (Patm + Phydrostatic)
    fCO2insitu[i_flag] <- fCO2[i_flag] * K0[i_flag] / K0insitu[i_flag]
    xCO2approx <- pCO2[i_flag]       #Results virtually insensitive to this (do not use fCO2insitu)
    fugfacinsitu <- exp((Patm[i_flag] + P[i_flag]/1.01325)*(B + 2*((1-xCO2approx)^2)*(57.7-0.118*tk))/(82.05736*tk))
    pCO2insitu[i_flag]    <- fCO2insitu[i_flag] / fugfacinsitu
    # 3) Potential pCO2: compute from potential T & surface P (Patm only)  
    tkp <- theta(S=SP[i_flag], T=InsT[i_flag], P=P[i_flag], Pref=0) + 273.15       #Potential temperature in Kelvin 
    Bpot   <- -1636.75 + 12.0408*tkp - 0.0327957*(tkp*tkp) + 0.0000316528*(tkp*tkp*tkp);
    fCO2pot[i_flag] <- fCO2[i_flag] * K0[i_flag] / K0pot[i_flag]
    fugfacpot <- exp((Patm[i_flag])*(Bpot + 2*((1-fCO2pot[i_flag])^2)*(57.7-0.118*tkp))/(82.05736*tkp))
    pCO2pot[i_flag] <- fCO2pot[i_flag] / fugfacpot

    # ------------ Compute fCO2 for cases 21 to 25 (pCO2, pCO2pot or pCO2insitu is an input variable)
    # ------------ Notice that there are three cases : insitu, potential or standard
    
    # Indices of flag elements where 21 <= flag <= 25
    i_flag <- which (flag >= 21 & flag <= 25)
    # gas option for those elements
    gas_option <- gas[i_flag]

    # ------------ Process "insitu" case

    # select in "i_flag" those indices that have "insitu" option
    i_flag_insitu <- i_flag[gas_option == "insitu"]

    #if there is any such case
    if (any(i_flag_insitu))
    {
        tk    <- TK[i_flag_insitu]
        tkp   <- theta(S=SP[i_flag_insitu], T=InsT[i_flag_insitu], P=P[i_flag_insitu], Pref=0) + 273.15       #Potential temperature in Kelvin
        B     <- -1636.75 + 12.0408*tk  - 0.0327957*(tk*tk)   + 0.0000316528*(tk*tk*tk);
        Bpot  <- -1636.75 + 12.0408*tkp - 0.0327957*(tkp*tkp) + 0.0000316528*(tkp*tkp*tkp);

        # From in situ pCO2 (in situ T, atm + hydro P), compute potential and standard fCO2 and pCO2
        # a) define input as pCO2insitu & compute fCO2insitu
        pCO2insitu[i_flag_insitu] <- var1[i_flag_insitu] * 1e-6
        # xCO2approx <- pCO2insitu[i_flag_insitu]  (this would be a very gross overestimate at say 5000 m)
        xCO2approx <- 0.    #This poor approximation has virtually no effect on computed result
        fugfacinsitu <- exp((Patm[i_flag_insitu]+P[i_flag_insitu]/1.01325)*(B    + 2*((1-xCO2approx)^2)*(57.7-0.118*tk)) /(82.05736*tk))
        fCO2insitu[i_flag_insitu] <- pCO2insitu[i_flag_insitu] * fugfacinsitu
        # b) compute fCO2pot & fCO2
        fCO2pot[i_flag_insitu]    <- fCO2insitu[i_flag_insitu] * K0insitu[i_flag_insitu] / K0pot[i_flag_insitu]
        fCO2[i_flag_insitu]       <- fCO2insitu[i_flag_insitu] * K0insitu[i_flag_insitu] / K0[i_flag_insitu]
        # c) compute pCO2pot & pCO2
        fugfac       <- exp((Patm[i_flag_insitu])                  *(B    + 2*((1-fCO2[i_flag_insitu])^2)      *(57.7-0.118*tk)) /(82.05736*tk))
        fugfacpot    <- exp((Patm[i_flag_insitu])                  *(Bpot + 2*((1-fCO2pot[i_flag_insitu])^2)      *(57.7-0.118*tkp))/(82.05736*tkp))
        pCO2pot[i_flag_insitu] <- fCO2pot[i_flag_insitu] / fugfacpot
        pCO2[i_flag_insitu]    <- fCO2[i_flag_insitu]    / fugfac
    }

    # ------------ Process "potential" case

    # select in "i_flag" those indices that have "potential" option
    i_flag_potential <- i_flag[gas_option == "potential"]

    #if there is any such case
    if (any(i_flag_potential))
    {
        tk    <- TK[i_flag_potential]
        tkp   <- theta(S=SP[i_flag_potential], T=InsT[i_flag_potential], P=P[i_flag_potential], Pref=0) + 273.15       #Potential temperature in Kelvin
        B     <- -1636.75 + 12.0408*tk  - 0.0327957*(tk*tk)   + 0.0000316528*(tk*tk*tk);
        Bpot  <- -1636.75 + 12.0408*tkp - 0.0327957*(tkp*tkp) + 0.0000316528*(tkp*tkp*tkp);

        # From potential pCO2 (potential T, atm P), compute standard and in situ fCO2 and pCO2
        # a) define input as pCO2pot & compute fCO2pot
        pCO2pot[i_flag_potential] <- var1[i_flag_potential] * 1e-6
        fugfacpot <- exp((Patm[i_flag_potential]                     )*(Bpot + 2*((1-pCO2pot[i_flag_potential])^2)*(57.7-0.118*tkp))/(82.05736*tkp))
        fCO2pot[i_flag_potential] <- pCO2pot[i_flag_potential] * fugfacpot
        # b) compute fCO2 & fCO2insitu
        fCO2[i_flag_potential]       <- fCO2pot[i_flag_potential] * K0pot[i_flag_potential] / K0[i_flag_potential]
        fCO2insitu[i_flag_potential] <- fCO2pot[i_flag_potential] * K0pot[i_flag_potential] / K0insitu[i_flag_potential]
        # c) Compute pCO2 & pCO2insitu
        fugfac       <- exp((Patm[i_flag_potential]                  )*(B    + 2*((1-fCO2[i_flag_potential])^2)*(57.7-0.118*tk))/(82.05736*tk))
        fugfacinsitu <- exp((Patm[i_flag_potential]+P[i_flag_potential]/1.01325)*(B    + 2*((1-fCO2[i_flag_potential])^2)*(57.7-0.118*tk))/(82.05736*tk))
        pCO2[i_flag_potential]      <- fCO2[i_flag_potential]      / fugfac
        pCO2insitu[i_flag_potential] <- fCO2insitu[i_flag_potential] / fugfacinsitu
    }

    # ------------ Process "standard" case

    # select in "i_flag" those indices that have "standard" option
    i_flag_standard <- i_flag[gas_option == "standard"]

    #if there is any such case
    if (any(i_flag_standard))
    {
        tk    <- TK[i_flag_standard]
        tkp   <- theta(S=SP[i_flag_standard], T=InsT[i_flag_standard], P=P[i_flag_standard], Pref=0) + 273.15       #Potential temperature in Kelvin
        B     <- -1636.75 + 12.0408*tk  - 0.0327957*(tk*tk)   + 0.0000316528*(tk*tk*tk);
        Bpot  <- -1636.75 + 12.0408*tkp - 0.0327957*(tkp*tkp) + 0.0000316528*(tkp*tkp*tkp);

        # From standard pCO2 (in situ T, atm P), compute potential and in situ fCO2 and pCO2
        # a) define input as pCO2 & compute fCO2
        pCO2[i_flag_standard] <- var1[i_flag_standard] * 1e-6
        fugfac <- exp((Patm[i_flag_standard]                        )*(B    + 2*((1-pCO2[i_flag_standard])^2)*(57.7-0.118*tk))/(82.05736*tk))
        fCO2[i_flag_standard] <- pCO2[i_flag_standard] * fugfac
        # b) compute fCO2pot & fCO2insitu
        fCO2pot[i_flag_standard]    <- fCO2[i_flag_standard] * K0[i_flag_standard] / K0pot[i_flag_standard]
        fCO2insitu[i_flag_standard] <- fCO2[i_flag_standard] * K0[i_flag_standard] / K0insitu[i_flag_standard]
        # c) Compute pCO2pot & pCO2insitu
        fugfacpot    <- exp((Patm[i_flag_standard]                  )*(Bpot + 2*((1-pCO2[i_flag_standard])^2)*(57.7-0.118*tkp))/(82.05736*tkp))
        fugfacinsitu <- exp((Patm[i_flag_standard]+P[i_flag_standard]/1.01325)*(B    + 2*((1-pCO2[i_flag_standard])^2)*(57.7-0.118*tk))/(82.05736*tk))
        pCO2pot[i_flag_standard]   <- fCO2pot[i_flag_standard]   / fugfacpot
        pCO2insitu[i_flag_standard] <- fCO2insitu[i_flag_standard] / fugfacinsitu
    }
    

    # ------------ case 21.) PH and pCO2 given ----
    # Indices of flag elements where flag = 21
    i_flag_21 <- which (flag == 21)
    PH[i_flag_21] <- var2[i_flag_21]
    h <- 10^(-PH[i_flag_21])
    H[i_flag_21] <- h
    CO2[i_flag_21]  <- K0[i_flag_21]*fCO2[i_flag_21]
    HCO3[i_flag_21] <- K1[i_flag_21]*CO2[i_flag_21]/h
    CO3[i_flag_21]  <- K2[i_flag_21]*HCO3[i_flag_21]/h
    DIC[i_flag_21]  <- CO2[i_flag_21] + HCO3[i_flag_21] + CO3[i_flag_21]

    # ------------ case 22.) HCO3 and pCO2 given ----
    # Indices of flag elements where flag = 22
    i_flag_22 <- which (flag == 22)
    HCO3[i_flag_22] <- var2[i_flag_22]
    CO2[i_flag_22]  <- fCO2[i_flag_22]*K0[i_flag_22]
    h <- CO2[i_flag_22]*K1[i_flag_22]/HCO3[i_flag_22]
    CO3[i_flag_22]  <- HCO3[i_flag_22]*K2[i_flag_22]/h
    DIC[i_flag_22]  <- CO2[i_flag_22] + HCO3[i_flag_22] + CO3[i_flag_22]
    PH[i_flag_22] <- -log10(h)
    H[i_flag_22] <- h

    # ------------ case 23.) CO3 and pCO2 given ----
    # Indices of flag elements where flag = 23
    i_flag_23 <- which (flag == 23)
    CO3[i_flag_23] <- var2[i_flag_23]
    h <- sqrt(K0[i_flag_23]*K1[i_flag_23]*K2[i_flag_23]*fCO2[i_flag_23]/CO3[i_flag_23])
    HCO3[i_flag_23] <- h*CO3[i_flag_23]/K2[i_flag_23]
    CO2[i_flag_23]  <- h*HCO3[i_flag_23]/K1[i_flag_23]
    DIC[i_flag_23] <- CO2[i_flag_23] + HCO3[i_flag_23] + CO3[i_flag_23]
    PH[i_flag_23] <- -log10(h)
    H[i_flag_23] <- h

    # ------------ case 24.) ALK and pCO2 given ----
    # Indices of flag elements where flag = 24
    i_flag_24 <- which (flag == 24)
    ALK[i_flag_24] <- var2[i_flag_24]
    CO2[i_flag_24] <- fCO2[i_flag_24]*K0[i_flag_24]

    # From this line on, this case is similar to case 4
    # Calculate [H+] from [CO2] and total alk
    h <- rep(NA, length(i_flag_24))
    j <- 1 
    for(i in (i_flag_24))
    {
        # Create a list containing all the dissociation constants
        dissoc = list(K1_DIC=K1[i], K2_DIC=K2[i], K_BT=Kb[i], 
                        K1_PO4 = K1p[i], K2_PO4 = K2p[i], K3_PO4 = K3p[i], K_Sil = Ksi[i], 
                        K2_Sil = K2si[i], K_NH4 = Kn[i], K_H2S = Khs[i],
                        K_HSO4 = Ks[i], K_HF = Kff[i], K_H2O = Kw[i])

        # use SolveSAPHE v2.0
        result <- solve_pH_from_AT(ALK[i], CO2[i], BOR[i], Pt[i], Sit[i], NH4t[i], HSt[i],
                                    ST[i], FLUO[i], pHscale[i], "CO2", FALSE, dissoc)
        h[j] <- result[1]
        j <- j + 1
    }
    
    HCO3[i_flag_24] <- K1[i_flag_24]*CO2[i_flag_24]/h
    CO3[i_flag_24]  <- K2[i_flag_24]*HCO3[i_flag_24]/h
    PH[i_flag_24]   <- -log10(h)
    H[i_flag_24] <- h
    DIC[i_flag_24]  <- CO2[i_flag_24] + HCO3[i_flag_24] + CO3[i_flag_24]

    # ------------ case 25.) DIC and pCO2 given ----
    # Indices of flag elements where flag = 25
    i_flag_25 <- which (flag == 25)
    DIC[i_flag_25] <- var2[i_flag_25]
    CO2[i_flag_25] <- K0[i_flag_25]*fCO2[i_flag_25]
    # Though case 25 is the same as case 5, computations are made in a different way
    K <- K1[i_flag_25]/K2[i_flag_25]
    b <- K*K0[i_flag_25]*fCO2[i_flag_25]
    c <- (K*K0[i_flag_25]*fCO2[i_flag_25]) * 
         (K0[i_flag_25]*fCO2[i_flag_25]-DIC[i_flag_25])
    D <- b*b - 4*c
    HCO3[i_flag_25] <- (1/2)*(-b + sqrt(D))
    CO3[i_flag_25]  <- DIC[i_flag_25] - CO2[i_flag_25] - HCO3[i_flag_25]
    h <- K1[i_flag_25]*CO2[i_flag_25]/HCO3[i_flag_25]
    PH[i_flag_25] <- -log10(h)
    H[i_flag_25] <- h

    # ------------ CALCULATION OF ALK in cases ----
    cases <- c(1, 2, 3, 5, 6, 7, 9, 10, 12, 14, 21, 22, 23, 25)
    # Indices of flag elements in these cases
    i_flag <- which (flag %in% cases)
    h <- H[i_flag]

    # HCO3[i_flag] <- DIC[i_flag]*h*K1[i_flag]/(h*h + K1[i_flag]*h + K1[i_flag]*K2[i_flag])
    # CO3[i_flag]  <- HCO3[i_flag] * K2[i_flag] / h
    boh4 <- BOR[i_flag]/(1+h/Kb[i_flag])
    oh <- Kw[i_flag]/h 
    temp <- h^3 + K1p[i_flag]*h^2 + K1p[i_flag]*K2p[i_flag]*h + K1p[i_flag]*K2p[i_flag]*K3p[i_flag]
    h3po4 <- Pt[i_flag]*(h^3) / temp
    hpo4  <- Pt[i_flag]*K1p[i_flag]*K2p[i_flag]*h / temp
    po4   <- Pt[i_flag]*K1p[i_flag]*K2p[i_flag]*K3p[i_flag] / temp

    if (fullresult) {
        # adapted to include second dissociation constant of silicate
        siooh3 <- Sit[i_flag]*h*Ksi[i_flag] / 
            (h*h + Ksi[i_flag]*h + Ksi[i_flag]*K2si[i_flag])
        sio2oh2 <- Sit[i_flag]*Ksi[i_flag]*K2si[i_flag] / 
            (h*h + Ksi[i_flag]*h + Ksi[i_flag]*K2si[i_flag])
        nh3 <- NH4t[i_flag]/(1+h/Kn[i_flag])
        hs <- HSt[i_flag]/(1+h/Khs[i_flag])
    }
    else {
        siooh3 <- Sit[i_flag]/(1+h/Ksi[i_flag])
    }
    
    ## calculate Hfree anf Htot
    hfree <- rep(NA, length(i_flag))
    htot  <- rep(NA, length(i_flag))    
    sc <- pHscale[i_flag] 
    st <- ST[i_flag]
    ks  <- Ks[i_flag]
    fluo <- FLUO[i_flag]
    kff  <- Kff[i_flag]
    
    # Where pHscale=="F", pHscale = free scale
    i_sc_F <- which (sc == "F") 
    hfree[i_sc_F] <- h[i_sc_F]  
    htot[i_sc_F]  <- h[i_sc_F] / (1+st[i_sc_F]/ks[i_sc_F]) 
    # Where pHscale=="T", pHscale = total scale
    i_sc_T <- which (sc == "T") 
    hfree[i_sc_T] <- h[i_sc_T] * (1+st[i_sc_T]/ks[i_sc_T])
    htot[i_sc_T]  <- h[i_sc_T]
    # Where pHscale=="SWS", pHscale = SW scale
    i_sc_S <- which (sc == "SWS") 
    hfree[i_sc_S] <- h[i_sc_S] * (1 + st[i_sc_S]/ks[i_sc_S] + fluo[i_sc_S]/kff[i_sc_S])
    htot[i_sc_S]  <- hfree[i_sc_S] / (1+st[i_sc_S]/ks[i_sc_S])

    hso4 <- st/(1+ks/hfree)
    hf   <- fluo/(1+Kf[i_flag]/htot)

    if (fullresult) {
        ALK[i_flag]  <- HCO3[i_flag] + 2*CO3[i_flag] + boh4+oh+hpo4+2*po4+siooh3+2*sio2oh2+nh3+hs-hfree-hso4-hf-h3po4
    }
    else {
        ALK[i_flag]  <- HCO3[i_flag] + 2*CO3[i_flag] + boh4+oh+hpo4+2*po4+siooh3-hfree-hso4-hf-h3po4
    }

    ##########################################################
    # CALCULATION OF ARAGONITE AND CALCITE SATURATION STATE  #
    ##########################################################

    Ca = (0.02128/40.078) * SP/1.80655  #Improved formula from Dickson et al. (2007) as discussed in Orr & Epitalon (2014, revised)
    Oa  <- (Ca*CO3)/Kspa
    Oc  <- (Ca*CO3)/Kspc

    #PCO2 and fCO2 - convert from atm to microatm
    pCO2       <- pCO2*1e6
    fCO2       <- fCO2*1e6    
    pCO2pot    <- pCO2pot*1e6
    fCO2pot    <- fCO2pot*1e6    
    pCO2insitu <- pCO2insitu*1e6
    fCO2insitu <- fCO2insitu*1e6    

    if (fullresult) {
        ####################################################
        # CALCULATION OF FULL ACID-BASE SPECIATION         #
        # not included: HNO3, NO3, HNO2, NO2, S(2-), H2SO4 #
        ####################################################
        NH4     <- (H) / (H + Kn) * NH4t      
        NH3     <- (Kn) / (H + Kn) * NH4t
        BOH3    <- (H) / (H + Kb) * BOR
        BOH4    <- (Kb) / (H + Kb) * BOR
        H3PO4   <- (H^3) / (H^3 + K1p * H^2 + K1p * K2p * H + K1p * K2p * K3p) * Pt
        H2PO4   <- (K1p * H^2) / (H^3 + K1p * H^2 + K1p * K2p * H + K1p * K2p * K3p) * Pt
        HPO4    <- (K1p * K2p * H) / (H^3 + K1p * H^2 + K1p * K2p * H + K1p * K2p * K3p) * Pt
        PO4     <- (K1p * K2p * K3p) / (H^3 + K1p * H^2 + K1p * K2p * H + K1p * K2p * K3p) * Pt
        H2S     <- (H) / (H + Khs) * HSt
        HS      <- (Khs) / (H + Khs) * HSt
        SiOH4   <- (H^2) / (H^2 + Ksi * H + Ksi * K2si) * Sit
        SiOOH3  <- (Ksi * H) / (H^2 + Ksi * H + Ksi * K2si) * Sit
        SiO2OH2 <- (Ksi * K2si) / (H^2 + Ksi * H + Ksi * K2si) * Sit
        HF      <- (H) / (H + Kf) * FLUO
        F       <- (Kf) / (H + Kf) * FLUO
        HSO4    <- (H) / (H + Ks) * ST
        SO4     <- (Ks) / (H + Ks) * ST
        OH      <- Kw/H
        
        RES <- data.frame(flag, S, T, Patm, P, PH, CO2, fCO2, pCO2, fCO2pot, pCO2pot, fCO2insitu, pCO2insitu, HCO3, CO3, DIC, ALK, Oa, Oc,
                        NH4, NH3, BOH3, BOH4, H3PO4, H2PO4, HPO4, PO4, H2S, HS, SiOH4, SiOOH3, SiO2OH2, HF, F, HSO4, SO4, H, OH, NH4t, BOR, Pt, HSt, Sit, FLUO, ST)
        names(RES) <- c("flag", "S", "T", "Patm", "P", "pH", "CO2", "fCO2", "pCO2", "fCO2pot", "pCO2pot", "fCO2insitu", "pCO2insitu", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite",
                        "NH4", "NH3", "BOH3", "BOH4", "H3PO4", "H2PO4", "HPO4", "PO4", "H2S", "HS", "SiOH4", "SiOOH3", "SiO2OH2", "HF", "F", "HSO4", "SO4", "H", "OH", "NH4t", "BOR", "Pt", "HSt", "Sit",
                        "FLUO", "ST")
    }
    else {
        RES <- data.frame(flag, S, T, Patm, P, PH, CO2, fCO2, pCO2, fCO2pot, pCO2pot, fCO2insitu, pCO2insitu, HCO3, CO3, DIC, ALK, Oa, Oc)
        names(RES) <- c("flag", "S", "T", "Patm", "P", "pH", "CO2", "fCO2", "pCO2", "fCO2pot", "pCO2pot", "fCO2insitu", "pCO2insitu", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite")
    }
    return(RES)
}
