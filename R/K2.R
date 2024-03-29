# Copyright (C) 2008 Jean-Pierre Gattuso and Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye
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
"K2" <-
function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale="x",ktotal2SWS_P0="x",warn="y")
{

    ##---------- check k1k2 parameter
    k1k2_valid_options <- c("cw", "l", "m02", "m06", "m10", "mp2", "p18", "r", "s20", "sb21", "w14", "x")
    test_valid <- is.element(k1k2, k1k2_valid_options)
    
    if (! all (test_valid))
    {
        # If k1k2 is a vector, find indexes of all invalid values in it
        invalids <- which(test_valid == FALSE)
        # Display the first invalid value
        stop ("Invalid parameter k1k2 = ", k1k2[invalids[1]])
    }
    
    ##-------- Create vectors for all input (if vectorial)

    nK <- max(length(S), length(T), length(P), length(k1k2), length(pHscale), length(kSWS2scale) ,length(ktotal2SWS_P0))

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(k1k2)!=nK){k1k2 <- rep(k1k2[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}
    
    ##----------  Initialisation of vectors
    pHsc <- rep(NA,nK) # pHsc : this vector note the actual pHscale because it can change during processing
    pHlabel <- rep(NA,nK)   
    K2 <- rep(NA, nK)

    ##----------Check the validity of the method regarding the T/S range

    is_x <- k1k2 == 'x'
    is_outrange <- T>35 | T<2 | S<19 | S>43
    k1k2[is_x] <- 'l'  ## luecker by default
    k1k2[is_x & is_outrange] <- "w14"  # Waters et al. (2014) if outrange
    
    #-------Constants----------------

    #---- issues of equic----
    tk <- 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK <- T + tk;           # T [C]; T [K]
    
    # --------------------- K2 lueker et al. 2000 ------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Mehrbach et al. (1973) refit by Lueker et al. (2000).
    #
    #   (Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
    #   Dickson, Sabin and Christian , 2007, Chapter 5, p. 14)
    #
    #   pH-scale: 'total'. mol/kg-soln

    is_l <- k1k2 == "l"
    logK2lue <- -471.78/TK[is_l] - 25.9290 + 3.16967*log(TK[is_l]) + 0.01781*S[is_l] - 0.0001122*S[is_l]*S[is_l]
    K2[is_l]<- 10^logK2lue
    pHsc[is_l] <- "T"

    
    # --------------------- K2 Sulpis et al. 2020----------------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #     Papadimitriou et al (2018)  in Geochimica et Cosmochimica Acta 220 (2018) 55–70
    #
    #   pH-scale: 'total'. mol/kg-soln
    is_p18 <- k1k2 == "p18"
    Sp18 <- S[is_p18]
    pK2 <- -323.52692 + 27.557655*sqrt(Sp18) + 0.154922*Sp18 - 2.48396e-4*Sp18*Sp18 +
           (14763.287 - 1014.819*sqrt(Sp18) - 14.35223*Sp18)/TK[is_p18] +
           (50.385807 - 4.4630415*sqrt(Sp18)) * log(TK[is_p18])
    K2[is_p18]<- 10^(-pK2)
    pHsc[is_p18] <- "T"

    # --------------------- K2 Shockman & Byrne 2020----------------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #     Shockman & Byrne (2020)  in Geochimica et Cosmochimica Acta, 2020
    #
    #   pH-scale: 'total'. mol/kg-soln
    is_sb21 <- k1k2 == "sb21"
    Ssb21 <- S[is_sb21]
    pK2 <- 116.8067 - 3655.02/TK[is_sb21] - 16.45817*log(TK[is_sb21]) + 0.04523*Ssb21 - 0.615*sqrt(Ssb21) -
                0.0002799*Ssb21*Ssb21 + 4.969*Ssb21/TK[is_sb21]
    K2[is_sb21]<- 10^(-pK2)
    pHsc[is_sb21] <- "T"


    # --------------------- K2 Sulpis et al. 2020----------------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #     Sulpis et al (2020)  in Ocean Science, 16, 847–862, 2020
    #
    #   pH-scale: 'total'. mol/kg-soln
    is_s20 <- k1k2 == "s20"
    pK2 <- 4226.23/TK[is_s20] - 59.4636 + 9.60817*log(TK[is_s20]) - 0.01781*S[is_s20] + 0.0001122*S[is_s20]*S[is_s20]
    K2[is_s20]<- 10^(-pK2)
    pHsc[is_s20] <- "T"


    # --------------------- K2 Roy et al. 1993----------------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
    #   pH-scale: 'total'. mol/kg-soln
    is_r <- k1k2 == "r"	
    tmp1 = -9.226508 - 3351.6106 / TK[is_r] - 0.2005743 * log(TK[is_r]);
    tmp2 = (-0.106901773 - 23.9722 / TK[is_r]) * sqrt(S[is_r]);
    tmp3 = 0.1130822 * S[is_r] - 0.00846934 * S[is_r]^1.5 + log(1 - 0.001005 * S[is_r]);
    lnK2roy = tmp1 + tmp2 + tmp3;
    K2[is_r] <- exp(lnK2roy)
    pHsc[is_r] <- "T"


    # --------------------- K2 Cai & Wang, 1998----------------------------------------
    #
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Cai and Wang (1998)
    #   pH-scale: 'total'. mol/kg-soln
    is_cw <- k1k2 == "cw"	
    F2  <- -129.24/TK[is_cw] + 1.4381
    pK2[is_cw] <- 2902.39/TK[is_cw] + 0.02379*TK[is_cw] - 6.4980 - 0.3191*F2*S[is_cw]^0.5 + 0.0198*S[is_cw]
    K2[is_cw] <- 10.0^(-pK2)         # this is on the NBS scale
    # Convert from NBS to SWS scale using combined activity coefficient fH 
    # Takahashi (1982, GEOSECS Pacific Expedition, Chap 3, p. 80), who say:
    # "fH is the total activity coeff., which includes contributions from HSO4- and HF [as well as H+].
    #  Culberson and Pytkowicz (28) determined fH as a function of temperature and salinity, and
    #  their results can be approximated by:"
    fH = (1.2948 - 0.002036*TK[is_cw] + (0.0004607 - 0.000001475*TK[is_cw])*S[is_cw]^2)
    K2[is_cw] <- K2[is_cw]/fH[is_cw]  # Convert from NBS to SWS scale
    pHsc[is_cw] <- "SWS"
    # CO2SYS-MATLAB and PyCO2SYS mention that this NBS -> SWS conversion is uncertain at low S due to junction potential
    # Culberson, C. H., & Pytkowicz, R. M. (1973). Ionization of water in seawater. Marine Chemistry, 1(4), 309-316.


    # --------------------- K2 ---------------------------------------
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Mojica Prieto and Millero (2002)
    #   pH-scale: 'seawater'. mol/kg-soln
    is_mp2 <- k1k2 == "mp2"
    Smp2 <- S[is_mp2] 
    pK2 <- -452.0940 + 13.142162 * Smp2 - 8.101e-4 * Smp2*Smp2 + 21263.61/TK[is_mp2] + 68.483143 * log(TK[is_mp2]) +
          (-581.4428 * Smp2 + 0.259601 * Smp2*Smp2) / TK[is_mp2] - 1.967035 * Smp2 * log(TK[is_mp2])
    K2[is_mp2]<- 10^(-pK2)
    pHsc[is_mp2] <- "SWS"

     
    # --------------------- K2 ---------------------------------------
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    # Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
    # Calculated from overdetermined WOCE-era field measurements
    # This is from page 1715
    #
    #   pH-scale: 'seawater'. mol/kg-soln
    is_m02 <- k1k2 == "m02"
    pK2 =  9.867 - 0.01314 * S[is_m02] - 0.01904 * T[is_m02] + 2.448e-5 * T[is_m02]*T[is_m02]
    K2[is_m02]<- 10^(-pK2)
    pHsc[is_m02] <- "SWS"
    
        
    # --------------------- K2 ---------------------------------------
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Millero et al. 2006 Marine Chemistry
    #   pH-scale: 'SWS scale'. mol/kg-soln
    
    is_m06 <- k1k2 == "m06"
    pK2o <- -90.18333 + 5143.692/TK[is_m06] + 14.613358*log(TK[is_m06])
    Sm06 <- S[is_m06] 
    A2 <- 21.0894*Sm06^(0.5) + 0.1248*Sm06 - (3.687e-4)*Sm06^2
    B2 <- -772.483*Sm06^(0.5) - 20.051*Sm06
    C2 <- -3.3336*Sm06^(0.5)
    pK2 <- pK2o + A2 + B2/TK[is_m06] + C2*log(TK[is_m06])
    K2[is_m06] <- 10^(-pK2)
    pHsc[is_m06] <- "SWS"
    
    # --------------------- K2 ---------------------------------------
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Millero 2010 Marine and Fresh water research
    
    is_m10 <- k1k2=="m10"
    
    #   pH-scale: 'SWS scale'. mol/kg-soln
    is_m10_SWS <- is_m10 & (pHscale=="SWS" | P>0)
    A2 <- 21.3728*S[is_m10_SWS]^(0.5) + 0.1218*S[is_m10_SWS] - (3.688e-4)*S[is_m10_SWS]^2
    B2 <- -788.289*S[is_m10_SWS]^(0.5) - 19.189*S[is_m10_SWS]
    C2 <- -3.374*S[is_m10_SWS]^(0.5)
    pK2o <- 5143.692/TK[is_m10_SWS] + 14.613358*log(TK[is_m10_SWS]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_m10_SWS] + C2*log(TK[is_m10_SWS])
    K2[is_m10_SWS] <- 10^(-pK2)
    pHsc[is_m10_SWS] <- "SWS"
    
    #   pH-scale: 'Total scale'. mol/kg-soln
    is_m10_T <- is_m10 & (pHscale=="T" & P==0)
    A2 <- 21.5724*S[is_m10_T]^(0.5) + 0.1212*S[is_m10_T] - (3.714e-4)*S[is_m10_T]^2
    B2 <- -798.292*S[is_m10_T]^(0.5) - 18.951*S[is_m10_T]
    C2 <- -3.403*S[is_m10_T]^(0.5)
    pK2o <- 5143.692/TK[is_m10_T] + 14.613358*log(TK[is_m10_T]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_m10_T] + C2*log(TK[is_m10_T])
    K2[is_m10_T] <- 10^(-pK2)
    pHsc[is_m10_T] <- "T"
    
    #   pH-scale: 'Free scale'. mol/kg-soln
    is_m10_F <- is_m10 & (pHscale=="F" & P==0)
    A2 <- 11.0637*S[is_m10_F]^(0.5) + 0.1379*S[is_m10_F] - (3.688e-4)*S[is_m10_F]^2
    B2 <- -366.178*S[is_m10_F]^(0.5) - 23.288*S[is_m10_F]
    C2 <- -1.810*S[is_m10_F]^(0.5)
    pK2o <- 5143.692/TK[is_m10_F] + 14.613358*log(TK[is_m10_F]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_m10_F] + C2*log(TK[is_m10_F])
    K2[is_m10_F] <- 10^(-pK2)
    pHsc[is_m10_F] <- "F"

    # --------------------- K2 ---------------------------------------
    #   second acidity constant:
    #   [H^+] [CO_3^--] / [HCO_3^-] = K_2
    #
    #   Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014)
    
    is_w14 <- k1k2=="w14"
    
    #   pH-scale: 'SWS scale'. mol/kg-soln
    is_w14_SWS <- is_w14 & (pHscale=="SWS" | P>0)
    A2 <- 21.225890*S[is_w14_SWS]^(0.5) + 0.12450870*S[is_w14_SWS] - (3.7243e-4)*S[is_w14_SWS]^2
    B2 <- -779.3444*S[is_w14_SWS]^(0.5) - 19.91739*S[is_w14_SWS]
    C2 <- -3.3534679*S[is_w14_SWS]^(0.5)
    pK2o <- 5143.692/TK[is_w14_SWS] + 14.613358*log(TK[is_w14_SWS]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_w14_SWS] + C2*log(TK[is_w14_SWS])
    K2[is_w14_SWS] <- 10^(-pK2)
    pHsc[is_w14_SWS] <- "SWS"
    
    #   pH-scale: 'Total scale'. mol/kg-soln
    is_w14_T <- is_w14 & (pHscale=="T" & P==0)
    A2 <- 21.389248*S[is_w14_T]^(0.5) + 0.12452358*S[is_w14_T] - (3.7447e-4)*S[is_w14_T]^2
    B2 <- -787.3736*S[is_w14_T]^(0.5) - 19.84233*S[is_w14_T]
    C2 <- -3.3773006*S[is_w14_T]^(0.5)
    pK2o <- 5143.692/TK[is_w14_T] + 14.613358*log(TK[is_w14_T]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_w14_T] + C2*log(TK[is_w14_T])
    K2[is_w14_T] <- 10^(-pK2)
    pHsc[is_w14_T] <- "T"
    
    #   pH-scale: 'Free scale'. mol/kg-soln
    is_w14_F <- is_w14 & (pHscale=="F" & P==0)
    A2 <- 13.396949*S[is_w14_F]^(0.5) + 0.12193009*S[is_w14_F] - (3.8362e-4)*S[is_w14_F]^2
    B2 <- -472.8633*S[is_w14_F]^(0.5) - 19.03634*S[is_w14_F]
    C2 <- -2.1563270*S[is_w14_F]^(0.5)
    pK2o <- 5143.692/TK[is_w14_F] + 14.613358*log(TK[is_w14_F]) -90.18333
    pK2 <- pK2o + A2 + B2/TK[is_w14_F] + C2*log(TK[is_w14_F])
    K2[is_w14_F] <- 10^(-pK2)
    pHsc[is_w14_F] <- "F"

    ##----------------- Conversion from total to SWS scale
    ##                  if pressure correction needed
    ##                  or pH scale conversion required anyway
    ##                 (in which case SWS may be an intermediate stage of conversion)
    convert <- (pHsc == "T") & ((P > 0) | (pHscale != pHsc))
    if (any (convert))
    {
        ##------------- Convert from total to SWS scale
        # if correction factor (from Total scale to seawater at P=0) not given
        if (missing(ktotal2SWS_P0) || ktotal2SWS_P0[1] == "x")
        {
            # Compute it
            ktotal2SWS_P0  <- rep(1.0,nK)
            ktotal2SWS_P0 <- kconv(S=S[convert], T=T[convert], P=0)$ktotal2SWS
        }
        else
        {
            # Check its length
            if(length(ktotal2SWS_P0)!=nK) ktotal2SWS_P0 <- rep(ktotal2SWS_P0[1], nK)
            # Filter
            ktotal2SWS_P0 <- ktotal2SWS_P0[convert]
        }
        K2[convert] <- K2[convert] * ktotal2SWS_P0
        pHsc[convert] <- "SWS"

        # Just computed constant is on Free scale only if required scale is free scale
        # and if no pressure correction needed  
        # --> No need to convert from free to other scale
        # --> No need to determine conversion factor from free to SWS scale
    }

    # ------------------- Pressure effect --------------------------------
    i_press <- which (P > 0)
    if (length(i_press) > 0)
    {
        # Call Pcorrect() on SWS scale
        K2[i_press] <- Pcorrect(Kvalue=K2[i_press], Ktype="K2", T=T[i_press], 
            S=S[i_press], P=P[i_press], pHscale=pHsc[i_press], 1., 1., warn=warn)
    }

    ## ------------------- Last conversion in the require pHscale ------------------

    convert <- (pHscale != pHsc)
    # At this stage, if any conversion is needed, it is from SWS scale
    # Which pH scales are required ?
    is_total <- convert & (pHscale=="T")
    is_free  <- convert & (pHscale=="F")

    # if any pH scale correction required (from SWS scale)
    if (any(convert))
    {
        # if pH scale correction factor not given
        if (missing(kSWS2scale) || kSWS2scale[1] == "x")
        {
            # Compute it
            kSWS2scale <- rep(1.0,nK)
            if (any(is_total)){ kSWS2scale[is_total] <- kconv(S=S[is_total], T=T[is_total], P=P[is_total])$kSWS2total }
            if (any(is_free)){  kSWS2scale[is_free]  <- kconv(S=S[is_free],  T=T[is_free],  P=P[is_free])$kSWS2free }
        }
        else
            # Check its length
            if (length(kSWS2scale)!=nK) kSWS2scale <- rep(kSWS2scale[1], nK)
        # Apply pH scale correction
        K2[convert] <- K2[convert] * kSWS2scale[convert]
    }

    # Return full name of pH scale
    is_total <- (pHscale=="T")
    is_free  <- (pHscale=="F")
    pHlabel[is_total] <- "total scale"
    pHlabel[is_free]  <- "free scale"
    pHlabel[!is_total & !is_free] <- "seawater scale"

    ##------------Warnings

    is_w <- warn == "y"

    if (   any (is_w & is_l & (T>35 | T<2 | S<19 | S>43))  
        || any (is_w & is_r & (T<0 | T>45 | S<5 | S>45)) 
        || any (is_w & is_mp2 & (T<0 | T>45 | S<5 | S>42)) 
        || any (is_w & is_m02 & (T<(-1.6) | T>35 | S<34 | S>37))
        || any (is_w & is_m06 & (T<1 | T>50 | S<0.1 | S>50)) 
        || any (is_w & (is_m10 | is_w14) & (T<0 | T>50 | S<1 | S>50))
        || any (is_w & is_cw & (T<0.2 | T>30 | S<0.0 | S>40)) 
        || any (is_w & is_p18 & (T<(-6) | T>25 | S<33 | S>100)) 
        || any (is_w & is_s20 & (T<(-1.7) | T>31.8 | S<30.7 | S>37.6)) 
        || any (is_w & is_sb21 & (T<15 | T>35 | S<19.6 | S>41)) )
        warning("S and/or T is outside the range of validity of the formulation chosen for K2.")

    if (any (is_w & (T>50 | S>50))) 
        warning("S and/or T is outside the range of validity of the formulations available for K2 in seacarb.")

    ##---------------Attributes
    method <- rep(NA, nK)
    method[is_cw]   <- "Cai and Wang (1998)"
    method[is_mp2]  <- "Mojica Prieto et al. (2002)"
    method[is_m02]  <- "Millero et al. (2002)"
    method[is_m06]  <- "Millero et al. (2006)"
    method[is_m10]  <- "Millero (2010)"
    method[is_w14]  <- "Waters et al. (2014)"
    method[is_r]    <- "Roy et al. (1993)"
    method[is_p18]  <- "Papadimitriou et al. (2018)"
    method[is_s20]  <- "Sulpis et al. (2020)"
    method[is_sb21] <- "Shockman & Byrne (2021)"
    method[! (is_cw | is_mp2 | is_m02 | is_m06 | is_m10 | is_w14 | is_r | is_p18 | is_s20 | is_sb21) ] <- "Luecker et al. (2000)"
    
    attr(K2,"unit")     = "mol/kg-soln"
    attr(K2,"pH scale") = pHlabel
    attr(K2,"method") = method
    return(K2)
}
