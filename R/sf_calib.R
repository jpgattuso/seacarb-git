# Copyright (C) 2017 Samir Alliouane, Lydia Kapsenberg, Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"sf_calib" <-
function(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calSal=35.6){

nK <- max(length(calEint), length(calEext), length(calpH), length(calT), length(calSal))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)
    
    if(length(calEint)!=nK){calEint <- rep(calEint[1], nK)}
    if(length(calEext)!=nK){calEext <- rep(calEext[1], nK)}
    if(length(calpH)!=nK){calpH <- rep(calpH[1], nK)}
    if(length(calT)!=nK){calT <- rep(calT[1], nK)}
    if(length(calSal)!=nK){calSal <- rep(calSal[1], nK)}
    
    # --------------------- sf_calib ----------------------------
    
    # Gas constant, Faraday constant, 
    R <- 8.3145
    FF <- 96487 
    # Temperature dependence oF standard potentials, Martz et al. 2010
    dE0int <- -0.001101
    dE0ext <- -0.001048
    # See Martz et al. 2010 For greater detail
    tempK <- calT + 273.15 # Convert temp From C to K
    S_T <- (R*tempK)/FF*log(10) # Nernst temp dependence
    
    E0int <- calEint-S_T*calpH # Calc E0int From Nernst & pH @ calibration point
    E0int25 <- E0int+dE0int*(25-calT)
    
    Z <- 19.924*calSal/(1000-1.005*calSal) # Ionic strength, Dickson et al. 2007
    SO4_tot <- (0.14/96.062)*(calSal/1.80655)  # Total conservative sulFate
    cCl <- 0.99889/35.453*calSal/1.80655 # Conservative chloride
    mCl <- cCl*1000/(1000-calSal*35.165/35) # mol/kg-H2O
    K_HSO4 <- exp(-4276.1/tempK+141.328-23.093*log(tempK) +
                    (-13856/tempK+324.57-47.986*log(tempK))*Z^0.5 + 
                    (35474/tempK-771.54+114.723*log(tempK))*Z-2698/tempK*Z^1.5 +
                    1776/tempK*Z^2+log(1-0.001005*calSal)) # BisulFate equilibrium const., Dickson et al. 2007
    pHint_free <- calpH+log10(1+SO4_tot/K_HSO4)
    cHfree <- 10^(-pHint_free) # mol/kg-sw
    pHint_free <- pHint_free+log10((1000-calSal*35.165/35)/1000) # mol/kg-H2O
    mHfree <- 10^(-pHint_free) # mol/kg-H2O
    DHconst <- 0.00000343*calT^2+0.00067524*calT+0.49172143 # Debye-Huckel, Khoo et al. 1977
    log10gamma_HCl <- 2*(-DHconst*sqrt(Z)/(1+1.394*sqrt(Z))+(0.08885-0.000111*calT)*Z)
    aHfree_aCl <- mHfree*mCl*10^(log10gamma_HCl) 
    E0ext <- calEext+S_T*log10(aHfree_aCl)
    E0ext25 <- E0ext+dE0ext*(25-calT)
    
    E0 <- data.frame(E0int25, E0ext25)
    return(E0)

    attr(sf_calib,"unit") = "Volts"
    return(sf_calib)
}
