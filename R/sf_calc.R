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
"sf_calc" <-
function(calEint=0.0865, calEext= -0.93, E0int25 =-0.39, E0ext25=-1.46,calT=16.2, calSal=35.6){

nK <- max(length(calEint), length(calEext), length(E0int25), length(E0ext25), length(calT), length(calSal))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)
    
    if(length(calEint)!=nK){S <- rep(S[1], nK)}
    if(length(calEext)!=nK){T <- rep(T[1], nK)}
    if(length(E0int25)!=nK){P <- rep(P[1], nK)}
    if(length(E0ext25)!=nK){P <- rep(P[1], nK)}
    if(length(calT)!=nK){P <- rep(P[1], nK)}
    if(length(calSal)!=nK){P <- rep(P[1], nK)}
    
    # --------------------- sf_calc ----------------------------
    
    # Gas constant, Faraday constant
    FF <- 96487 
    R <- 8.3145 
    # Temperature dependence oF standard potentials, Martz et al. 2010
    dE0Int <- -0.001101 
    dE0Ext <- -0.001048
    # See Martz et al. 2010 For greater detail
    tempK <- calT+273.15 # Convert temp From C to K
    S_T <- (R*tempK)/FF*log(10) # Nernst temp dependence
    pHint_tot <- (calEint-(E0int25+dE0Int*(calT-25)))/S_T # Calc pHint From Nernst
    Z <- 19.924*calSal/(1000-1.005*calSal) # Ionic strength, Dickson et al. 2007
    SO4_tot <- (0.14/96.062)*(calSal/1.80655) # Total conservative sulFate
    cCl <- 0.99889/35.453*calSal/1.80655 # Conservative chloride
    mCl <- cCl*1000/(1000-calSal*35.165/35) # mol/kg-H2O
    K_HSO4 <- exp(-4276.1/tempK+141.328-23.093*log(tempK) +
                    (-13856/tempK+324.57-47.986*log(tempK))*Z^0.5 +
                    (35474/tempK-771.54+114.723*log(tempK))*Z-2698/tempK*Z^1.5 +
                    1776/tempK*Z^2+log(1-0.001005*calSal)) # BisulFate equilibrium const., Dickson et al. 2007
    pHint_free <- pHint_tot+log10(1+SO4_tot/K_HSO4) # free scale mol/kg-sw
    DHconst <- 0.00000343*calT^2+0.00067524*calT+0.49172143 # Debye-Huckel, Khoo et al. 1977
    log10gamma_HCl <- 2*(-DHconst*sqrt(Z)/(1+1.394*sqrt(Z))+(0.08885-0.000111*calT)*Z)
    pHext_free <- -(((E0ext25+dE0Ext*(calT-25))-calEext)-S_T*(log10(mCl)+log10gamma_HCl))/S_T # mol/kg-H2O
    pHext_free <- pHext_free-log10((1000-calSal*35.165/35)/1000) # mol/kg-sw
    pHext_tot <- pHext_free-log10(1+SO4_tot/K_HSO4)
    
    calc <- data.frame(pHint_tot, pHext_tot)
    return(calc)

    attr(sf_calc,"unit") = "Total scale"
    return(sf_calc)
}
