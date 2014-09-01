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
"K0" <-
function(S=35,T=25,P=0){

    nK <- max(length(S), length(T), length(P))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}

    #-------Constantes----------------

    #---- issues de equic----
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TC = T + tk;           # TC [C]; T[K]


    #---------------------------------------------------------------------
    #---------------------- K0 (K Henry) ---------------------------------
    #
    #               CO2(g) <-> CO2(aq.)
    #               K0      = [CO2]/ p CO2
    #
    #   Weiss (1974)   [mol/kg/atm]
    #
    #                             

    tmp = 9345.17 / TC - 60.2409 + 23.3585 * log(TC/100);
    nK0we74 = tmp + S*(0.023517-0.00023656*TC+0.0047036e-4*TC*TC);

    ##------------Warnings

    if (any (T>45 | S>45 | T<(-1)) ) {warning("S and/or T is outside the range of validity of the formulation available for K0 in seacarb.")}


    K0= exp(nK0we74);
    attr(K0,"unit") = "mol/kg"	
    return(K0)
}
