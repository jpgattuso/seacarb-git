# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


"Kn" <- 
function (S=35, T=25, P=0, pHscale="T", kSWS2scale="x", warn="y")

#--------------------------------------------------------------
# Dissociation constant of ammonium 
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1],nK)}

#-------Constantes----------------

tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

#--------------------------------------------------------------
# Dissociation constant of ammonium on seawater scale - Millero 1995
#--------------------------------------------------------------
lnK <- -6285.33/TC+0.0001635*TC-0.25444+(0.46532-123.7184/TC)* sqrt(S)+
      (-0.01992+3.17556/TC)*S
Kn <- exp(lnK)

# ----------------- Pressure Correction ------------------	
Kn <- Pcorrect(Kvalue=Kn, Ktype="Kn", T=T, S=S, P=P, pHscale="SWS", warn=warn)


###----------------pH scale corrections
# ----------------- pH scale corrections -----------------
# Kn is SWS-native.  kSWS2scale is the SWS -> chosen-scale factor.  It is supplied by
# the caller (as calculate_carb() does for Kw, Ksi, K1p, K2p, K3p) so that Kn is put on
# the same pH scale, by the same route, as every other constant in the calculation.
#
# If it is not supplied we fall back to computing it here, which calls kconv() with its
# DEFAULT kf="x".  Kf.R then does
#     is_outrange <- T>33 | T<10 | S<10 | S>40
#     kf[is_x] <- "pf" ; kf[is_x & is_outrange] <- "dg"
# i.e. Perez & Fraga IN range and Dickson & Riley outside.  A caller who asked for
# kf="dg" would then silently get a Kn built on a DIFFERENT Kf than K1, K2, Kw, ... ,
# an inconsistency of up to ~0.9%.  Worse, kf="x" is a STEP function of T and S, so Kn
# on the total scale is discontinuous across T=10, T=33, S=10 and S=40, and any
# derivative taken through it inherits that discontinuity.  Pass kSWS2scale to avoid it.
if (missing(kSWS2scale) || (length(kSWS2scale)==1 && kSWS2scale[1]=="x")) {
    factor <- rep(NA,nK)
    for (i in 1:nK) {
        if (pHscale[i]=="T")   factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total
        if (pHscale[i]=="F")   factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free
        if (pHscale[i]=="SWS") factor[i] <- 1
    }
} else {
    factor <- kSWS2scale
    if (length(factor)!=nK) factor <- rep(factor[1], nK)
}
pHsc <- rep(NA,nK)
pHsc[pHscale=="T"]   <- "total scale"
pHsc[pHscale=="F"]   <- "free scale"
pHsc[pHscale=="SWS"] <- "seawater scale"
Kn <- Kn * factor

##------------Warnings

    is_w <- warn == "y"

    if (any (is_w & (T>45 | T<0 | S<0 | S>45)) ) {warning("S and/or T is outside the range of validity of the formulation chosen for Kn.")}

  attr(Kn,"pH scale") = pHsc
  attr(Kn,"unit")     = "mol/kg-soln"
  return(Kn)

}  ## END Kn
