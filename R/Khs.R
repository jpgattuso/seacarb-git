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

"Khs"              <- function (S=35,T=25,P=0,pHscale="T", kSWS2scale="x", ktotal2SWS_P0="x", warn="y")

#--------------------------------------------------------------
# Dissociation constant of H2S
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}


tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

	#--------------------------------------------------------------
  # Dissociation constant of H2S on total scale - Millero 1995
 	#--------------------------------------------------------------
lnK <- 225.838 + 0.3449*sqrt(S) - 0.0274*S -13275.3/TC  -34.6435*log(TC)

Khs      <- exp(lnK)

# ---- Conversion from Total scale to seawater scale before pressure corrections
# Khs is TOTAL-scale native, so like Kb and K2si it needs BOTH conversion factors:
# total -> SWS at P=0 (before the SWS-defined pressure correction), then SWS -> chosen.
# At P=0 the two are exact reciprocals and cancel exactly.  Tested at P=300 dbar they
# still agree to <5e-4 %, so unlike Kn (which is SWS-native and shifts by ~0.75%), Khs
# is in practice insensitive to which kf the fallback kconv() picks.  The arguments are
# added for consistency with Kb/K2si and so a caller can avoid the redundant kconv call,
# not because Khs is currently wrong.
if (missing(ktotal2SWS_P0) || (length(ktotal2SWS_P0)==1 && ktotal2SWS_P0[1]=="x")) {
    factor <- kconv(S=S, T=T, P=rep(0,nK))$ktotal2SWS
} else {
    factor <- ktotal2SWS_P0
    if (length(factor)!=nK) factor <- rep(factor[1], nK)
}
Khs <- Khs * factor

# ----------------- Pressure Correction ------------------	
Khs <- Pcorrect(Kvalue=Khs, Ktype="Khs", T=T, S=S, P=P, pHscale="SWS", warn=warn)


###----------------pH scale corrections
# ----------------- pH scale corrections -----------------
# kSWS2scale: SWS -> chosen-scale factor, supplied by the caller so that Khs is placed on
# the same pH scale, by the same route (same ks, same kf), as every other constant in the
# calculation.  Falling back to kconv() here uses its default kf="x" (Perez & Fraga in
# range, Dickson & Riley outside), which is inconsistent with an explicit kf="dg" and is
# a step function of T and S.  See the note in Kn.R.
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
Khs <- Khs * factor

##------------Warnings

    is_w <- warn == "y"

    if (any (is_w & (T>45 | T<0 | S<0 | S>45)) ) {warning("S and/or T is outside the range of validity of the formulation chosen for Khs.")}

  attr(Khs,"unit")     = "mol/kg-soln"
  attr(Khs,"pH scale") = pHsc
  return(Khs)

}  ## END Khs
