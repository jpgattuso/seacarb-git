# Copyright (C) 2016 Mathilde Hagens (M.Hagens@uu.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

K2si <- function (S=35, T=25, P=0, pHscale="T",
                  kSWS2scale="x", ktotal2SWS_P0="x") 
{
  
  nK <- max(length(S), length(T), length(P), length(pHscale), 
            length(kSWS2scale), length(ktotal2SWS_P0))
  
  ##-------- Creation of vector for all entries
  
  if (length(S) != nK) {S <- rep(S[1], nK)}
  if (length(T) != nK) {T <- rep(T[1], nK)}
  if (length(P) != nK) {P <- rep(P[1], nK)}
  if (length(pHscale) != nK) {pHscale <- rep(pHscale[1], nK)}
  
  #-------Constants----------------
  
  tk = 273.15       # [K] (for conversion [deg C] <-> [K])
  TC = T + tk       # TC [C]; T[K]
  
  #--------------------------------------------------------------
  # Dissociation constant of SiO(OH)3 on total scale
  # Source: Wischmeyer et al. (2003), including corrections by 
  # D. Wolf-Gladrow and M. Hagens
  #--------------------------------------------------------------
  
  logK <- 8.9613 + (-4465.18) / TC + (-0.021952) * TC
  K2si <- 10^(logK)
  
  ## ---- Conversion from Total scale to seawater scale before pressure corrections
  
  # if correction factor (from Total scale to seawater at P=0) not given
  if (missing(ktotal2SWS_P0) || ktotal2SWS_P0 == "x")
  {
    # Compute it
    ktotal2SWS_P0 <- kconv(S=S, T=T, P=0)$ktotal2SWS
  }
  else
    # Check its length
    if(length(ktotal2SWS_P0)!=nK) ktotal2SWS_P0 <- rep(ktotal2SWS_P0[1], nK)
  K2si <- K2si * ktotal2SWS_P0
  
  # If needed, pressure correction. Conversion factors for pressure
  # correction are calculated within Pcorrect
  
  if (any(P != 0)) {
    K2si <- Pcorrect(Kvalue = K2si, Ktype = "K2si", T = T, S = S, 
                     P = P, pHscale = "SWS") }
  
  ###----------------pH scale corrections
  
  # Which pH scales are required ?
  is_total <- pHscale=="T"
  is_free   <- pHscale=="F"
  
  # if any pH scale correction required (from total scale)
  if (any(is_total) || any(is_free))
  {
    # if pH scale correction factor not given
    if (missing(kSWS2scale) || kSWS2scale == "x")
    {
      # Compute it
      kSWS2scale <- rep(1.0,nK)
      if (any(is_total)) 
        {kSWS2scale[is_total] <- kconv(S=S[is_total], T=T[is_total], P=P[is_total])$kSWS2total}
      if (any(is_free))  
        {kSWS2scale[is_free]  <- kconv(S=S[is_free], T=T[is_free], P=P[is_free])$kSWS2free}
    }
    else
      # Check its length
      if(length(kSWS2scale)!=nK){kSWS2scale <- rep(kSWS2scale[1], nK)}
    # Apply pH scale correction
    K2si <- K2si*kSWS2scale
  }
  
  # Return full name of pH scale
  pHsc <- rep(NA,nK)
  pHsc[is_total] <- "total scale"
  pHsc[is_free]  <- "free scale"
  pHsc[!is_total & !is_free] <- "seawater scale"
  
  ##------------Warnings
  ## No warning is given; range of validity for this equation not clear from Wischmeyer et al (2003)
  
  attr(K2si, "unit") = "mol/kg-soln"
  attr(K2si, "pH scale") = pHsc
  return(K2si)
}
