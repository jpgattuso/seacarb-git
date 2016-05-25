# Copyright (C) 2015 Mathilde Hagens (M.Hagens@uu.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##############################################################################################
# Function for the calculation of the second dissociation constant for silicate
# Source: modified (corrected) version from Hofmann et al (Aquatic Geochemistry, 2010)
#############################################################################################
Ksi2 <- function (S = 35, T = 25, P = 0, pHscale = "T", kSWS2scale = 0) 
{
  nK <- max(length(S), length(T), length(P), length(pHscale), 
            length(kSWS2scale))
  if (length(S) != nK) {
    S <- rep(S[1], nK)
  }
  if (length(T) != nK) {
    T <- rep(T[1], nK)
  }
  if (length(P) != nK) {
    P <- rep(P[1], nK)
  }
  if (length(pHscale) != nK) {
    pHscale <- rep(pHscale[1], nK)
  }
  tk = 273.15
  TC = T + tk
  lnK <- 8.9613 + (-4465.18) / TC + (-0.021952) * TC
  Ksi2 <- 10^(lnK)
  if (any(P != 0)) 
    Ksi2 <- Pcorrect(Kvalue = Ksi2, Ktype = "Ksi", T = T, S = S, 
                     P = P, pHscale = "SWS", 1, 1)
  is_total <- pHscale == "T"
  is_free <- pHscale == "F"
  if (any(is_total) || any(is_free)) {
    if (missing(kSWS2scale)) {
      kSWS2scale <- rep(1, nK)
      if (any(is_total)) 
        kSWS2scale[is_total] <- kconv(S = S[is_total], 
                                      T = T[is_total], P = P[is_total])$kSWS2total
      if (any(is_free)) 
        kSWS2scale[is_free] <- kconv(S = S[is_free], 
                                     T = T[is_free], P = P[is_free])$kSWS2free
    }
    else if (length(kSWS2scale) != nK) {
      kSWS2scale <- rep(kSWS2scale[1], nK)
    }
    Ksi2 <- Ksi2 * kSWS2scale
  }
  pHsc <- rep(NA, nK)
  pHsc[is_total] <- "total scale"
  pHsc[is_free] <- "free scale"
  pHsc[!is_total & !is_free] <- "seawater scale"
  if (any(T > 45 | S > 45 | T < 0 | S < 0)) {
    warning("S and/or T is outside the range of validity of the formulation available for Ksi2 in seacarb.")
  }
  attr(Ksi2, "unit") = "mol/kg-soln"
  attr(Ksi2, "pH scale") = pHsc
  return(Ksi2)
}
