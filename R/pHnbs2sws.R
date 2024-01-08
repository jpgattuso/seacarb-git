# Copyright (C) 2024 James Orr
# This file is part of seacarb.
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Convert pHNBS to pHSWS using the total hydrogen ion activity coefficient, which depends on T and S

pHnbs2sws <- function(pHNBS, S=35, T=25){
    # --------------------------------------------------------------------------------------------------------------
    # pHNBS = -log10(ah)       #pH on the NBS scale  (-log10 of the activity, not the concentration)
    # phSWS = -log10([H+]sws)  #pH on the SWS scale  (-log10 of the concentration)
    # 
    # Note:
    # -----
    # ah = fH * [H+]sws
    # where ah  is the activity of hydrogen ion,
    #       fH  is the total activity coefficient, and
    #       [H+]sws = [H+] + [HSO4-] + [HF], or in other words "the hydrogen ion conccentration on the seawater scale"
    # --------------------------------------------------------------------------------------------------------------
    #
    # The NBS-to-SWS conversion is done with the total activity coefficient fH (combined activity coeff for H+, HSO4-, and HF)
    # from Takahashi (1982) based on data from Culberson and Pytkowicz (1973). The approach is old and full of uncertainty.
    # Newer approaches are more complicated (Pitzer equations) and big uncertainties remain (Marion et al., 2011; Pilson, 2013).
    #
    # Culberson, CH, & Pytkowicz, RM (1973). Ionization of water in seawater. Marine Chemistry, 1(4), 309-316.
    #
    # Marion GM, Millero FJ, Camoes MF, Spitzer P, Feistel R, Chen CTA. 2011. pH of seawater. Marine Chemistry 126: 89-96
    #
    # Pilson MEQ. 2013. An introduction to the chemistry of the sea, 2 edn. Cambridge, UK: Cambridge University Press.
    #
    # Takahashi, T. et al (1982). Carbonate chemistry. GEOSECS Pacific Expedition, Volume 3, Hydrographic Data 1973-1974, 77-83.

    ah    <-  10^(-pHNBS)      
    hSWS  <-  ah / fH(S=S, T=T)
    pHSWS <- -log10(hSWS) 

return(pHSWS)
}
