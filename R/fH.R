# Copyright (C) 2024 James Orr
# This file is part of seacarb.
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Compute total hydrogen ion activity coefficient to convert between SWS and NBS scales

fH <- function(S=35, T=25){
    # --------------------------------------------------------------------------------------------------------------
    # ah = fH * [H+]sws
    # where ah  is the activity of hydrogen ion,
    #       fH  is the total activity coefficient, and
    #       [H+]sws = [H+] + [HSO4-] + [HF], or in other words "the hydrogen ion conccentration on the seawater scale"
    #
    # Note:
    # -----
    # pHnbs = -log10(ah)       #pH on the NBS scale
    # phsws = -log10([H+]sws)  #pH on the SWS scale
    # --------------------------------------------------------------------------------------------------------------

    # The activity coefficient (fH) is used to convert from H+ conccentration on SWS scale 
    # to H+ activity (ah), as used for NBS scale and vice versa. 
    # Here, fH is taken from Takahashi et al (1982, GEOSECS Pacific Expedition, Chap 3, p. 80) who say:
    # "fH is the total activity coeff., which includes contributions from HSO4- and HF [as well as H+].
    #  Culberson & Pytkowicz (28) determined fH as a function of temperature and salinity, and
    #  their results can be approximated by:"

    f = (1.2948 - 0.002036*T + (0.0004607 - 0.000001475*T)*S^2)

    # The approach used to compute fH is old. Its use to convert between the NBS and SWS scales is
    # full of uncertainty.  Newer approaches are more complicated (Pitzer equations) but big uncertainties
    # remain (Marion et al., 2011; Pilson, 2013).
    #
    # Culberson, CH, & Pytkowicz, RM (1973). Ionization of water in seawater. Marine Chemistry, 1(4), 309-316.
    #
    # Marion GM, Millero FJ, Camoes MF, Spitzer P, Feistel R, Chen CTA. 2011. pH of seawater. Marine Chemistry 126: 89-96
    #
    # Pilson MEQ. 2013. An introduction to the chemistry of the sea, 2 edn. Cambridge, UK: Cambridge University Press.
    #
    # Takahashi, T. et al (1982). Carbonate chemistry. GEOSECS Pacific Expedition, Volume 3, Hydrographic Data 1973-1974, 77-83.

return(f)
}
