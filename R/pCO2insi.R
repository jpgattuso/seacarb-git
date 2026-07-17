# Copyright (C) 2021  Jean-Pierre Gattuso
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
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
#
pCO2insi <- function(pCO2lab = 400, Tlab = 20, SST = 19){
  #Takahashi (1993) recommended by Pierrot et al. (2009)
  pCO2_insi <- pCO2lab * exp(0.0423 * (SST - Tlab))
return(pCO2_insi)
}
