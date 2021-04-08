# Copyright (C) 2008 Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye and Jean-Pierre Gattuso  
# with a most valuable contribution of Bernard Gentili
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
carbb<-
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", badd=0, warn="y", eos="eos80", long=1.e20, lat=1.e20){
    
    RES <- calculate_carb(flag, var1, var2, S, T, Patm, P, Pt, Sit, NH4t=0, HSt=0, k1k2, kf, ks, pHscale, b, gas, badd, warn, eos, long, lat, fullresult=FALSE)
    return(RES)
}
