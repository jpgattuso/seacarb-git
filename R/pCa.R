# Copyright (C) 2008 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"pCa" <- function(flag,  var1,  var2, Ca, S=35, T=20, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", eos="eos80", long=1.e20, lat=1.e20){
  # if the concentrations of total silicate and total phosphate are NA
  # they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  
# Only two options for eos
if (eos != "teos10" && eos != "eos80")
    stop ("invalid parameter eos: ", eos)

# if use of EOS-10 standard
if (eos == "teos10")
{
    # Must convert temperature and salinity from TEOS-10 to EOS-80
    # convert temperature: from Conservative (CT) to in-situ temperature
    # and salinity from Absolute to Practical (SP)
    STeos <- teos2eos_geo (S, T, P, lon, lat)
    InsT <- STeos$T
    SP <- STeos$SP
}
else
{
    InsT <- T
    SP <- S
}

# pass to carb() original S and T to have them insert in the result table
ci <- carb(flag, var1= var1, var2= var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks="d", pHscale=pHscale, b=b, eos=eos, long=long, lat=lat) # initial carbonate chemistry
	#0.01028*(S/35) # calcium concentration (Dickson et al., 2007)
	#co <- ci
	#if (length(Ca) > 1) {for(i in 1:(length(Ca)-1)) {co <- rbind(co,c)}}
	Oa <- Ca * ci$CO3 / Kspa(S=SP, T=InsT, P=P)
	Oc <- Ca * ci$CO3 / Kspc(S=SP, T=InsT, P=P)
	cf <- ci
	if (length(Ca) > 1) {for(i in 1:(length(Ca)-1)) {cf <- rbind(cf,ci)}}
	co <- as.data.frame(c("pCa-initial", rep("pCa-final", nrow(cf))))
	cf$OmegaAragonite <- Oa
	cf$OmegaCalcite <- Oc
	out <- rbind(ci, cf)
	out <- cbind(co, out)
	names(out)[1] <- "comment"
	return(out[,1:20])
}
