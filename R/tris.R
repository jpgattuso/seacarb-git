# Copyright (C) 2009 Jean-Pierre Gattuso, 2018 updated by Jens Daniel Mueller
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

"tris" <-
  
function(S=35, T=25, k="m18", b=0.04)
  
  {
  
  n <- max(length(S), length(T), length(k), length(b))
  
  if(length(S)!=n){S <- rep(S[1],n)}
  if(length(T)!=n){T <- rep(T[1],n)}
  if(length(k)!=n){k <- rep(k[1],n)}
  if(length(b)!=n){b <- rep(b[1],n)}

  
  # Initialise output vector  
  tris <- rep(NA, n)
  
  
  for(i in (1:n)){  
  
	delvalls1998 = ((11911.08-18.2499*S[i]-0.039336*S[i]^2)*1/(T[i]+273.15))-
	  366.27059+0.53993607*S[i]+0.00016329*S[i]^2+
	  ((64.52243-0.084041*S[i])*log(T[i]+273.15))-
	  (0.11149858*(T[i]+273.15))
	
	mueller2018 =
	  -327.3307 -
	    2.400270 * S[i] +
	    8.124630e-2 * S[i]^2 -
	    9.635344e-4 * S[i]^3 -
	    
	    9.103207e-2 * T[i] -
	    1.963311e-3 * S[i] * T[i] +
	    6.430229e-5 * S[i]^2 * T[i] -
	    7.510992e-7 * S[i]^3 * T[i] +
	    
	    56.92797 * log(T[i]) +
	    5.235889e-1 * S[i] * log(T[i]) -
	    1.7602e-2 * S[i]^2 * log(T[i]) +
	    2.082387e-4 * S[i]^3 * log(T[i]) +
	    
	    11382.97 * (1/T[i]) -
	    
	    2.417045 * b[i] +
	    7.645221e-2 * b[i] * S[i] +
	    1.122392e-2 * b[i] * T[i] -
	    3.248381e-4 * b[i] * S[i] * T[i] -
	    
	    4.161537 * b[i]^2 +
	    6.143395e-2 * b[i]^2 * S[i]
	
	

	# flag 1 : seawater to total
	if (k[i]=="d98") {tris[i] <- delvalls1998[i]} 
	if (k[i]=="m18") {tris[i] <- mueller2018[i]} 
	
	
	
	attr(tris,"unit") <- "mol/kg"
	return(tris)
	
  }
	
	
}
