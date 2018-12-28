# Copyright (C) 2009 Jean-Pierre Gattuso
# Copyright (C) 2018 Jens Daniel Mueller, update
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"tris" <-
function(S=35, T=25, b=0.04, k="m18", warn="y")
  {
  #-------Harmonize input vector length-----
  
  n <- max(length(S), length(T), length(k), length(b))
  
  if(length(S)!=n){S <- rep(S[1],n)}
  if(length(T)!=n){T <- rep(T[1],n)}
  if(length(b)!=n){b <- rep(b[1],n)}
  if(length(k)!=n){k <- rep(k[1],n)}


  #-------Initialise output vector------- 
  
  tris <- rep(NA, n)
  method <- rep(NA, n)
  
  #-------Change temperature scale from [deg C] to [K] ----
  
  TK = T + 273.15;           # T [C]; TK [K]
  
  #-------Calculate tris buffer pH based on input values-------
  
  #--------------------------------------------------------------
  #-------tris characterization by DelValls and Dickson 1998----
    #
    #       doi: 10.1016/S0967-0637(98)00019-3
    #     
    #       pH of TRIS buffered artifical seawater solutions
    #
    #       pH-scale: total scale
    #        
    #       correct for T range : 0 - 45 °C
    #       correct for S range : 20 - 40
    #       correct for equmolal TRIS/TRISH+ molality b: 0.04 mol/kg-H20
  #--------------------------------------------------------------  
  
  tris_d98 <- ((11911.08-18.2499*S-0.039336*S^2)*1/(TK))-
	  366.27059+0.53993607*S+0.00016329*S^2+
	  ((64.52243-0.084041*S)*log(TK))-
	  (0.11149858*(TK))
  
	#--------------------------------------------------------------
	#-------tris characterization by Mueller et al 2018-----------
	#
	#       doi: 10.3389/fmars.2018.00176    
	#     
	#       pH of TRIS buffered artifical seawater solutions
	#       extend to low salinities and different TRIS/TRISH+ molalities
	#
	#       pH-scale: total scale
	#        
	#       correct for T range : 5 - 45 °C
	#       correct for S range : 5 - 40
	#       correct for equmolal TRIS/TRISH+ molalities b: 0.01 - 0.04 mol/kg-H20 for S in 5-20
	#       correct for equmolal TRIS/TRISH+ molalities b: 0.04 mol/kg-H20 for S in 20-40
	#--------------------------------------------------------------
	
  tris_m18 <-
	  -327.3307 -                 
	  2.400270 * S +    
	  8.124630e-2 * S^2 -   
	  9.635344e-4 * S^3 -   
	  
	  9.103207e-2 * TK -
	  1.963311e-3 * S * TK +
	  6.430229e-5 * S^2 * TK -
	  7.510992e-7 * S^3 * TK +
	  
	  56.92797 * log(TK) +
	  5.235889e-1 * S * log(TK) -
	  1.7602e-2 * S^2 * log(TK) +
	  2.082387e-4 * S^3 * log(TK) +
	  
	  11382.97 * (1/TK) -
	  
	  2.417045 * b +
	  7.645221e-2 * b * S +
	  1.122392e-2 * b * TK -
	  3.248381e-4 * b * S * TK -
	  
	  4.161537 * b^2 +
	  6.143395e-2 * b^2 * S
	
  #-------------------------------------------------------------------
	#-------Choose between tris characterizations (=method)----------
	
  is_d98 <- (k=='d98')
  is_m18 <- (k=='m18')
	
  tris[is_d98] <- tris_d98[is_d98]
  tris[is_m18] <- tris_m18[is_m18]
	
	#-------Assign method -------------------
	
	method[is_d98] <- "DelValls and Dickson (1998)"
  method[is_m18] <- "Mueller et al. (2018)"
	
	#-------Set warnings-----------
	is_w <- warn == "y"
	if (any(is_w & is_d98 & (T>45 | T<0 | S>40 | S<20 | b!=0.04))) 
	  {warning("S, T, and/or b is outside the range of validity for the TRIS buffer pH formulation by DelValls and Dickson (1998).")}

	if (any(is_w & is_m18 & (T>45 | T<5 | S>40 | b>0.04 | b<0.01))) 
	  {warning("S, T, and/or b is outside the range of validity for the TRIS buffer pH formulation by Mueller et al. (2018).")}
	
	#-------Assign attributes and define return value-----------
	
	attr(tris, "method") = method
	attr(tris, "pH scale") = "total scale"
	attr(tris,"unit") <- "mol/kg-soln"
	return(tris)

}

