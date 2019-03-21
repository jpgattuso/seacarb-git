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


"pHspec" <-
function(S=35, T=25, R=1, d="mCP", k="m18", warn="y")
  {
  #-------Harmonize input vector length-----
  
  n <- max(length(S), length(T), length(d), length(k))
  
  if(length(S)!=n){S <- rep(S[1],n)}
  if(length(T)!=n){T <- rep(T[1],n)}
  if(length(R)!=n){R <- rep(T[1],n)}
  if(length(d)!=n){d <- rep(d[1],n)}
  if(length(k)!=n){k <- rep(k[1],n)}


  #-------Initialise output vector------- 
  
  pHspec <- rep(NA, n)
  method <- rep(NA, n)
  
  #-------Change temperature scale from [deg C] to [K] ----
  
  TK = T + 273.15;           # T [C]; TK [K]
  
  #-------Calculate sample pH based on input values-------
  
 	#--------------------------------------------------------------
	#-------mCP characterization by Mueller and Rehder 2018-----------
	#
	#       https://doi.org/10.3389/fmars.2018.00177  
	#     
	#       Calculates pH of seawater, brackish, and freshwater samples
	#       based on absorbance ratios R obtained from spectrophotometric 
	#       measurements with the dye m-Cresol purple 
	#
	#       pH-scale: total scale
	#        
	#       correct for T range : 5 - 35 °C
	#       correct for S range : 0 - 40
	#--------------------------------------------------------------
	
  pHspec_m18 <-
    
    #firsTK seTK of coefficienTKs defines pK2e2 = f(Sal, TKem)
    
    1.08071477e+03                      -
    1.35394946e-01  *S^0.5            -   
    1.98063716e+02  *S^1.5            +
    6.31924397e+01  *S^2              -
    5.18141866e+00  *S^2.5            -
    2.66457425e+04  *TK^-1             +
    5.08796578e+03  *S^1.5 * TK^-1   -
    1.62454827e+03  *S^2 * TK^-1     +
    1.33276788e+02  *S^2.5 * TK^-1   -
    1.89671212e+02  *log(TK)           +
    3.49038762e+01  *S^1.5 * log(TK) -
    1.11336508e+01  *S^2 * log(TK)   +
    9.12761930e-01  *S^2.5 * log(TK) +
    3.27430677e-01  *TK              -
    7.51448528e-04  *S^0.5 * TK      +
    3.94838229e-04  *S * TK          -
    6.00237876e-02  *S^1.5 * TK      +
    1.90997693e-02  *S^2 * TK        -
    1.56396488e-03  *S^2.5 * TK      +
    
    #second seTK of coefficienTKs includes TKhe definiTKion of mCP absorpTKiviTKy raTKios e1 and e3/e3
    #as deTKermined by Liu eTK al. (2011) and defines TKhe log-TKerm calculaTKion 
    
    log10(
      (R -
         (-0.007762 + 4.5174e-5*TK)) /
        (1 - (R *  (- 0.020813 + 2.60262e-4*TK + 1.0436e-4*(S-35))))
    )
  
  
  
	
  #-------------------------------------------------------------------
	#-------Choose between mCP characterizations (=method)----------
	
  is_m18 <- (k=='m18')
	pHspec[is_m18] <- pHspec_m18[is_m18]
	
	#-------Assign method -------------------
	
  method[is_m18] <- "Mueller and Rehder (2018)"
	
	#-------Set warnings-----------
	is_w <- warn == "y"
	if (any(is_w & is_m18 & (T>35 | T<5 | S>40))) 
	  {warning("S, and/or T is outside the range of validity for the mCP characterization by Mueller and Rehder (2018).")}
	
	#-------Assign attributes and define return value-----------
	
	attr(pHspec, "method") = method
	attr(pHspec, "pH scale") = "total scale"
	attr(pHspec,"unit") <- "mol/kg-soln"
	return(pHspec)

}

