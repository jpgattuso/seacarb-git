# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
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
"rho" <-
function(S=35,T=25,P=0){	

	#------------ Density of pure water
	
	rhow = 999.842594 + 6.793952e-2*T -9.095290e-3*T^2 + 1.001685e-4*T^3 -1.120083e-6*T^4 + 6.536332e-9*T^5;
	
	#------------ Density of seawater at 1 atm, P=0
	
	A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T^2 - 8.2467e-7*T^3 + 5.3875e-9*T^4;
	    
	B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T^2; 
	
	C = 4.8314e-4;   
	
	rho0 = rhow + A*S + B*S^(3/2) + C*S^2;

	#-------------- Secant bulk modulus of pure water 
	#
	# The secant bulk modulus is the average change in pressure 
	# divided by the total change in volume per unit of initial volume.
	
	
	Ksbmw = 19652.21 + 148.4206*T - 2.327105*T^2 + 1.360477e-2*T^3 - 5.155288e-5*T^4;
	
	#-------------- Secant bulk modulus of seawater at 1 atm
	
	Ksbm0 = Ksbmw + S*( 54.6746 - 0.603459*T + 1.09987e-2*T^2 - 6.1670e-5*T^3) + S^(3/2)*( 7.944e-2 + 1.6483e-2*T - 5.3009e-4*T^2);
	
	
	#-------------- Secant bulk modulus of seawater at S,T,P
		
	Ksbm = Ksbm0 + P*( 3.239908 + 1.43713e-3*T + 1.16092e-4*T^2 - 5.77905e-7*T^3) + P*S*( 2.2838e-3 - 1.0981e-5*T - 1.6078e-6*T^2) + P*S^(3/2)*1.91075e-4 + P*P*(8.50935e-5 - 6.12293e-6*T + 5.2787e-8*T^2) + P^2*S*(-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T^2);
	 	
	
	#------------- Density of seawater at S,T,P
	rho = rho0/(1-P/Ksbm);
	attr(rho,"unit") = "(kg/m3)"
    return(rho)
#	print(rho);
	}
