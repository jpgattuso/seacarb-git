# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye # with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr> 
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr> 
# # This file is part of seacarb. 
# # Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version. 
# # Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 
# # You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA # #  

buffesm <-  function(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", b="l10"){  
  # if the concentrations of total silicate and total phosphate are NA
  # they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  
Carb <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b) 

RES <- data.frame()  
n <- nrow(Carb) 
if(length(b)!=n){ b <- rep(b[1],n)}  #transforme b en vecteur tel que carb le fait
 
k1k2.V <- k1k2 
if(length(k1k2.V)!=n){k1k2.V <- rep(k1k2.V[1], n)} 
kf.V <- kf 
if(length(kf.V)!=n){kf.V <- rep(kf.V[1], n)} 
ks.V <- ks 
if(length(ks.V)!=n){ks.V <- rep(ks.V[1], n)} 
pHscale.V <- pHscale 
if(length(pHscale.V)!=n){pHscale.V <- rep(pHscale.V[1], n)}  	
for (i in (1:nrow(Carb))){  	
S <- Carb[i,2] 	
T <- Carb[i,3] 	
P <- Carb[i,4] 	
PH <- Carb[i,5] 	
h <- 10^(-PH) 	
CO2 <- Carb[i,6] 	
pCO2 <- Carb[i,7] 	
fCO2 <- Carb[i,8] 	
HCO3 <- Carb[i,9] 	
CO3 <- Carb[i,10] 	
DIC <- Carb[i,11] 	
ALK <- Carb[i,12] 	
Oa <- Carb[i,13] 	
Oc <- Carb[i,14]  	
k1k2 <- k1k2.V[i] 	
kf <- kf.V[i] 	
ks <- ks.V[i] 	
pHscale <- pHscale.V[i]  


#-------Constantes----------------  

tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TK = T + tk;           # TK [K]; T[C]  

#---- issues de equic---- 
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille) 
cl3 = Cl^(1/3); 
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength 
iom0 = 19.924*S/(1000-1.005*S); 
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate bor = bor(S=S, b=b[i]);   # (mol/kg), boron total 
fluo = (7*(S/35))*1e-5        # (mol/kg), DOE94 fluoride total  
bor = bor(S=S, b=b[i]);   # (mol/kg), DOE94 boron total

#--------------------------------------------------------------------- 
#--------------------- calcul des K ---------------------------------- #---------------------------------------------------------------------  
K1 <- K1(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2)    
K2 <- K2(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2) 
Ks <- Ks(S=S, T=T, P=P, ks=ks) 
Kf <- Kf(S=S, T=T, P=P, pHscale=pHscale, kf=kf) 
Kw <- Kw(S=S, T=T, P=P, pHscale=pHscale) 
K0 <- K0(S=S, T=T, P=P) 
Kb <- Kb(S=S, T=T, P=P, pHscale=pHscale) 
K1p <- K1p(S=S, T=T, P=P, pHscale=pHscale) 
K2p <- K2p(S=S, T=T, P=P, pHscale=pHscale) 
K3p <- K3p(S=S, T=T, P=P, pHscale=pHscale) 
Ksi <- Ksi(S=S, T=T, P=P, pHscale=pHscale) 
Kspa <- Kspa(S=S, T=T, P=P) 
Kspc <- Kspc(S=S, T=T, P=P) 	 
rho <- rho(S=S,T=T,P=P)  

#--------------------------------------------------------------------- 
#--------------------    buffer effects    --------------------------- #---------------------------------------------------------------------  
#   	Buffer Factors from Egleston et al. (2010), Glob. Biogeochem. Cycles, 24, GB1002, doi:10.1029/2008GB003407  
#       Std definitions needed to comput buffer factors 

Alkc = (2*CO3 + HCO3) 	
Borate  = (bor / (1 + h/Kb))  
# could also use equivalent formulation of Borate = bor * Kb / (Kb + h) 	
oh = Kw / h  
#   	Special definitions needed for buffer-factor calculations 
#       - originally from Table 1 of Egleston et al;  
#       - later modified to comply w/ formulas in Excel sheet of Chris Sabine (23 Aug 2010) 	
#Segle   = (HCO3 + 4*CO3 + (h*Borate/(Kb + h)) + h - oh)  
#Formula from Table 1 (Egleston et al.)         !Incorrect (inverse last sign) 	 
Segle   = (HCO3 + 4*CO3 + (h*Borate/(Kb + h)) + h + oh)  
#Formula from Sabine (Excel sheet, 23 Aug 2010) !Correct 	 
Pegle  = (2*CO2 + HCO3)                                  
#Formula from Table 1 (Egleston et al.)         !Correct  	
#Qegle  = (HCO3 - (h*Borate/(Kb + h)) - h - oh)           
#Formula from Sabine (Excel sheet, 23 Aug 2010) !Correct  	 
Qegle  = 2*Alkc - Segle                                 
 #Formula from derivation by J. Orr (same result as line just above) !Correct 
# #       Compute 6 buffer factors: 
#       *** NOTE - units of buffer factors (mol/kg) 
#                - to convert to mmol/kg (units shown by Egleston), multiply each factor by 1000 after calling this routine 
	
gammaDIC = (DIC - (Alkc*Alkc)/Segle) 	
gammaALK = ( (Alkc*Alkc - DIC*Segle) / Alkc ) 	
betaDIC  = ( (DIC*Segle - Alkc*Alkc) / Alkc ) 	
betaALK  = (Alkc*Alkc/DIC - Segle) 	
#omegaDIC = ( DIC - (Alkc*Pegle/HCO3) )                  
 #INCORRECT - Formula from Table 1 (Egleston et al.) - replace HCO3 by Qegle 	
#omegaALK = (Alkc - DIC*HCO3/Pegle)                       
#INCORRECT - Formula from Table 1 (Egleston et al.) - replace HCO3 by Qegle 	
omegaDIC = ( DIC - (Alkc*Pegle/Qegle) )                  
 #CORRECT   - Formula from Sabine (Excel sheet, 23 Aug 2010) 	
omegaALK = ( Alkc - DIC*Qegle/Pegle )                    
 #CORRECT   - Formula from Sabine (Excel sheet, 23 Aug 2010)  
#	Revelle Factor (agrees exactly with that computed in seacarb's "buffer" function (BetaD, from Frankignoulle, 1994)  
R = (DIC/gammaDIC)  	
col <- c("gammaDIC", "betaDIC", "omegaDIC", "gammaALK", "betaALK", "omegaALK", "R") 	
res <- data.frame(gammaDIC, betaDIC, omegaDIC, gammaALK, betaALK, omegaALK, R) 	
RES <- rbind(RES, res) 	}
return(RES) }   
