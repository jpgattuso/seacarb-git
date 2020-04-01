# Copyright (C) 2020 Kimberlee Baldry
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#

sir <- function(flag, var1, var2, S = 35, T = 25, Patm = 1, P = 0, 
               Pt = 0, Sit = 0, k1k2 = "x", kf = "x", ks = "d", pHscale = "T", 
               b = "u74", gas = "potential", warn = "y", eos = "eos80", 
               long = 1e+20, lat = 1e+20){
  Carb_out = carb(flag, var1, var2, S , T , Patm , P , 
                   Pt, Sit , k1k2 , kf , ks , pHscale , 
                   b , gas , warn , eos , 
                   long , lat )
  ### pH conversions
  # carb allows multiple pH scale inputs
  # set up flag vectors
  free_flag <- NA
  total_flag <- NA
  sws_flag <- NA
  for(ph in 1:length(pHscale)){
    if(pHscale[ph] == "F"){
      free_flag[ph] <- NA
      total_flag[ph] <- 2
      sws_flag[ph] <- 6
          }
    if(pHscale[ph] == "T"){
      free_flag[ph] <- 4
      total_flag[ph] <- NA
      sws_flag[ph] <- 3
          }
    if(pHscale[ph] == "SWS"){
      free_flag[ph] <- 5
      total_flag[ph] <- 1
      sws_flag[ph] <- NA
          }}
  pH_f <- Carb_out$pH
  pH_t <- Carb_out$pH
  pH_sws <- Carb_out$pH
  #identify those that need conversion
  idx_f <- is.finite(free_flag)
  idx_t <- is.finite(total_flag)
  idx_sws <- is.finite(sws_flag)
  #convert only those not in scale
  if(length(which(idx_f)) != 0){pH_f[idx_f] = pHconv(flag = free_flag[idx_f], pH=Carb_out$pH[idx_f], S=Carb_out$S[idx_f], T=Carb_out$T[idx_f], P=Carb_out$P[idx_f], ks = ks)}
  if(length(which(idx_t)) != 0){pH_t[idx_t] = pHconv(flag = total_flag[idx_t], pH=Carb_out$pH[idx_t], S=Carb_out$S[idx_t], T=Carb_out$T[idx_t], P=Carb_out$P[idx_t], ks = ks)}
  if(length(which(idx_sws)) != 0){pH_sws[idx_sws] = pHconv(flag = sws_flag[idx_sws], pH=Carb_out$pH[idx_sws], S=Carb_out$S[idx_sws], T=Carb_out$T[idx_sws], P=Carb_out$P[idx_sws], ks = ks)}
  
  ### calculate H ions on different scales
  H_f <- (10^(-pH_f)*1000000) * 1000 / rho(S=Carb_out$S, T=Carb_out$T, P=Carb_out$P) #mol/L to mol/kg (devide by kg/L)
  H_t <- (10^(-pH_t)*1000000) * 1000 / rho(S=Carb_out$S, T=Carb_out$T, P=Carb_out$P) #mol/L to mol/kg (devide by kg/L)
  H_sws <- (10^(-pH_sws)*1000000) * 1000 / rho(S=Carb_out$S, T=Carb_out$T, P=Carb_out$P) #mol/L to mol/kg (devide by kg/L)
  
  ### SIR calculation
  sir <- (Carb_out$HCO3) / H_f
  
  
  ### output
  RES <- data.frame(flag, Carb_out$S, Carb_out$T, Carb_out$Patm, Carb_out$P, Carb_out$pH, Carb_out$CO2, Carb_out$fCO2, Carb_out$pCO2, 
                    Carb_out$fCO2pot, Carb_out$pCO2pot, Carb_out$fCO2insitu, Carb_out$pCO2insitu, Carb_out$HCO3, Carb_out$CO3, 
                    Carb_out$DIC, Carb_out$ALK, Carb_out$OmegaAragonite, Carb_out$OmegaCalcite, sir, H_f,H_t,H_sws,pH_f,pH_sws,pH_t)
  names(RES) <- c("flag", "S", "T", "Patm", "P", "pH", "CO2", 
                  "fCO2", "pCO2", "fCO2pot", "pCO2pot", "fCO2insitu", "pCO2insitu", 
                  "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite","SIR","H_free","H_sws","H_t","pH_free","pH_sws","pH_t")
  return(RES)
}
