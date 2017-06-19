# Copyright (C) 2016 Mathilde Hagens (M.Hagens@uu.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

buffergen <- function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", 
                      gas="potential", NH4t=0, HSt=0) {
  
  n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), 
           length(ks), length(b), length(gas), length(NH4t), length(HSt))
  if(length(flag)!=n){flag <- rep(flag[1],n)}
  if(length(var1)!=n){var1 <- rep(var1[1],n)}
  if(length(var2)!=n){var2 <- rep(var2[1],n)}
  if(length(S)!=n){S <- rep(S[1],n)}
  if(length(T)!=n){T <- rep(T[1],n)}
  if(length(Patm)!=n){Patm <- rep(Patm[1],n)}
  if(length(P)!=n){P <- rep(P[1],n)}
  if(length(Pt)!=n){Pt <- rep(Pt[1],n)}
  if(length(Sit)!=n){Sit <- rep(Sit[1],n)}
  if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
  if(length(kf)!=n){kf <- rep(kf[1],n)}
  if(length(ks)!=n){ks <- rep(ks[1],n)}
  if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
  if(length(b)!=n){b <- rep(b[1],n)}
  #if(length(gas)!=n){gas <- rep(gas[1],n)}
  if(length(NH4t)!=n){NH4t <- rep(NH4t[1],n)}
  if(length(HSt)!=n){HSt <- rep(HSt[1],n)}
  
  # if the concentrations of total silicate, total phosphate, total nitrite,
  # total ammonium, and total hydrogen sulfide are NA, they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  NH4t[is.na(NH4t)] <- 0
  HSt[is.na(HSt)] <- 0
  
  # Calculate acid-base speciation
  Carbfull <- carbfull(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b,
                   gas=gas, NH4t=NH4t, HSt=HSt)

  #---------------------------------------------------------------------
  #--------------------    buffer effects    ---------------------------
  #---------------------------------------------------------------------
  
  # Define a function to calculate the buffer factors
  bufferfunc <- function(Carbfull, species=species){
  with(as.list(Carbfull),{

    # Calculation of acid-base fractions
    n1  <- NH4 / NH4t
    n2  <- NH3 / NH4t
    b1  <- BOH3 / BOR
    b2  <- BOH4 / BOR
    c1  <- CO2 / DIC
    c2  <- HCO3 / DIC
    c3  <- CO3 / DIC
    p1  <- H3PO4 / Pt
    p2  <- H2PO4 / Pt
    p3  <- HPO4 / Pt
    p4  <- PO4 / Pt
    s1  <- H2S / HSt
    s2  <- HS / HSt
    si1 <- SiOH4 / Sit
    si2 <- SiOOH3 / Sit
    si3 <- SiO2OH2 / Sit
    f1  <- HF / FLUO
    f2  <- F / FLUO
    so1 <- HSO4 / ST
    so2 <- SO4 / ST
    
    # Ceation of output vectors
    dALK.dH <- dtotX.dH <- dALK.dX <- dtotX.dX <- dALK.dpH <- dtotX.dpH <- 
      dH.dALK <- dX.dALK <- dX.dtotX <- dpH.dALK <- dH.dtotX <- dpH.dtotX <-
      rep(NA, length.out = length(species))
    
    # Assign names to output vectors
    names(dALK.dH) <- names(dtotX.dH) <- names(dALK.dX) <- 
      names(dtotX.dX) <- names(dALK.dpH) <- names(dtotX.dpH) <-
      names(dH.dALK) <- names(dX.dALK) <- names(dX.dtotX) <-
      names(dpH.dALK) <- names(dH.dtotX) <- names(dpH.dtotX) <- species
    
    # Define function for calculating ALK.X
    ALK.species <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
         X.name == "DIC") {HCO3 + 2*CO3} else
           if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
             NH3} else
               if(X.name == "H2PO4" || X.name == "H3PO4" || 
                  X.name == "HPO4" || X.name == "PO4" || 
                  X.name == "Pt") {-H3PO4 + HPO4 + 2*PO4} else
                    if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {HS} else
                        if(X.name == "SiOH4" || X.name == "SiOOH3" || 
                           X.name == "SiO2OH2" || X.name == "Sit") {
                          SiOOH3 + 2*SiO2OH2} else
                            if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
                              BOH4} else
                                 if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
                                   -HF} else
                                     if(X.name == "SO4" || X.name == "HSO4" || 
                                        X.name == "ST") {-HSO4} else
                                         NULL
    }
    
    # Define function for determining N
    N.function <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "BOH3" || X.name == "NH4" || 
         X.name == "H2PO4" || X.name == "H2S" || X.name == "SiOH4" || X.name == "F" || 
         X.name == "SO4" || X.name == "DIC" || X.name == "BOR" || 
         X.name == "NH4t" || X.name == "Pt" || 
         X.name == "HSt" || X.name == "Sit" || 
         X.name == "FLUO" || X.name == "ST") {N <- 0} else
           
           if(X.name == "HCO3" || X.name == "NH3" || 
              X.name == "HPO4" || X.name == "HS" ||
              X.name == "SiOOH3" || X.name == "BOH4") {N <- 1} else
                
                if(X.name == "CO3" || X.name == "PO4" || X.name == "SiO2OH2") {N <- 2} else
                     
                     if(X.name == "H3PO4" || X.name == "HF" || X.name == "HSO4") {
                       n <- -1} else
                         
                         NULL
    }
    
    # Define function for calculating dALKdH  
    dALKdH.function <- function(X.name)
    {
      N <- N.function(X.name)
      ALK.X <- ALK.species(X.name)
      if(DIC!=0) {
        if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
           X.name == "DIC") {
          dALKdH.CO2 <- (-1/H)*(-N*ALK.X + HCO3 + 4*CO3)}  else {
            dALKdH.CO2 <- (-1/H)*(HCO3*(c1-c3) + 2*CO3*(2*c1+c2))}
      } else dALKdH.CO2 <- 0
      if(NH4t!=0) {
        if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
          dALKdH.NH4 <- (-1/H)*(-N*ALK.X + NH3)} else {
            dALKdH.NH4 <- (-1/H)*(NH3*n1)}
      } else dALKdH.NH4 <- 0
      if(Pt!=0) {
        if(X.name == "H2PO4" || X.name == "H3PO4" || X.name == "HPO4" || 
           X.name == "PO4" || X.name == "Pt") {
          dALKdH.H2PO4 <- (-1/H)*(-N*ALK.X + H3PO4 + HPO4 + 4*PO4)}  else {
            dALKdH.H2PO4 <- (-1/H)*(-H3PO4*(-p2-2*p3-3*p4) + 
                                     HPO4*(2*p1+p2-p4) + 2*PO4*(3*p1+2*p2+p3))}
      } else dALKdH.H2PO4 <- 0
      if(HSt!=0) {
        if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {
          dALKdH.H2S <- (-1/H)*(-N*ALK.X + HS)}  else {
            dALKdH.H2S <- (-1/H)*(HS*s1)} 
      } else dALKdH.H2S <- 0
      if(Sit!=0) {
        if(X.name == "SiOH4" || X.name == "SiOOH3" || 
           X.name == "SiO2OH2" || X.name == "Sit") {
          dALKdH.SiOH4 <- (-1/H)*(-N*ALK.X + SiOOH3 + 4*SiO2OH2)} else {
            dALKdH.SiOH4 <- (-1/H)*(SiOOH3*(si1-si3) + 
                                     2*SiO2OH2*(2*si1+si2))} 
      } else dALKdH.SiOH4 <- 0
      if(BOR!=0) {
        if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
          dALKdH.BOH3 <- (-1/H)*(-N*ALK.X + BOH4)} else {
            dALKdH.BOH3 <-(-1/H)*(BOH4*b1)}
      } else dALKdH.BOH3 <- 0
      if(FLUO!=0) {
        if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
          dALKdH.F <- (-1/H)*(-N*ALK.X + HF)} else {
            dALKdH.F <- (-1/H)*(-HF*f2)}
      } else dALKdH.F <- 0  
      if(ST!=0) {
        if(X.name == "SO4" || X.name == "HSO4" || X.name == "ST") {
          dALKdH.SO4 <- (-1/H)*(-N*ALK.X + HSO4)} else {
            dALKdH.SO4 <- (-1/H)*(-HSO4*so2)}
      } else dALKdH.SO4 <- 0  
      dALKdH.H <- (-1/H)*(OH+H)       # Internal enhancement of buffering
      return(dALKdH.CO2 + dALKdH.NH4 + dALKdH.H2PO4 + dALKdH.H2S + dALKdH.SiOH4 + dALKdH.BOH3 + 
               dALKdH.F + dALKdH.SO4 + dALKdH.H)
    }

    # Define function for determining totX
    totX.function <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
         X.name == "DIC") {DIC} else
           if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
             NH4t} else
               if(X.name == "H2PO4" || X.name == "H3PO4" || 
                  X.name == "HPO4" || X.name == "PO4" || 
                  X.name == "Pt") {Pt} else
                    if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {HSt} else
                      if(X.name == "SiOH4" || X.name == "SiOOH3" || X.name == "SiO2OH2" || 
                         X.name == "Sit") {Sit} else
                           if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
                             BOR} else
                               if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
                                 FLUO} else
                                   if(X.name == "SO4" || X.name == "HSO4" || 
                                      X.name == "ST") {ST} else
                                        NULL
    }
    
    # Define function to calculate beta.H
    beta.H_func <- function()
    {
      if(DIC!=0) {dALKdH.CO2 <- (-1/H)*(HCO3*(c1-c3) + 2*CO3*(2*c1+c2))} else {
        dALKdH.CO2 <- 0}
      if(NH4t!=0) {dALKdH.NH4 <- (-1/H)*(NH3*n1)} else {dALKdH.NH4 <- 0}
      if(Pt!=0) {dALKdH.H2PO4 <- (-1/H)*(-H3PO4*(-p2-2*p3-3*p4) + 
                                        HPO4*(2*p1+p2-p4) + 2*PO4*(3*p1+2*p2+p3))} else {
        dALKdH.H2PO4 <- 0}
      if(HSt!=0) {dALKdH.H2S <- (-1/H)*(HS*s1)} else {dALKdH.H2S <- 0}
      if(Sit!=0) {dALKdH.SiOH4 <- (-1/H)*(SiOOH3*(si1-si3) + 
                                         2*SiO2OH2*(2*si1+si2))} else {dALKdH.SiOH4 <- 0}
      if(BOR!=0) {dALKdH.BOH3 <-(-1/H)*(BOH4*b1)} else {dALKdH.BOH3 <- 0}
      if(FLUO!=0) {dALKdH.F <- (-1/H)*(-HF*f2)} else {dALKdH.F <- 0}  
      if(ST!=0) {dALKdH.SO4 <- (-1/H)*(-HSO4*so2)} else {dALKdH.SO4 <- 0}  
      dALKdH.H <- (-1/H)*(OH+H)       # Internal enhancement of buffering
      return(dALKdH.CO2 + dALKdH.NH4 + dALKdH.H2PO4 + dALKdH.H2S + dALKdH.SiOH4 + dALKdH.BOH3 + 
               dALKdH.F + dALKdH.SO4 + dALKdH.H)
    }
    
    # Define a function to calculate the Revelle factor 
    Revelle_func <- function()
    {
      X.name   <- "CO2"
      X        <- get(X.name)
      N        <- N.function(X.name)
      ALK.X    <- ALK.species(X.name)
      totX     <- totX.function(X.name)
      
      dALK.dH   <- dALKdH.function(X.name)
      dX.dtotX <- X * H * dALK.dH / (ALK.X^2 + (H*dALK.dH - N*ALK.X)*totX)
      RF <- as.numeric(DIC / CO2 * dX.dtotX)
      
      return(RF)
    }
    
    # Loop that calculates output species per variable
    for (i in 1:length(species))
    {
      
      # In case a change in the total concentration is specified, 
      # change X to the reference species for ALK
      X.name <- species[i]
      
      if(X.name == "DIC") X <- get("CO2") else                      
        if(X.name == "NH4t") X <- get("NH4") else                     
          if(X.name == "Pt") X <- get("H2PO4") else
            if(X.name == "HSt") X <- get("H2S") else 
              if(X.name == "Sit") X <- get("SiOH4") else
                if(X.name == "BOR") X <- get("BOH3") else
                  if(X.name == "FLUO") X <- get("F") else
                    if(X.name == "ST") X <- get("SO4") else
                      X <- get(species[i])
                        
      N              <- N.function(X.name)
      ALK.X          <- ALK.species(X.name)
      totX           <- totX.function(X.name)
      
      if(totX!=0) {
      # Calculation of factors in Table 3a
      dALK.dH[i]     <- dALKdH.function(X.name)
      dtotX.dH[i]    <- (N*totX - ALK.X) / H
      dALK.dX[i]     <- ALK.X / X
      dtotX.dX[i]    <- totX / X
      
      # Transformation from dH to dpH
      dH.dpH         <- -log(10)*H
      
      dALK.dpH[i]    <- dH.dpH * dALK.dH[i]
      dtotX.dpH[i]   <- dH.dpH * dtotX.dH[i]
      
      # Calculation of factors in Table 3b
      dH.dALK[i]     <- totX*H  / (ALK.X^2 + (H*dALK.dH[i] - N*ALK.X)*totX)
      dH.dtotX[i]    <- -ALK.X*H / (ALK.X^2 + (H*dALK.dH[i] - N*ALK.X)*totX)
      dX.dALK[i]     <- -X * (N*totX - ALK.X)  / (ALK.X^2 + (H*dALK.dH[i] - 
                                                             N*ALK.X)*totX)
      dX.dtotX[i]    <- X * H * dALK.dH[i] / (ALK.X^2 + (H*dALK.dH[i] - 
                                                         N*ALK.X)*totX)
      
      # Transformation from dH to dpH
      dpH.dH         <- 1/(-log(10)*H)
      
      dpH.dALK[i]     <- dpH.dH * dH.dALK[i]
      dpH.dtotX[i]   <- dpH.dH * dH.dtotX[i]
      }
    }  
    
    # Calculate the Revelle factor
    RF <- Revelle_func()
    
    # Calculate beta.H
    beta.H <- as.numeric(beta.H_func())
      
    # Combine the results in a list
    result <- list(Carbfull=as.matrix(Carbfull), dALK.dH=dALK.dH, dtotX.dH=dtotX.dH, dALK.dX=dALK.dX, 
                   dtotX.dX=dtotX.dX, dALK.dpH=dALK.dpH, dtotX.dpH=dtotX.dpH, dH.dALK=dH.dALK, 
                   dH.dtotX=dH.dtotX, dX.dALK=dX.dALK, dX.dtotX=dX.dtotX, dpH.dALK=dpH.dALK, 
                   dpH.dtotX=dpH.dtotX, beta.H=beta.H, RF=RF)
    
    # Return result
    return(result)
  })
}
  
  # Define all species for which bufferfunc needs to be calculated
  species <- names(Carbfull)[c(7,14:15,20:36,39:40,16,41:45)]
  
  # Create an empty output list
  output <- list()
  
  # Loop buffer function over results of Carbfull
  for(j in 1:dim(Carbfull)[1]){
  
    # Run function
    result <- bufferfunc(Carbfull[j,], species)
  
    # Transfer the result to an output
    for(k in 1:length(result)){
      if(j==1) {
        output[[k]] <- result[[k]]
        names(output)[k] <- names(result)[k]} else {
          output[[k]] <- rbind(output[[k]],result[[k]])}
    }
    }
  
  # Return statement
  return(output)
  
}