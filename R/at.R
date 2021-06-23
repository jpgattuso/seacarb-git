# Copyright (C) 2021 May-Linn Paulsen and Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Updates by May-Linn Paulsen (2021) for a) air-buoyancy correction, b) option to use NaCl-fortified or non-fortified HCl, and c) choice of pH scale. 

# "at" <-
  at <- function(S=35, T=25, C=0.1, HClunit="mol/L",CNaCl=0, HCltemp=20, pHtot="yes", pHTris=NULL, ETris=NULL, weightSample, E, volume){
    
    # S : salinity in practical salinity units (constant)
    # T : temperature in deg. C (vector or constant)
    # C : concentration of acid in (units defined) (constant)
    # d : density of acid in g/cm^3 (constant)
    # pHTris : pH used for the calibration of the electrode  (constant)
    # ETris : voltage used for the calibration of the electrode, in mV (constant)
    # E : voltage recorded during the titration in mV (vector)
    # volume : volume of acid in mL (vector)
    # weightSample : weight of the sample g (constant)
    
    R <- 8.31447215        #constant
    F <- 96485.339924 
    
    # check that T is given as a vector
    test <- length(T) != length(E)
    if(test) { T <- rep(T[1], length(E)) }
    
    Tk <- T + 273.15
    
    if (E[1]<= 1){
      E = E*1000
    }
    
    # creation of a table p
    p <- data.frame(E=E, volume=volume, Tk=Tk)
    z <- p
    # transform mV in pH (total scale)
    if(!is.null(pHTris)&!is.null(ETris)){
      pH <- pHTris + (ETris/1000-E/1000)/(R*Tk*log(10)/F)     
      
      # creation of a table p
      p <- data.frame(p, pH=pH)
      #to use only the value between pH 3.5 and 3 (according to Dickson, 2007):
      iii <- which((3<= p$pH)&(p$pH<=3.5))
      z<- p[iii,]
    }
    
    options(digits=9)

    #convert from weight to mass using sample and HCl density
    dSample <- seacarb::rho(S,T[1],0)/1000
    MMHCl <- 35.45+1.01 # molar mass of HCl
    if(HClunit == 'mol/L' | HClunit == 'mol L-1' | HClunit == 'mol/l' | HClunit == 'mol l-1'){
      C <- C/(1000+C*MMHCl)*1000 }# converting from mol/L to mol/kg-sol

    # calculate titrant density based on NaCl fortification and titrant temperature 
    if(CNaCl == 0){ #based on Handbook of CO2 meas, Dickson et al 2007, chapter 5, paragraph 4.4 
      d = 1.0}
      else {
        MMNaCl<- 22.99 + 35.45
        mHCl  <- C/(1000 - C*MMHCl - CNaCl*MMNaCl)*1000
        mNaCl <- CNaCl/(1000 - C*MMNaCl - CNaCl*MMNaCl)*1000
        gHCl  <- 17.854 + 1.460*sqrt(mHCl) - 0.307*mHCl
        gNaCl <- 16.613 + 1.811*sqrt(mNaCl) + 0.094*mNaCl
        gmix  <- ((mHCl*gHCl) + (mNaCl*gNaCl))/(mHCl + mNaCl)
        mT    <- ((MMHCl*mHCl) + (MMNaCl*mNaCl))/(mHCl + mNaCl)
        rhoW  <- (999.84847 + 6.337563*10**-2*HCltemp - 8.523829*10**-3*HCltemp**2 + 6.943248*10**-5*HCltemp**3 - 3.821216*10**-7*HCltemp**4)/1000
        d     <- (rhoW*(10^3 + mT*(mHCl + mNaCl)))/(10^3 + gmix*(mHCl + mNaCl)*rhoW)
      }
    m <- z$volume*d*((1-(0.0012013/8))/(1-(0.0012013/d)))	# Mass of acid
    m0 <- weightSample*((1-(0.0012013/8))/(1-(0.0012013/dSample)))		# Mass of the sample

    #linear estimation of the total alkalinity (gran function):
    #####################################################
    F1<-(m0+m)*exp((z$E/1000)/(R*(z$Tk)/F))
    f<-lm(m~F1)
    TA<-f$coefficients[1]*C/m0[1]
    E0est <- z$E/1000-(R*z$Tk/F)*log((-m0*TA+m*C)/(m0+m))
    Hprime <- exp((z$E/1000-E0est)/(R*z$Tk/F))

    #non linear estimation:
    
    St <- (0.14/96.062)*(S/1.80655)
    Ks <- Ks(S,T[1], 0)
    Z <- 1+ St/Ks
    Ft <- (0.000067/18.998)*(S/1.80655)
    Kf <- exp(874/z$Tk-9.68+0.111*S**(1/2))
    y <- (m/m0)
    regr <- nls(y ~ ((At + (St/(1 + Ks*Z/(f*Hprime)))+(Ft/ (1+Kf/(f*Hprime)))+(f*Hprime/Z))/(C-f*Hprime/Z)),start=list(At=TA, f=1))

    ALK <- summary(regr)$parameters[1]
    attr(ALK, "unit") <- "mol/kg-soln"
    attr(ALK, "name") <- "Total Alkalinity"

    ## At= Total Alkalinity (mol/kg)
    E0 <- (mean(E0est - (R*z$Tk/F)*log(summary(regr)$parameters[2])))*1000
    attr(E0,"unit") <- "mV"
    ALK_E0 <- list("AT" = ALK, "E0" = E0)
    return(ALK_E0)
  }
