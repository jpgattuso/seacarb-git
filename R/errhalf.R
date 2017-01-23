# errhalf()
# This subroutine determines the halfway curve (in coordinates of sigma1, sigma2) where the total error from
# the measurements (both members of input pair) is equal to the error from all the constants. Midway on that 
# curve (at theta=45°) is the halfway point, where the error from the first member of the input pair (e1) 
# is equal to that from the second (e2).
#
# Input parameters :
#   - epK            :  standard error (uncertainty) in all 7 dissociation constants (a vector)
#   - others         :  same as input of subroutine carb(): scalar or vectors
#
# All parameters may be scalars or vectors except epK, and gas.
#   * gas must be scalar
#   * epK must be vector of seven values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#     these errors are assumed to be equal for all input data.
#
# Returns a 2-dimensional data-frame, with the folowing columns :
#   - H      total error to [H+] concentration  (mol/kg)
#   - pH     total error to pH
#   - CO2    total error to CO2 concentration (mol/kg)
#   - pCO2   total error to "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (µatm)
#   - fCO2   total error to "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (µatm)
#   - HCO3   total error to HCO3 concentration (mol/kg)
#   - CO3    total error to CO3 concentration (mol/kg)
#   - DIC    total error to DIC concentration (mol/kg)
#   - ALK    total error to ALK, total alkalinity (mol/kg)
#   - OmegaAragonite  total error of Omega aragonite (aragonite saturation state)
#   - OmegaCalcite    total error of Omega calcite   (calcite saturation state)
#
errhalf <- 
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, 
         epK=c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02, 0.01), 
         k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", warn="y")
{

#   0. Setup routine matters
#   -------------------------

    # Replicate input given as individual values to have length equal to the number of input data sets
    # Note: evar1 & evar2 are not replicated because they are independant of the input data and each other
    # ----------------------------------------------------------------------------------------------------
    n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), 
             length(k1k2), length(kf), length(pHscale), length(ks), length(b))
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

    # Constant table :  names of input pair variables sorted by flag number
    varnames  = rbind (
        c("H", "CO2"),            # flag = 1    
        c("CO2", "HCO3"),         # flag = 2    
        c("CO2", "CO3"),          # flag = 3    
        c("CO2", "ALK"),          # flag = 4    
        c("CO2", "DIC"),          # flag = 5    
        c("H", "HCO3"),           # flag = 6    
        c("H", "CO3"),            # flag = 7    
        c("H", "ALK"),            # flag = 8    
        c("H", "DIC"),            # flag = 9    
        c("HCO3", "CO3"),         # flag = 10   
        c("HCO3", "ALK"),         # flag = 11   
        c("HCO3", "DIC"),         # flag = 12   
        c("CO3", "ALK"),          # flag = 13   
        c("CO3", "DIC"),          # flag = 14   
        c("ALK", "DIC"),          # flag = 15 
        c("", ""),
        c("", ""),
        c("", ""),
        c("", ""),
        c("", ""),
        c("pCO2", "H"),           # flag = 21   
        c("pCO2", "HCO3"),        # flag = 22   
        c("pCO2", "CO3"),         # flag = 23   
        c("pCO2", "ALK"),         # flag = 24   
        c("pCO2", "DIC")          # flag = 25   
    )
    var1name <- varnames[flag,1]
    var2name <- varnames[flag,2]

    # if epK=NULL, set all pK errors to zero
    if(is.null(epK)) {epK = c(0, 0, 0, 0, 0, 0, 0, 0)}

    # Default value for epK
    if (missing(epK))
    {
        epK <- c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02, 0.01)
    }
    else
    {
        # Check validity of epK
        if (length(epK) == 1 && epK == 0)
        {
            # this means that the caller does not want to account for errors on dissoc. constants
            epK <- rep(0, 8)
        }
        else if (length(epK) != 8)
            stop ("invalid parameter epK: ", epK)
        else
        {
            # Check sign of given epK
            neg_epK <- epK < 0
            epK[neg_epK] <- -epK[neg_epK] 
        }
    }

#   If concentrations of total silicate and phosphate are NA, set them to zero
    Sit[is.na(Sit)] <- 0
    Pt[is.na(Pt)]   <- 0
  
#  1. Compute base values of computed variables
#  --------------------------------------------
   vars <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, 
                k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn)
   H <- 10**(-1.0*vars$pH)

#  Need conditional here to vary the columns saved with the flag value
## vars <- data.frame(vars[,c('CO2','fCO2','pCO2','HCO3','CO3', 'DIC', 'OmegaAragonite','OmegaCalcite')] )

#  Compute H+ and HCO3-/H+, then add latter to data frame
   HCO3H  <- vars$HCO3 / H
   vars <- data.frame(vars, HCO3H) 

#  2. Compute errors from all constants (Kall):
#  --------------------------------------------
   eKall <- errors  (flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit,  
                         evar1=0, evar2=0, eS=0, eT=0, ePt=0, eSit=0, epK=epK,
                         k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn)
#  eKall[isH]$H <- var1*0.0

#  Make logical variable 'Hpair' that is true when pH is a member of the input pair
   Hpair <- (var1name == 'H' || var2name == 'H') 

#  Add HCO3/H ratio to eKall
#  - if pH is an input variable contribution of H is zero (.e., there is NO error in pH from the constants).
#  - for other pairs (without pH) contribution of H must be considered
#  So, if eKall has no column for H+, add it (all zeros) for calc on following line 
   if (Hpair) { eKall$H <- eKall$HCO3*0.0  }
#  just below, the last term in parentheses is zero when pH is a member of the input pair
   eKallHCO3Hr <- sqrt( (eKall$HCO3/vars$HCO3)^2 + (eKall$H/H)^2 )
   HCO3H       <- vars$HCO3H * eKallHCO3Hr
   eKall       <- data.frame(eKall, HCO3H)
   if (Hpair) eKall$H <- NULL

#  3. Define angle
## ---------------
   theta <- c(0, 1, seq(5,85,5), 89, 90)  # angle in degrees
   thetar <- theta * 2 * pi / 360         # in radians
   nang <- length(theta)

#  4. Compute derivatives w.r.t. var1 and var2
#  --------------------------------------------

#  Derivatives of computed variables w.r.t. 'var2' of the input pair (e.g., Alk with the pH-Alk pair).
   dd2 <- derivnum  ('var2', flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, 
                      k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn)

#  Derivatives of computed variables w.r.t. 'var1' of the input pair (e.g., H+ with the pH-Alk pair).
   dd1 <- derivnum  ('var1', flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit,
                      k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn)

#  If pH is a member of the input pair, then add column for d/dH (0 or 1) to dd1 and dd2 data frames
   if (Hpair) {
       # Then one member of the input pair is pH
       if (var1name == 'H') {
           # var1 is pH
           # derivative in last term in numerator in dd1 equation below is one
           dd1$H <- dd1$HCO3 * 0.0 + 1
           # derivative in last term in numerator in dd2 equation below is zero 
           dd2$H <- dd2$HCO3 * 0.0 
       } else if (var2name == 'H') {
           # var2 is pH
           # derivative in last term in numerator in dd1 equation below is zero 
           dd1$H <- dd1$HCO3 * 0.0 
           # derivative in last term in numerator in dd2 equation below is one
           dd2$H <- dd2$HCO3 * 0.0 + 1
       }
   } 
   # Otherwise H is already a column in dd1 and in dd2

#  Add to dd1 the correspondig derivative of the HCO3/H+ ratio (SIR)
#  From derivative quotient rule: dSIR/dH = (H * dHCO3/dH  - HCO3*dH/dH)/H^2  
   HCO3H <- (H*dd1$HCO3 - vars$HCO3*dd1$H) / H^2
   dd1 <- data.frame(dd1, HCO3H)
   if (Hpair) dd1$H <- NULL

#  Add to dd2 the corresponding derivative of the HCO3/H+ ratio (SIR)
#  From derivative quotient rule: dSIR/dALK = (H * dHCO3/dALK  - HCO3*dH/dALK)/H^2
   HCO3H <- (H*dd2$HCO3 - vars$HCO3*dd2$H) / H^2
   dd2 <- data.frame(dd2, HCO3H)
   if (Hpair) dd2$H <- NULL

#  6. Reshape 'eKall' & 'dd' arrays: duplicate each array's single row to match shape of e1 and e2 arrays
#  ------------------------------------------------------------------------------------------------------
#  Duplicate rows in 'eKall', 'dd1', and 'dd2' to each have the same number of rows as 'theta
   eKall <- eKall[rep(row.names(eKall), nang), ]
   dd1  <- dd1[rep(row.names(dd1), nang), ]
   dd2  <- dd2[rep(row.names(dd2), nang), ]

#  Label rows with value of appropriate input error (this step is not really necessary)
   rownames(eKall) <- theta
   rownames(dd1)   <- theta
   rownames(dd2)   <- theta

#  Replicate columns of theta, to multiply that by the data frames in the next section (part 7.)
#  Function to replicate columns: see https://www.r-bloggers.com/a-quick-way-to-do-row-repeat-and-col-repeat-rep-row-rep-col/
   rep.col<-function(x,n){
      matrix(rep(x,each=n), ncol=n, byrow=TRUE)
   }
#  Apply function:
   ncol = dim(eKall)[2]
   thetar <- rep.col(thetar, ncol)
   thetar <- data.frame(thetar)

#  7. Compute critical curve
#  -------------------------
#  Note: sigma1 and sigma2 come in pairs; they are NOT independent, although they are stored in different data frames.
#        Together they define the critical curve where e1^2 + e2^2 = Sum (eKi)^2
   sigma1 <- eKall * cos(thetar) / abs(dd1)    
   sigma2 <- eKall * sin(thetar) / abs(dd2)    

#  8. Compute sigmay (error in computed variable on halfway curve) from e1, e2, and eKall
#  ---------------------------------------------------------------------------------------
   e1 = abs(dd1) * sigma1
   e2 = abs(dd2) * sigma2
   sigmay = sqrt(e1^2 + e2^2 + eKall^2)

#  10. Compute halfway point (where e1 = e2) on halfway curve
#  -----------------------------------------------------------
   sigma1hp <- sqrt(0.5 * (sigmay^2 - eKall^2) / dd1^2 )
   sigma2hp <- sqrt(0.5 * (sigmay^2 - eKall^2) / dd2^2 )

#  11. When pH is a member of the input pair, convert appropriate sigma in H+ to sigma in pH
#  -----------------------------------------------------------------------------------------
   if (Hpair) {
       # Then one member of the input pair is pH
       if (var1name == 'H') {
           # var1 is pH
           sigma1   <- (1/log(10)) * sigma1 / H   #Convert sigma1   from delta H+ to delta pH
           sigma1hp <- (1/log(10)) * sigma1hp / H #Convert sigma1hp from delta H+ to delta pH
       } else if (var2name == 'H') {
           # var2 is pH
           sigma2   <- (1/log(10)) * sigma2 / H   #Convert sigma2   from delta H+ to delta pH
           sigma2hp <- (1/log(10)) * sigma2hp / H #Convert sigma2hp from delta H+ to delta pH
       }
   } 

#  12. Add theta to output data frames
#  ----------------------------------
#  Critical point for error in var1 (given corresponding input of evar2)
   sigma1 <- data.frame(theta, sigma1)  # in mol/kg
   rownames(sigma1) <- theta

#  Critical point for error in var2 (given corresponding input of evar1)
   sigma2 <- data.frame(theta, sigma2)          
   rownames(sigma2) <- theta

#  Error in computed variable along the citical curve
   sigmay <- data.frame(sigmay)          
#  rownames(sigmay) <- theta

   sigma1hp <- data.frame(sigma1hp)          
   sigma2hp <- data.frame(sigma2hp)          

#  Keep only 1st row (all others are identical)
   sigma1hp <- sigma1hp[1,]
   sigma2hp <- sigma2hp[1,]
   sigmay <- sigmay[1,]

   return (list(sigma1, sigma2, sigmay, sigma1hp, sigma2hp))
}
