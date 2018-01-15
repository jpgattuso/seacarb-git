# errors()
# This subroutine propagates errors to computed carbonate system variables 
# from errors (uncertainties) in six input variables 
#  - pair of carbonate system variables 
#  - nutrients (silicate and phosphate concentrations)
#  - temperature and salinity
# and from errors in dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#
# It propagates error from input to output variables using one of two methods: 
#    * Gaussian (the method of moments without covariance term)
#      It is a very general technique for estimating the second moment of a variable z
#      (variance or standard deviation) based on a first-order approximation to z.
#      This is the default method.
#
#    * Monte Carlo method
#      It relies on repeated random sampling on input errors to obtain numerical results
#      from which standard error is estimated
#
# Input parameters :
#   - evar1, evar2   :  standard error (uncertainty) in var1 and var2 of input pair of carbonate system variables
#   - eS, eT         :  standard error (uncertainty) in Salinity and Temperature (default value = 0.01)
#   - ePt, eSit      :  standard error (uncertainty) in Phosphorus and Silicon total inorganic concentrations
#   - epK            :  standard error (uncertainty) in all 7 dissociation constants (a vector)
#   - method         :  case insensitive character string : "ga", "mo", or "mc"
#                       default is "ga" (gaussian)
#   - r              :  correlation coefficient (from -1 to 1) between var1 and var2 (only for Method of moments')
#   - runs           :  number of random samples (for Monte Carlo method only); default is 10000
#   - others         :  same as input of subroutine carb(): scalar or vectors
#
# All parameters may be scalars or vectors except epK, method, runs and gas.
#   * method, runs and gas must be scalar
#   * epK must be vector of seven values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#     these errors are assumed to be equal for all input data.
#
# In constrast, evar1, evar2, eS, eT, ePt and eSit, 
#   - if vectors, are errors associated with each data point
#   - if scalars, are one error value associated to all data points
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
# Details :
#   Computation time depends on the chosen method, with the Monte Carlo approach being more computationally intensive
#
#   The Monte Carlo method computation time is proportional to the number of runs.
#   But more runs means more accurate results :
#      runs = 10000 appears to be the minimum needed to obtain results with an accuracy of less than 1%
#      Accuracy is inversely proportional to the number of runs.
#
#   Computation time also depends on the type of input pair of variables.
#     For example, the input pair DIC and Alk (flag=15) requires much more computational time 
#     than does the input pair pH and Alkalinity (flag=8).
#
errors <- 
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, 
         evar1=0, evar2=0, eS=0.01, eT=0.01, ePt=0, eSit=0, epK=c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02),
         eBt=0.01,
         method="ga", r=0, runs=10000, 
         k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", warn="y", eos="eos80", long=1.e20, lat=1.e20)
{
  # if the concentrations of total silicate and total phosphate are NA
  # they are set to 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  
  # Input checking
    # --------------
    
    if (! method %in% c("ga", "mo", "mc"))
        stop ("Invalid input parameter: ", method)

    # Only two options for eos
    if (eos != "teos10" && eos != "eos80")
        stop ("invalid parameter eos: ", eos)

    # Input conditioning
    # -------------------
    
    n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), 
             length(evar1), length(evar2), length(r), length(eS), length(eT), length(ePt), length(eSit),
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
    if(length(evar1)!=n){evar1 <- rep(evar1[1],n)}
    if(length(evar2)!=n){evar2 <- rep(evar2[1],n)}
    if(length(r)!=n){r <- rep(r[1],n)}
    if(length(eS)!=n){eS <- rep(eS[1],n)}
    if(length(eT)!=n){eT <- rep(eT[1],n)}
    if(length(ePt)!=n){ePt <- rep(ePt[1],n)}
    if(length(eSit)!=n){eSit <- rep(eSit[1],n)}
    if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
    if(length(kf)!=n){kf <- rep(kf[1],n)}
    if(length(ks)!=n){ks <- rep(ks[1],n)}
    if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
    if(length(b)!=n){b <- rep(b[1],n)}

    # Check sign of all input errors
    neg_evar1 <- evar1 < 0
    evar1[neg_evar1] <- -evar1[neg_evar1]
    neg_evar2 <- evar2 < 0
    evar2[neg_evar2] <- -evar2[neg_evar2]
    neg_eS <- eS < 0
    eS[neg_eS] <- -eS[neg_eS]
    neg_eT <- eT < 0
    eT[neg_eT] <- -eT[neg_eT]
    neg_ePt <- ePt < 0
    ePt[neg_ePt] <- -ePt[neg_ePt]
    neg_eSit <- eSit < 0
    eSit[neg_eSit] <- -eSit[neg_eSit]
    
    # if eBt=NULL, set eBt equal to zero
    if(is.null(eBt)) {eBt = 0.0}

    # if epK=NULL, set all pK errors to zero
    if(is.null(epK)) {epK = rep(0, 7)}
  
    # Default value for epK
    if (missing(epK))
    {
        epK <- c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02)
    }
    else
    {
        # Check validity of epK
        if (length(epK) == 1 && epK == 0)
        {
            # this means that the caller does not want to account for errors on dissoc. constants
            epK <- rep(0, 7)
        }
        else if (length(epK) != 7)
            stop ("invalid parameter epK: ", epK)
        else
        {
            # Check sign of given epK
            neg_epK <- epK < 0
            epK[neg_epK] <- -epK[neg_epK] 
        }
    }
    
    if (method == "ga")
    {
        r0 <- r*0.0
        errs <- .errors_ga (flag, var1, var2, S, T, Patm, P, Pt, Sit, evar1, evar2, r0, eS, eT, ePt, eSit,
                epK, eBt, k1k2, kf, ks, pHscale, b, gas, warn, eos=eos, long=long, lat=lat)
    }
    else if (method == "mo")
    {
        errs <- .errors_ga (flag, var1, var2, S, T, Patm, P, Pt, Sit, evar1, evar2, r, eS, eT, ePt, eSit,
                epK, eBt, k1k2, kf, ks, pHscale, b, gas, warn, eos=eos, long=long, lat=lat)
    }
    else if (method == "mc")
    {
        errs <- .errors_mc (flag, var1, var2, S, T, Patm, P, Pt, Sit, evar1, evar2, eS, eT, ePt, eSit,
                epK, eBt, k1k2, kf, ks, pHscale, b, gas, runs, warn, eos=eos, long=long, lat=lat)
    }

    return (errs)
}

#==============================================================================================
#                                                                                             #
#    1st method : Gaussian (r=0) or Method of moments (r nonzero, i.e. with covariance terms) #
#                                                                                             #
#==============================================================================================

# .errors_ga()
#
# This routine esiimates uncertainties in computed carbonate system variables 
# by propagating errors (uncertainties) in the six input variables, including 
#  - the pair of carbonate system variables, 
#  - the 2 nutrients (total dissolved inorganic silicon and phosphorus concentrations), and
#  - temperature and salinity, as well as
# the errors in dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc.
#
# It computes numerical derivatives then applies them to propagate errors using the Gaussian method
# The Gaussian method is the standard technique for estimating the second moment of a computed variable z
# (its variance or standard deviation) based on a first-order approximation to z.
# The Gaussian method is the same as the method of moments but without covariance terms.
#
# Input parameters :
#   - evar1, evar2   :  standard error (uncertainty) in var1 and in var2 of input pair of carbonate system variables
#   - r              :  correlation coefficient between var1 and var2 
#   - eS, eT         :  standard error (uncertainty) in salinity and in temperature
#   - ePt, eSit      :  standard error (uncertainty) in total dissolved inorganic phosphorus and in total dissolved inorganic silicon concentrations
#   - epK            :  standard error (uncertainty) in the 7 dissociation constants mentioned above (a vector)
#   - others         :  same as input of subroutine  carb() : scalar or vectors
#
# All parameters may be scalars or vectors except epK and gas.
#   * gas must be a scalar
#   * epK must be vector of 7 values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#     these errors are assumed to be the same for all input data.
#
# In constrast, evar1, evar2, eS, eT, ePt and eSit, 
#   - if vectors, are errors associated with each data point
#   - if scalars, are one error value associated to all data points
#
# Returns a 2-dimensional data-frame, with the folowing columns :
# - H      total error to [H+] concentration  (mol/kg)
# - pH     total error to pH
# - CO2    total error to CO2 concentration (mol/kg)
# - fCO2   total error to "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (µatm)
# - pCO2   total error to "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (µatm)
# - HCO3   total error to HCO3 concentration (mol/kg)
# - CO3    total error to CO3 concentration (mol/kg)
# - DIC    total error to DIC concentration (mol/kg)
# - ALK    total error to ALK, total alkalinity (mol/kg)
# - OmegaAragonite  total error of Omega aragonite (aragonite saturation state)
# - OmegaCalcite    total error of Omega calcite   (calcite saturation state)

.errors_ga <- 
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, evar1=0, evar2=0, r=0, eS=0.01, eT=0.01,
         ePt=0, eSit=0, epK=NULL, eBt=NULL, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential",
         warn="y",
         eos="eos80", long=1.e20, lat=1.e20)
{
    # names of dissociation constants
    Knames <- c ('K0','K1','K2','Kb','Kw','Kspa', 'Kspc')

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

    # Convert error on pH to error on [H+] concentration
    isH <- (var1name == 'H')
    r[isH]  <- -1.0*r[isH]       # Inverse sign of 'r' if var1 is pH
    pH <- var1[isH]
    epH <- evar1[isH]       # Error on pH
    H  <- 10**(-pH)         # H+ concentration

    # dpH = d(-log10[H])
    #     = d(- ln[H] / ln[10] )
    #     = -(1/ln[10]) * d (ln[H])
    #     = -(1/ln[10]) * (dH / H)
    # Thus dH = - ln[1O] * [H] dpH
    eH <-  log(10) * H * epH   # the formula is del H = -log(10) H epH, but we drop minus sign here (errors are always positive)
    evar1[isH] <- eH

    # Same conversion for second variable
    isH <- (var2name == 'H')
    r[isH]  <- -1.0*r[isH]       # Inverse sign of 'r' if var2 is pH
    pH <- var2[isH]
    epH <- evar2[isH]       # Error on pH
    H  <- 10**(-pH)         # H+ concentration
        
    eH <-  log(10) * H * epH   
    evar2[isH] <- eH

    # initialise total square error
    n <- length(flag)
    sq_err <- numeric(n)
        
    # Contribution of var1 to squared standard error
    if (any (evar1 != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv1 <- derivnum ('1', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                            pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv1 * evar1
        sq_err <- sq_err + err * err
    }

    # Contribution of var2 to squared standard error
    if (any (evar2 != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv2 <- derivnum ('2', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                            pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv2 * evar2
        sq_err <- sq_err + err * err
    }

    # Contribution of covariance of var1 and var2 to squared standard error
    if (any (r != 0.0) & any (evar1 != 0.0) & any (evar2 != 0.0))
    {
        # Compute covariance from correlation coeff. and standard deviations
        covariance <- r * evar1 * evar2
        # Contribution to squared error
        err2 <- 2 * deriv1 * deriv2 * covariance
        sq_err <- sq_err + err2
    }

    # Contribution of Silicon (total dissolved inorganic concentration) to squared standard error
    if (any (eSit != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv <- derivnum ('sil', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                           pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv * eSit
        sq_err <- sq_err + err * err
    }

    # Contribution of Phosphorus (total dissolved inorganic concentration) to squared standard error
    if (any (ePt != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv <- derivnum ('phos', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                           pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv * ePt
        sq_err <- sq_err + err * err
    }

    # Contribution of T (temperature) to squared standard error
    if (any (eT != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv <- derivnum ('T', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                           pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv * eT
        sq_err <- sq_err + err * err
    }

    # Contribution of S (salinity) to squared standard error
    if (any (eS != 0.0))
    {
        # Compute sensitivities (partial derivatives)
        deriv <- derivnum ('S', flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                           pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
        err <- deriv * eS
        sq_err <- sq_err + err * err
    }
    
    # Salinity and Temperature converted to EOS-80 (if necessary)
    SP    <- rep(NA, n)
    InsT  <- rep(NA, n)
 
    # if use of TEOS-10 standard
    if (eos == "teos10")
    {
        # Must convert temperature and salinity from TEOS-10 to EOS-80
        # convert temperature: from Conservative (CT) to in-situ temperature
        # and salinity from Absolute to Practical (SP)
        STeos <- teos2eos_geo (S, T, P, long, lat)
        InsT <- STeos$T
        SP <- STeos$SP
    }
    else
    {
        InsT <- T
        SP <- S
    }

    # Preliminary calculations for dissociation constants
    if (any (epK != 0))
    {

        # Ks (free pH scale) at zero pressure and given pressure
        Ks_P0 <- Ks(S=SP, T=InsT, P=0, ks=ks, warn=warn)
        Ks    <- Ks(S=SP, T=InsT, P=P, ks=ks, warn=warn)

        # Kf on free pH scale
        Kff_P0 <- Kf(S=SP, T=InsT, P=0, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
        Kff <- Kf(S=SP, T=InsT, P=P, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
        # Kf on given pH scale
        Kf <- Kf(S=SP, T=InsT, P=P, pHscale=pHscale, kf=kf, Ks_P0, Ks, warn=warn)

        # Conversion factor from total to SWS pH scale at zero pressure
        ktotal2SWS_P0 <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks_P0,Kff=Kff_P0,warn=warn)$ktotal2SWS

        # Conversion factor from SWS to chosen pH scale
        conv <- kconv(S=SP,T=InsT,P=P,kf=kf,Ks=Ks,Kff=Kff,warn=warn)
        kSWS2chosen <- rep(1.,n)
        kSWS2chosen [pHscale == "T"] <- conv$kSWS2total [pHscale == "T"]
        kSWS2chosen [pHscale == "F"] <- conv$kSWS2free [pHscale == "F"]  
    }
        
    # Contribution of all pKi to squared standard error
    for (i in 1:length(epK))
    {
        # if error on Ki is given
        if (epK[i] != 0.0)
        {
            # Compute Ki
            Ki <- switch (i,
                          K0(S=SP, T=InsT, Patm=Patm, P=0, warn=warn),
                          K1(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn),
                          K2(S=SP, T=InsT, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0, warn=warn),
                          Kb(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0, warn=warn),
                          Kw(S=SP, T=InsT, P=P, pHscale=pHscale, kSWS2chosen, warn=warn),
                          Kspa(S=SP, T=InsT, P=P, warn=warn),
                          Kspc(S=SP, T=InsT, P=P, warn=warn)
                          )
            # compute error on Ki from that on pKi
            eKi <- - epK[i] * Ki * log(10)

            # Compute sensitivities (partial derivatives)
            # No need to pass option "eos" since conversion to eos has been done already
            deriv <- derivnum (Knames[i], flag, var1, var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                pHscale=pHscale, b=b, gas=gas, warn=warn)
            err <- deriv * eKi
            sq_err <- sq_err + err * err
        }
    }

    # Contribution of Total Boron (computed from Bt/S ratio) to squared standard error
    if (eBt != 0.0)
    {
        # Compute sensitivities (partial derivatives)
        deriv <- derivnum ('bor', flag, var1, var2, S=SP, T=InsT, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, 
                           pHscale=pHscale, b=b, gas=gas, warn=warn)
        # err is derivative * absolute error in boron (i.e., eBt=0.01 is a 1% error)
        err <- deriv * bor(S=S, b=b) * eBt  
        sq_err <- sq_err + err * err
    }

    
    # Compute resulting total error (or uncertainty)
    error <- sqrt (sq_err)

    return (error)
}


#===========================================================================================
#                                                                                          #
#    2nd method : Monte Carlo                                                              #
#                                                                                          #
#===========================================================================================

# Function that computes standard deviation 
# on first dimension of a given 3D matrix "x"
.Sd_3D_1stD <- function(x)
{          
    first_dim <- dim(x)[1]
    # Compute mean over first dimension
    x_mean <- colMeans(x)
    # x_mean is a 2D array
    # Extend array x_mean to 3D by replication
    new_dim <- c(dim(x_mean), first_dim)
    x_mean <- array (x_mean, new_dim)
    # Transpose extended x_mean so that it conforms to array x
    x_mean <- aperm (x_mean, c(3,1,2))

    # Compute mean-centered array
    y <- x - x_mean
    # Compute variance estimator
    var <- colMeans(y^2)*(first_dim/(first_dim-1))

    return (sqrt(var))
}


# Function that generates deviate values for dissoc. constants Kx from given error on pKx (and same for Bt from eBt)
#
# Special case for K0 :
#    This function generates a set of small deltas departing from 0 whose distribution is close to normal
#    These deltas will then be used to generate deviate samples of K0, from standard value of K0
#
# Other cases : K1, K2, Kb, Kw, Kspa and Kspc
#    This function generates a set of deviate values for dissociation constant Kx
#    Distribution of that set is close to a Normal distribution centered on standard value of Kx
#
# Input parameters :
#   - epK      :  standard error (or uncertainty) on all seven dissociation constants (a vector of length 7)
#   - S, T, P  :  Salinity, Temperature and Pressure, vectors of same length n
#   - Patm     :  Surface atmospheric pressure in atm, vector of length n
#   - pHscale  :  pH scale, character ("T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale)
#                 vector of length n
#   - k1k2     :  option for dissociation constants K1 and K2, vector of length n
#   - kf, ks   :  options for dissociation constants Kf and Ks, vectors of length n
#   - runs     :  number of runs of Monte Carlo (= number of deviates per sample)
#
# Returns a 2-dimensional data-frame, with 7 columns: K0, K1, K2, Kb, Kw, Kspa and Kspc 
#   each column contains deviate delta values of one dissociation constant
#   column length is n * runs
#
.gen_delta_Kx <- function (epK, eBt, S, T, P, Patm, pHscale, k1k2, kf, ks, b, runs, warn="y")
{
    n <- length(S)

    # For convenience, add eBt (fractional absolute err in Bt) to
    # end of epK vector (errors on K values in terms of pK)
    epKplus <- c(epK, eBt)

    # names of dissociation constants
    Kplusnames <- c ('K0','K1','K2','Kb','Kw','Kspa', 'Kspc', 'bor')

    # Initalise output data frame
    nrows <- n * runs
    delta_Kx <- data.frame (matrix( numeric(0), nrow=nrows, ncol=0))

    # Devise a function that generates simulation samples
    # for one variable with central value "val" and standard error "std_err"
    gen_sim <- function (val, std_err)
    {
        # if no standard error"
        if (all(std_err == 0))
            # Replicate variable value
            sim_array <- rep(val, times=runs)
        else
        {
            # Generate deviates from value using Normal distribution
            sim_array <- rnorm(runs, val, std_err)
            # Set to zero all negative values
            sim_array[sim_array < 0] = 0.0
        }
        return(sim_array)
    }

    # Preliminary calculations for dissociation constants

    # Ks (free pH scale) at zero pressure and given pressure
    Ks_P0 <- Ks(S=S, T=T, P=0, ks=ks, warn=warn)
    Ks    <- Ks(S=S, T=T, P=P, ks=ks, warn=warn)

    # Kf on free pH scale
    Kff_P0 <- Kf(S=S, T=T, P=0, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
    Kff <- Kf(S=S, T=T, P=P, pHscale="F", kf=kf, Ks_P0, Ks, warn=warn)
    # Kf on given pH scale
    Kf <- Kf(S=S, T=T, P=P, pHscale=pHscale, kf=kf, Ks_P0, Ks, warn=warn)

    # Conversion factor from total to SWS pH scale at zero pressure
    ktotal2SWS_P0 <- kconv(S=S,T=T,P=P,kf=kf,Ks=Ks_P0,Kff=Kff_P0,warn=warn)$ktotal2SWS

    # Conversion factor from SWS to chosen pH scale
    conv <- kconv(S=S,T=T,P=P,kf=kf,Ks=Ks,Kff=Kff, warn=warn)
    kSWS2chosen <- rep(1.,n)
    kSWS2chosen [pHscale == "T"] <- conv$kSWS2total [pHscale == "T"]
    kSWS2chosen [pHscale == "F"] <- conv$kSWS2free [pHscale == "F"]  

    # Convert error on pKi to error on Ki
    for (i in 1:length(epKplus))
    {
        # Compute Ki
        Ki <- switch (i,
                      seacarb::K0(S=S, T=T, Patm=Patm, P=0),
                      seacarb::K1(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0),
                      seacarb::K2(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0),
                      seacarb::Kb(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0),
                      seacarb::Kw(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen),
                      seacarb::Kspa(S=S, T=T, P=P),
                      seacarb::Kspc(S=S, T=T, P=P),
                      seacarb::bor(S=S, b=b)
                      )
        if (i == 1)
            center_value = 0.0    # special case for K0
        else
            center_value = Ki

        # if error on Ki is given
        if (epKplus[i] != 0.0)
        {
            # compute error (not signed) on Ki from that on pKi
            if ( i == 8 ) 
            {
                eKi <-  epKplus[i] * Ki
            } else {
                eKi <-  epKplus[i] * Ki * log(10)
            }
            # Generate deviate values for Ki or deltas for K0
            spl_Ki <- mapply (gen_sim, center_value, eKi)
            # Reshape samples from 2D matrix to vector
            dim(spl_Ki) <- c(n*runs)
        }
        else
        {
            # Simply replicate value of Ki or  replicate 0.0 for K0
            spl_Ki <- rep (center_value, each=runs)
        }

        # Store (deviate) values for Ki in a data frame column
        Kname <- Kplusnames[i]
        delta_Kx[[Kname]] <- spl_Ki
    }
    return (delta_Kx)
}


# errors()
# This subroutine does error propagation on the computation of carbonate system variables 
# from errors (or uncertainties) on six input 
#  - pair of carbonate system variables 
#  - nutrients (total inorganic silicon and phosphorus concentrations)
#  - temperature and salinity
# plus errors on dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#
# It computes standard error on carbonate system output variables using Monte Carlo method.
#
# Input parameters :
#   - evar0, evar2   :  standard error (or uncertainty) on var1 and var2 of input pair of carbonate system variables
#   - eS, eT         :  standard error (or uncertainty) on Salinity and Temperature
#   - ePt, eSit      :  standard error (or uncertainty) on Phosphorus and Silicon total inorganic concentrations
#   - epK            :  standard error (or uncertainty) on all seven dissociation constants (a vector)
#   - eBt            :  standard error (or uncertainty) on total boron (relative factional error, e.g. eBt=0.01 is a 1% error in total boron

#   - runs           :  number of runs of Monte Carlo (= number of simulated samples)
#                       default is 10000
#   - others         :  same as input of subroutine  carb() : scalar or vectors
#
# All parameters may be scalars or vectors except epK, eBt, method, runs, and gas.
#   * runs must be a scalar
#   * gas and method must each be a character string
#   * epK must be a vector of seven values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
#   * eBt must be a scalar : error in total boron, e.g. eBt=0.01 is a 1% error 
#     These 3 types of errors are assumed to be the same for all input data.
#
# In constrast, for evar1, evar2, eS, eT, ePt and eSit:
#   - if they are vectors, they represent errors associated with each data point
#   - if they are scalars, they represent one error value each associated to all data points
#
# Returns a 2-dimensional data-frame, with the folowing columns :
# - H      total error to [H+] concentration  (mol/kg)
# - pH     total error to pH
# - CO2    total error to CO2 concentration (mol/kg)
# - fCO2   total error to "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (µatm)
# - pCO2   total error to "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (µatm)
# - HCO3   total error to HCO3 concentration (mol/kg)
# - CO3    total error to CO3 concentration (mol/kg)
# - DIC    total error to DIC concentration (mol/kg)
# - ALK    total error to ALK, total alkalinity (mol/kg)
# - OmegaAragonite  total error of Omega aragonite (aragonite saturation state)
# - OmegaCalcite    total error of Omega calcite   (calcite saturation state)
#
# Remarks :
#   Accuracy is strongly depending on number of runs. Default run number is 10,000, 
#   You may want to choose a number of 100,000 runs or more, but beware that such calculations may take several seconds per data point.
#
.errors_mc <- 
function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, evar1=0, evar2=0, eS=0.01, eT=0.01, ePt=0, eSit=0,
         epK=NULL, eBt=NULL, k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", runs=10000, warn="y",
         eos="eos80", long=1.e20, lat=1.e20)
{
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

    # if use of EOS-10 standard
    if (eos == "teos10")
    {
        # Must convert temperature and salinity from TEOS-10 to EOS-80
        # convert temperature: from Conservative (CT) to in-situ temperature
        # and salinity from Absolute to Practical (SP)
        STeos <- teos2eos_geo (S, T, P, long, lat)
        InsT <- STeos$T
        SP <- STeos$SP
    }
    else
    {
        InsT <- T
        SP <- S
    }
    
    # Devise a function that generates simulation samples 
    # for one variable with central value "val" and standard error "std_err"
    gen_sim <- function (val, std_err)
    {
        # if no standard error"
        if (all(std_err == 0))
            # Replicate variable value
            sim_array <- rep(val, times=runs)
        else
        {
            # Generate deviates from value using Normal distribution
            sim_array <- rnorm(runs, val, std_err)
            # Set to zero all negative values
            sim_array[sim_array < 0] = 0.0
        }
        return(sim_array)
    }

    # Generate deviate sample values for var1
    spl_var1 <- mapply (gen_sim, var1, evar1)
    # Same for var2, Salinity, Temperature, total inorganic Phosphorus and Silicon
    spl_var2 <- mapply (gen_sim, var2, evar2)
    spl_S    <- mapply (gen_sim, S, eS)
    spl_T    <- mapply (gen_sim, T, eT)
    spl_Pt   <- mapply (gen_sim, Pt, ePt)
    spl_Sit  <- mapply (gen_sim, Sit, eSit)

    # Reshape samples from 2D matrix to vector
    n <- length(var1)
    dim(spl_var1) <- c(n*runs)
    dim(spl_var2) <- c(n*runs)
    dim(spl_S)    <- c(n*runs)
    dim(spl_T)    <- c(n*runs)
    dim(spl_Pt)   <- c(n*runs)
    dim(spl_Sit)  <- c(n*runs)
        
    # Generate deviate delta values for K0
    # Generate deviate values for other dissoc. constants Kx
    # Note : use of Salinity and Temperature converted to EOS-80 as they are when routine carb() computes dissociation constants  
    spl_Kx <- .gen_delta_Kx (epK, eBt, S=SP, T=InsT, P, Patm, pHscale, k1k2, kf, ks, b, runs, warn=warn)

    # All other parameters and variables
    spl_flag <- rep (flag, each=runs)
    spl_Patm <- rep (Patm, each=runs)
    spl_P    <- rep (P, each=runs)
    spl_k1k2 <- rep (k1k2, each=runs)
    spl_kf   <- rep (kf, each=runs)
    spl_ks   <- rep (ks, each=runs)
    spl_pHscale <- rep (pHscale, each=runs)
    spl_b    <- rep (b, each=runs)

    # Define local functions K0, K1, K2, ...
    # --------------------------------------
    # Special case for K0
    #
    # It computes original K0 values then add precalculated deltas to generate 
    # and return deviate values
    K0 <- function(S=35,T=25,P=0,Patm=1,...)
    #K0 <- function(S=35,T=25,P=0,Patm=1,warn='n')
    {
        # Call original K0 function
       #out <- seacarb::K0(S, T, P, Patm=Patm, warn=warn)
        out <- seacarb::K0(S, T, P, Patm=Patm)
        # perturb value of K0 by adding deltas
        out = out + spl_Kx$K0
        return (out)
    }
    # General case : all other
    # 'function' returns precalculated deviate values
    K1 <- function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0, warn="y")  spl_Kx$K1
    K2 <- function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0, warn="y")  spl_Kx$K2
    Kw <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0, warn="y")  spl_Kx$Kw
    Kb <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0, warn="y")  spl_Kx$Kb
    Kspa <- function(S=35,T=25,P=0, warn="y")  spl_Kx$Kspa
    Kspc <- function(S=35,T=25,P=0, warn="y")  spl_Kx$Kspc
    bor <- function(S=35,b="u74")  spl_Kx$bor

    # Note : in the general case (K1, K2,...) the function Kx is called once by the function carb()
    #        We can calculate in anticipation deviate values and substitute real values with deviate ones
    #        at run time.
    #        In the special case of K0, function K0() is called twice by function carb() 
    #        with two different parameter sets so that we cannot precalculate deviate values.
    #        We must generate deviate values at run time, when K0() is called.
    #        This is the purpose of locally defined K0() function
    
    # Save then change execution environment of function carb()
    saved_env <- environment(carb)
    # Using environment(K0) will hide original seacarb::K0, seacarb::K1, ... functions
    environment(carb) <- environment(NULL)
    # Note :  environment(NULL) is the environment created when this function .errors_mc() is executed
    #         It then contains all locally defined functions, like K0, K1, K2, ...
    #         Then, when carb() is executed and looks for a name (variable or function) that is not a local variable,
    #         it will look into evironment(NULL), which has been attached to it, 
    #         before looking into global or seacarb module environment.

    # Compute output carbonate system variables
    seacarb = carb(spl_flag, spl_var1, spl_var2, S=spl_S, T=spl_T, Patm=spl_Patm, P=spl_P, Pt=spl_Pt, Sit=spl_Sit, 
                spl_k1k2, spl_kf, spl_ks, spl_pHscale, spl_b, warn=warn, eos=eos, long=long, lat=lat)

    # Restore environment of carb()
    environment(carb) <- saved_env

    # Add one column for [H+]
    H <- 10^(-seacarb$pH)
    seacarb <- cbind(H,seacarb)

    # Drop unnecessary columns
    drops=c('flag', 'S', 'T', 'P', 'Patm', 'fCO2pot', 'pCO2pot', 'fCO2insitu', 'pCO2insitu')
    seacarb <- seacarb[,!(names(seacarb) %in% drops)]
        
    # if all input pairs are of same type
    if (all(flag == flag[1]))
    {
        # Drop input pair
        var1name <- varnames[flag[1],1]
        var2name <- varnames[flag[1],2]
        drops=c(var1name, var2name)
        # if to drop H, drop also pH
        if (var1name == 'H' || var2name == 'H')   drops = c(drops, 'pH')
        seacarb <- seacarb[,!(names(seacarb) %in% drops)]
    }
    
    # if more than one data point
    if (n > 1)
    {
        # Reshape results to a matrix
        seacarb_3D <- data.matrix(seacarb)
        # Reshape to a 3D array
        nvars <- dim(seacarb)[2]
        dim(seacarb_3D) <- c(runs, n, nvars)

        # Compute standard deviation on all simulated samples
        std_dev <- .Sd_3D_1stD (seacarb_3D)
        # Convert to data frame
        std_dev <- data.frame (std_dev)
        colnames(std_dev) <- colnames(seacarb)
    }
    else
    {
        # devise a function that computes standard deviation on columns of a given 2D matrix "x"
        colSd_2D <- function(x)sqrt(rowMeans((t(x)-colMeans(x))^2)*((dim(x)[1])/(dim(x)[1]-1)))

        # Compute standard deviation on all simulated samples
        std_dev <- data.frame(t(colSd_2D(seacarb)))
    }
        
    return (std_dev)
}
