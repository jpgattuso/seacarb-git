# derivnum()
# This subroutine computes partial derivatives of output carbonate variables 
# with respect to input variables (two), plus nutrients (two), temperature and salinity,
# and dissociation constants.
#
# It uses the central differences, which consists of  
# - introducing a small perturbation delta (plus or minus) in one input
# - and computing the induced delta in output variables
#
# The ratio between delta and input value is chosen as follows: 
#    for perturbing input pair :  a constant which depends on the type of input pair of variables.
#    for perturbing nutrients  :  1.e-6
#    for perturbing T and S    :  1.e-3
#
# Input parameters :
#   - varid  :  (scalar) identifier of variable with respect to which derivative is requested
#                = variable length, case insensitive, character code
#                  case '1' or 'var1'  :  Variable 1 of the input pair (This is TAlk if flag is 15)
#                  case '2' or 'var2'  :  Variable 2 of the input pair (This is DIC  if flag is 15)
#                  case 'sil', 'silt', 'tsil' or 'silicate'      : Total dissolved inorganic silicon concentration
#                  case 'phos', 'phost', 'tphos' or 'phosphate'  : Total dissolved inorganic phosphorus concentration
#                  case 't', 'temp' or 'temperature' : temperature
#                  case 's', 'sal' or 'salinity'     : salinity
#                  case 'K0','K1','K2','Kb','Kw','Kspa' or 'Kspc' : dissociation constant
#                  case 'bor' : total boron
#
#   - others     :  same à input of subroutine  carb() : scalar or vectors
#
# Returns one set of partial derivatives
derivnum <- 
function(varid, flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, 
         k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", warn="y",  eos="eos80", long=1.e20, lat=1.e20)
{
    # Input conditionning
    # -------------------
    
    n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), length(ks), length(b))
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

    # Only two options for eos
    if (eos != "teos10" && eos != "eos80")
        stop ("invalid parameter eos: ", eos)

    varid <-toupper(varid)
    
    # Constants: names of variables, sorted by flag number
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
    
    # Relative deltas (perturbations) sorted by flag number
    # default value is 1.e-3
    const_deltas = rbind (
        c(1.e-8, 1.e-3),         # flag = 1    
        c(1.e-3, 1.e-3),         # flag = 2    
        c(1.e-3, 1.e-3),         # flag = 3    
        c(1.e-6, 1.e-6),         # flag = 4    
        c(1.e-5, 1.e-5),         # flag = 5    
        c(1.e-8, 1.e-3),         # flag = 6    
        c(1.e-8, 1.e-3),         # flag = 7    
        c(1.e-1, 1.e-5),         # flag = 8    
        c(1.e-1, 1.e-1),         # flag = 9    
        c(1.e-3, 1.e-3),         # flag = 10   
        c(1.e-3, 1.e-3),         # flag = 11   
        c(1.e-3, 1.e-3),         # flag = 12   
        c(1.e-3, 1.e-3),         # flag = 13   
        c(1.e-3, 1.e-3),         # flag = 14   
        c(1.e-6, 1.e-6),         # flag = 15   
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(1.e-3, 1.e-8),         # flag = 21   
        c(1.e-3, 1.e-3),         # flag = 22   
        c(1.e-3, 1.e-3),         # flag = 23   
        c(1.e-6, 1.e-6),         # flag = 24   
        c(1.e-5, 1.e-5)          # flag = 25   
    )

    # Replace all values of const_deltas (above) with 1e-6 for consistency with CO2SYS-MATLAB
    const_deltas = rbind (
        c(1.e-6, 1.e-6),         # flag = 1    
        c(1.e-6, 1.e-6),         # flag = 2    
        c(1.e-6, 1.e-6),         # flag = 3    
        c(1.e-6, 1.e-6),         # flag = 4    
        c(1.e-6, 1.e-6),         # flag = 5    
        c(1.e-6, 1.e-6),         # flag = 6    
        c(1.e-6, 1.e-6),         # flag = 7    
        c(1.e-6, 1.e-6),         # flag = 8    
        c(1.e-6, 1.e-6),         # flag = 9    
        c(1.e-6, 1.e-6),         # flag = 10   
        c(1.e-6, 1.e-6),         # flag = 11   
        c(1.e-6, 1.e-6),         # flag = 12   
        c(1.e-6, 1.e-6),         # flag = 13   
        c(1.e-6, 1.e-6),         # flag = 14   
        c(1.e-6, 1.e-6),         # flag = 15   
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(1.e-6, 1.e-6),         # flag = 21   
        c(1.e-6, 1.e-6),         # flag = 22   
        c(1.e-6, 1.e-6),         # flag = 23   
        c(1.e-6, 1.e-6),         # flag = 24   
        c(1.e-6, 1.e-6)          # flag = 25   
    )
    
    # Reference 'concentrations': sorted by flag number
    const_refs  = rbind (
        c(1e-8, 10e-6),          # flag = 1    
        c(10e-6, 1800e-6),       # flag = 2    
        c(10e-6, 100e-6),        # flag = 3    
        c(10e-6, 2300e-6),       # flag = 4    
        c(10e-6, 2000e-6),       # flag = 5    
        c(1e-8, 1800e-6),        # flag = 6    
        c(1e-8, 100e-6),         # flag = 7    
        c(1e-8, 2300e-6),        # flag = 8    
        c(1e-8, 2000e-6),        # flag = 9    
        c(1800e-6, 100e-6),      # flag = 10   
        c(1800e-6, 2300e-6),     # flag = 11   
        c(1800e-6, 2000e-6),     # flag = 12   
        c(100e-6, 2300e-6),      # flag = 13   
        c(100e-6, 2000e-6),      # flag = 14   
        c(2300e-6, 2000e-6),     # flag = 15 
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(0, 0),
        c(400, 1e-8),            # flag = 21   
        c(400, 1800e-6),         # flag = 22   
        c(400, 100e-6),          # flag = 23   
        c(400, 2300e-6),         # flag = 24   
        c(400, 2000e-6)          # flag = 25   
    )

    # Compute two slightly different values for input
    # -----------------------------------------------

    # Default values : not perturbed
    # no change on var1 and var2
    var11 = var12 = var1
    var21 = var22 = var2
    # no change on Sil total and Phos total
    Sit1 = Sit2 = Sit
    Pt1  = Pt2  = Pt
    # no change on T and S
    T1 = T2 = T
    S1 = S2 = S

    # Initialise absolute deltas
    abs_dx <- rep(NA, n)

    # Uppercase names of dissociation constants and 'total boron'
    K_id <- c('K0', 'K1', 'K2', 'KB', 'KW', 'KSPA', 'KSPC', 'BOR')
    # Flag for dissociation constant as perturbed variable 
    flag_dissoc_K = varid %in% K_id
    
    # if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
    {
        # Approximate values for K0, K1, K2, Kb, Kw, Kspa, Kspc, and bor
        # They will be used to compute an absolute perturbation value on these constants
        K <- c(0.034, 1.2e-06, 8.3e-10, 2.1e-09, 3.1e-14, 6.7e-07, 4.3e-07, 4.156877e-4)
        # Choose value of absolute perturbation
        index <- which (K_id == varid)
        perturbation = K[index] * 1.e-3   # 0.1 percent of Kx value
        abs_dx = 2 * perturbation
    }
    # Var1 is perturbed
    else if (varid %in% c('1', 'VAR1'))
    {
        # Define relative deltas
        delta = const_deltas[flag, 1]
        xref  = const_refs[flag, 1]

        # where input variable is [H+] concentration
        FH = (var1name == 'H')
        # Change slightly [H+]
        H_var1 <- 10^(-var1[FH])
        H1 = H_var1 - xref[FH]*delta[FH]
        H2 = H_var1 + xref[FH]*delta[FH]
        var11[FH] = -log10(H1)
        var12[FH] = -log10(H2)
        abs_dx[FH] = H2 - H1

        # where it is not
        GH = ! FH
        # Change slightly var1
        var11[GH] = var1[GH] - xref[GH]*delta[GH]
        var12[GH] = var1[GH] + xref[GH]*delta[GH]
        abs_dx[GH] = var12[GH] - var11[GH]
    }
    else if (varid %in% c('2', 'VAR2'))    # var2 (second variable of input pair) is perturbed
    {
        # Define relative deltas
        delta = const_deltas[flag, 2]
        xref  = const_refs[flag, 2]

        # where second input variable is [H+] concentration
        FH = (var2name == 'H')
        # Change slightly H+]
        H_var2 <- 10^(-var2[FH])
        H1 = H_var2 - xref[FH]*delta[FH]
        H2 = H_var2 + xref[FH]*delta[FH]
        var21[FH] = -log10(H1)
        var22[FH] = -log10(H2)
        abs_dx[FH] = H2 - H1

        # where it is not
        GH = ! FH
        # Change slightly var2
        var21[GH] = var2[GH] - xref[GH]*delta[GH]
        var22[GH] = var2[GH] + xref[GH]*delta[GH]
        abs_dx[GH] = var22[GH] - var21[GH]
    }
    else if (varid %in% c('SIL', 'TSIL', 'SILT', 'SILICATE', 'SIT'))    # Sil total
    {
        # Define relative delta
        delta  = 1.e-3
        Sitref = 7.5e-6 # mol/kg (Global surface mean)

        center_value = Sit
        # Change slightly total silicon
        Sit1 = Sit - Sitref*delta
        Sit2 = Sit + Sitref*delta
        abs_dx = Sit2 - Sit1
    }
    else if (varid %in% c('PHOS', 'TPHOS', 'PHOST', 'PHOSPHATE', 'PT'))    # Phos total
    {
        # Define relative delta
        delta = 1.e-3
        Ptref = 0.5e-6   # mol/kg (Global surface mean)

        center_value = Pt
        # Change slightly total phosphorus
        Pt1 = Pt - Ptref*delta
        Pt2 = Pt + Ptref*delta
        abs_dx = Pt2 - Pt1
    }
    else if (varid %in% c('T', 'TEMP', 'TEMPERATURE'))    # Temperature
    {
        # Define relative delta
        delta = 1.e-4
        Tref = 18   # degrees Celcius
        
        center_value = T
        # Change slightly temperature
        T1 = T - Tref*delta
        T2 = T + Tref*delta
        abs_dx = T2 - T1
    }
    else if (varid %in% c('S', 'SAL', 'SALINITY'))   # Salinity
    {
        # Define relative delta
        delta = 1.e-4
        Sref = 35

        center_value = S
        # Change slightly salinity
        S1 = S - Sref*delta
        S2 = S + Sref*delta
        abs_dx = S2 - S1
    }
    else
    {
        stop ("Invalid input parameter: ", varid)
    }

    # Define local functions K0, K1, K2, ...
    # --------------------------------------

    # These function will replace original K0, K1, K2.... functions, at some computation time
    # They are intended to return slightly modified (perturbed) values of dissociation constants
    #
    # Their input (beside their regular parameters):
    #
    #   perturbation :  absolute value of (tiny) perturbation
    #   sign_factor  :  factor -1 or +1 gives the sign or perturbation
    if (varid == 'K0') 
    {
        K0 <- function(S=35,T=25,P=0,Patm=1,...)
       #K0 <- function(S=35,T=25,P=0,Patm=1,warn='n')
        {
            # Call original K0 function
            out <- seacarb::K0(S, T, P, Patm=Patm)
            #out <- K0(S, T, P, Patm=Patm, warn=warn)
            # perturb value of K0
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'K1') 
    {
        K1 <-function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0,...)
       #K1 <-function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0,warn='n')
        {
            # Call original K1 function
            out <- seacarb::K1(S, T, P, k1k2=k1k2, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0)
           #out <- K1(S, T, P, k1k2=k1k2, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0, warn=warn)
            # perturb value of K1
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'K2') 
    {
        K2 <- function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0,...)
       #K2 <- function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0,warn='n')
        {
            # Call original K2 function
            out <- seacarb::K2(S, T, P, k1k2=k1k2, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0)
           #out <- seacarb::K2(S, T, P, k1k2=k1k2, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0, warn=warn)
            # perturb value of K2
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'KW') 
    {
        Kw <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0,...)
       #Kw <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0, warn='n')
        {
            # Call original Kw function
            out <- seacarb::Kw(S, T, P, pHscale=pHscale, kSWS2scale=kSWS2scale)
           #out <- Kw(S, T, P, pHscale=pHscale, kSWS2scale=kSWS2scale, warn=warn)
            # perturb value of Kw
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'KB')
    {
        Kb <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0,...)
       #Kb <- function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0, warn='n')
        {
            # Call original Kb function
            out <- seacarb::Kb(S, T, P, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0)
           #out <- Kb(S, T, P, pHscale=pHscale, kSWS2scale=kSWS2scale, ktotal2SWS_P0=ktotal2SWS_P0, warn=warn)
            # perturb value of Kb
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'KSPA')
    {
        Kspa <- function(S=35,T=25,P=0,...)
       #Kspa <- function(S=35,T=25,P=0, warn='n')
        {
            # Call original Kspa function
            out <- seacarb::Kspa(S, T, P)
           #out <- Kspa(S, T, P, warn=warn)
            # perturb value of Kspa
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'KSPC')
    {
        Kspc <- function(S=35,T=25,P=0,...)
       #Kspc <- function(S=35,T=25,P=0, warn='n')
        {
            # Call original Kspc function
            out <- seacarb::Kspc(S, T, P)
           #out <- Kspc(S, T, P, warn=warn)
            # perturb value of Kspc
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    else if (varid == 'BOR')
    {
        bor <- function(S=35,...)
        {
            # Call original bor function
            out <- seacarb::bor(S)
            # perturb value of bor
            out = out + sign_factor * perturbation  # sign_factor is +1 or -1
            return (out)
        }
    }
    
    # PERTURBATION:
    # -------------
    
    # if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
    {
        # Save then change execution environment of function carb()
        saved_env <- environment(carb)
        # Using environment(NULL) will hide original seacarb::K0, seacarb::K1, ... function (depending on "varid")
        environment(carb) <- environment(NULL)
        # Note :  environment(NULL) is the environment created when the function derivnum() is executed
        #         It then contains all locally defined functions, like K0, K1, K2, ...
        #         Then, when carb() is executed and looks for a name (variable or function) that is not a local variable,
        #         it will look into evironment(NULL), which has been attached to it, 
        #         before looking into global or seacarb module environment.
        
        # Point 1: (one dissociation constant is somewhat smaller)
        sign_factor = -1.0
        cdel1 <- carb(flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)

        # Point 2: (one dissociation constant is somewhat bigger)
        sign_factor = 1.0
        cdel2 <- carb(flag, var1, var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)

        # Restore environment of carb()
        environment(carb) <- saved_env
    }
    else
    {
        # Point 1: (one of var1, var2, T or S is somewhat smaller)
        cdel1 <- carb(flag, var11, var21, S=S1, T=T1, Patm=Patm, P=P, Pt=Pt1, Sit=Sit1, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)

        # Point 2: (one of var1, var2, T or S is somewhat bigger)
        cdel2 <- carb(flag, var12, var22, S=S2, T=T2, Patm=Patm, P=P, Pt=Pt2, Sit=Sit2, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b, gas=gas, warn=warn, eos=eos, long=long, lat=lat)
    }

    # Compute [H+] concentration and add it to data-frame
    H <- 10^(-cdel1$pH)
    cdel1 <- cbind(H,cdel1)
    H <- 10^(-cdel2$pH)
    cdel2 <- cbind(H,cdel2)

    # Centered difference
    dy <-(cdel2 - cdel1)

    # Drop unnecessary columns
    # ------------------------
    
    # list fixed columns to drop
    drops=c('flag', 'S', 'T', 'P', 'Patm', 'fCO2pot', 'pCO2pot', 'fCO2insitu', 'pCO2insitu')
    # add input variable 1 to list
    if (all(var1name == var1name[1]))
    {
        drops = c(drops, var1name)
        # if to drop H, drop also pH
        if (var1name[1] == 'H')
            drops = c(drops, 'pH')
    }
    # add input variable 2 to list
    if (all(var2name == var2name[1]))
    {
        drops = c(drops, var2name)
        # if to drop H, drop also pH
        if (var2name[1] == 'H')
            drops = c(drops, 'pH')
    }
    # Drop columns
    dy <- dy[,!(names(dy) %in% drops)]

    # Compute ratio dy/dx
    dydx <- dy / abs_dx

    # There are 10 or  8 output variables 
    # return partial derivatives as a 8x1 array
    return (dydx)
}
