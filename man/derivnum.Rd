\encoding{latin1}
\name{derivnum}
\alias{derivnum}
\title{Numerical derivatives of seawater carbonate system variables}
\description{Returns numerical derivatives of the seawater carbonate system output variables with respect to input variables.}
\usage{derivnum(varid, flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0,
       k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential", warn="y", 
       eos = "eos80", long = 1e+20, lat = 1e+20)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{varid}{Variable length, case insensitive, character identificator of variable with respect to which derivatives are requested.
Possible values are:

'1' or 'var1'  :  Variable 1 of the input pair (This is TAlk if flag is 15)

'2' or 'var2'  :  Variable 2 of the input pair (This is DIC  if flag is 15)

'sil', 'silt', 'tsil' or 'silicate'      : Total silicate concentration

'phos', 'phost', 'tphos' or 'phosphate'  : Total phosphate concentration

't', 'temp' or 'temperature' : temperature

's', 'sal' or 'salinity'     : salinity

'K0','K1','K2','Kb','Kw','Kspa' or 'Kspc' : one of the dissociation constants

'bor' : total boron
}
\item{flag}{select the input pair of carbonate-system variables
  available by choosing the following flag:

flag = 1      pH and CO2 given

flag = 2      CO2 and HCO3 given

flag = 3      CO2 and CO3 given

flag = 4      CO2 and ALK given

flag = 5      CO2 and DIC given

flag = 6      pH and HCO3 given

flag = 7      pH and CO3 given

flag = 8      pH and ALK given

flag = 9      pH and DIC given

flag = 10     HCO3 and CO3 given

flag = 11     HCO3 and ALK given

flag = 12     HCO3 and DIC given

flag = 13     CO3 and ALK given

flag = 14     CO3 and DIC given

flag = 15     ALK and DIC given

flag = 21     pCO2 and pH given

flag = 22     pCO2 and HCO3 given

flag = 23     pCO2 and CO3 given

flag = 24     pCO2 and ALK given

flag = 25     pCO2 and DIC given
}
\item{var1}{Value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
\item{var2}{Value of the second  variable in mol/kg, except for pH}
\item{S}{Salinity}
\item{T}{Temperature in degrees Celsius}
 \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
 \item{P}{Hydrostatic pressure in bar (surface = 0)}
 \item{Pt}{Concentration of total phosphate in mol/kg; set to 0 if NA}
 \item{Sit}{Concentration of total silicate in mol/kg; set to 0 if NA}
 \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010), "w14" from Waters et al. (2014), and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
\item{ks}{"d" for using Ks from Dickson (1990) and "k" for using Ks from Khoo et al. (1977), default is "d"}
\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
\item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74" }
 \item{gas}{used to indicate the convention for INPUT pCO2, i.e., when it is an input variable (flags 21 to 25): "insitu" indicates it is referenced to in situ pressure and in situ temperature; "potential" indicates it is referenced to 1 atm pressure and potential temperature; and "standard" indicates it is referenced to 1 atm pressure and in situ temperature. All three options should give identical results at surface pressure. This option is not used when pCO2 is not an input variable (flags 1 to 15). The default is "potential".}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for constants; "n" to supress warnings. The default is "y".}
    \item{eos}{"teos10" to specify T and S according to Thermodynamic Equation Of Seawater - 2010 (TEOS-10); "eos80" to specify T and S according to EOS-80.}
  \item{long}{longitude of data point, used when eos parameter is "teos10" as a conversion parameter from absolute to practical salinity.}
  \item{lat}{latitude of data point, used when eos parameter is "teos10".}
}

\details{
This subroutine has same input parameters as subroutine carb(). For details on these parameters,
refer to documentation of 'carb'.

This subroutine computes partial derivatives of each output variable with respect to
each of the input variable (including each of the two chosen carbonate system
variables, each of the nutrients (total silicon and total phosphorus), temperature,
and salinity.

It computes these derivatives (dy/dx) using the method of central differences, i.e.,
\itemize{
  \item for dx, it adds a positive and negative perturbation, same and equal in magnitude, to each input variable, one at a time, and
  \item for dy, it then computes the corresponding induced change in output variables
}

All arguments but the first (varid), can be given as scalers or vectors. If the
lengths of the vectors differs, only the longest vector is retained and the other
arguments are set equal to the first value of the other vectors are used. Hence users
should use either vectors with the same dimension or one vector for one argument and
scalars for others; otherwise, results may not be as intended.

long and lat are used as conversion parameters from absolute to practical salinity: when seawater is not of standard composition, practical salinity alone is not sufficient to compute absolute salinity and vice-versa. One needs to know the density. When long and lat are given, density is inferred from WOA silicate concentration at given location. When they are not, an arbitrary geographic point is chosen: mid equatorial Atlantic. Note that this implies an error on computed salinity up to 0.02 g/kg.
}

\value{The function returns a data frame containing the following columns:
 \item{H}{derivative of [H+] concentration  (mol/kg/...)}
 \item{pH}{derivative of pH}
 \item{CO2}{derivative of CO2 concentration (mol/kg/...)}
 \item{pCO2}{"standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm/...)}
 \item{fCO2}{"standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm/...)}
 \item{HCO3}{derivative of HCO3 concentration (mol/kg/...)}
 \item{CO3}{derivative of CO3 concentration (mol/kg/...)}
 \item{DIC}{derivative of DIC concentration (mol/kg/...)}
 \item{ALK}{derivative of ALK, total alkalinity (mol/kg/...)}
 \item{OmegaAragonite}{derivative of Omega aragonite, aragonite saturation state}
 \item{OmegaCalcite}{derivative of Omega calcite, calcite saturation state}

If all input data have the same 'flag' value, returned data frame does not show
derivatives of input pair of carbonate system variables. For example, if all input
flags are 15, the input pair is DIC and ALK; hence, derivatives of DIC and ALK are
not returned.

Units of derivative dy/dx is unit(y)/unit(x) where unit(x) are as follows:
 \item{degree C}{when x is Temperature}
 \item{psu}{when x is Salinity}
 \item{\eqn{\mu}atm}{when x is pCO2}
 \item{mol/kg}{for all other cases}
}

\references{
Dickson A. G. and Riley J. P., 1979 The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water. \emph{Marine Chemistry} \bold{7}, 89-99.

Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen Ion Concentration in Seawater from 5 to 40oC: Standard Potentials at Salinities from 20 to 45. \emph{Analytical Chemistry} \bold{22}, vol49 29-34.

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.

Millero F. J., 2010. Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006. Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2015. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences} \bold{12}, 1483-1510.

Perez F. F. and Fraga F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.

Waters, J., Millero, F. J., and Woosley, R. J., 2014. Corrigendum to ``The free proton concentration scale for seawater pH'', [MARCHE: 149 (2013) 8-22], \emph{Marine  Chemistry}, \bold{165}, 66-67.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of a non-ideal gas, \emph{Marine  Chemistry}, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and seawater, \emph{Marine  Chemistry}, \bold{8}, 347-359.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}

\author{
Jean-Marie Epitalon, James Orr, and Jean-Pierre Gattuso\email{gattuso@obs-vlfr.fr}
}

\examples{

## 1) For the input pair ALK and DIC (var1 and var2 when flag=15)
##    compute derivatives of all output varialbes 
##    with respect to DIC (i.e., var2)
derivnum(varid='var2', flag=15, var1=2300e-6, var2=2000e-6, 
         S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
         pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 2) For the input pair pH and ALK (var1 and var2 when flag=8)
##    compute derivatives of all output variables 
##    with respect to [H+] concentration
derivnum(varid='var1', flag=8, var1=8.2, var2=0.00234, 
         S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
         pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 3) Using vectors as arguments and compute derivatives of all output 
##    variables with respect to temperature
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("u74", "u74", "u74")
derivnum(varid='T', flag=flag, var1=var1, var2=var2, S=S, T=T, P=P,
         Pt=Pt, Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

# For more examples of use of derivnum.R,
# consult the code of seacarb's errors routine.
}

\keyword{utilities}
