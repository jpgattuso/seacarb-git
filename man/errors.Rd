\encoding{latin1}
\name{errors}
\alias{errors}
\title{Error propagation for computed marine carbonate system variables}
\description{Estimates uncertainties in computed carbonate system variables by propagating standard error (uncertainty) in six input variables, including:
\itemize{
\item the input pair of carbonate system variables, 
\item the 2 input nutrients (silicate and phosphate concentrations),
\item temperature and salinity. It also accounts for
\item the errors in the key dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
\item the error in total boron
}
}
\usage{errors(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, 
              evar1=0, evar2=0, eS=0.01, eT=0.01, ePt=0, eSit=0, 
              epK=c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02),
              eBt=0.01, method = "ga", r=0.0, runs=10000,
              k1k2='x', kf='x', ks="d", pHscale="T", b="u74", gas="potential", 
              warn="y", eos = "eos80", long = 1e+20, lat = 1e+20)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{flag}{select the pair of carbonate system input variables. The flags to be used are as follows:

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
\item{var1}{Value of the first variable (in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm)}
\item{var2}{Value of the second variable (in mol/kg, except for pH)}
\item{S}{Salinity (practical salinity scale)}
\item{T}{Temperature in degrees Celsius}
\item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
\item{P}{Hydrostatic pressure in bar (surface = 0)}
\item{Pt}{Concentration of total dissolved inorganic phosphorus (mol/kg); set to 0 if NA}
\item{Sit}{Concentration of total dissolved inorganic silicon (mol/kg); set to 0 if NA}
\item{evar1}{Error (uncertainty) in var1 of input pair of carbonate system variables}
\item{evar2}{Error (uncertainty) in var2 of input pair of carbonate system variables}
\item{eS}{Error (uncertainty) in salinity; default is 0.01}
\item{eT}{Error (uncertainty) in temperature (degree C); default is 0.01}
\item{ePt}{Error (uncertainty) in total dissolved inorganic phosphorus concentration (mol/kg)}
\item{eSit}{Error (uncertainty) in total dissolved inorganic silicon concentration (mol/kg)}
\item{epK}{Error (uncertainty) in 7 key dissociation constants: pK0, pK1, pK2, pKb, pKw, pKspa and pKspc. This is a vector. The default is c(0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02).}
\item{eBt}{Error (uncertainty) in total boron, given as a
  relative fractional error. The default is 0.01, which equates to a 1\% error}
\item{method}{Case insensitive character string to choose the error-propagation method: 1) Gaussian, 
2) Method of Moments, or 3) Monte Carlo).\cr
These methods are specified using the 2-letter codes "ga", "mo", or "mc", respectively. The default is "ga" (Gaussian).}
\item{r}{Correlation coefficient between uncertainties of var1 and var2 (only useful with method="mo",
i.e., ignored for the 2 other methods, the default is r=0.0}
\item{runs}{Number of random samples (ignored unless method="mc"; the default is runs=10000}
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
This function requires users to specify each input uncertainty as either
the standard deviation or the standard error of the mean. The latter
implies much smaller propagated uncertainties, but is appropriate only
when interested in the error in the mean, not the error of a given
measurement. Beware that it is easy to fool oneself when using the
standard error of the mean rather than the standard deviation.

This function requires different types of uncertainties:
\itemize{
\item absolute uncertainties for evar1, evar2, eS, eT, ePt, eSit (same units as the input data, e.g., mol/kg);
\item uncertainties in pK units for epK; and
\item uncertainties in relative fractional units (between 0.0 and 1.0) for eBt.
}

This function propagates uncertainties from input to output variables using one of three methods: 
\itemize{
\item Gaussian:
     The Gaussian method is the standard technique for estimating a
     computed variable's (z) second moment (its variance or standard
     deviation) based on a first-order approximation to z.  More
     precisely, we use here the basic 1st order, 2nd moment
     uncertainty analysis (a type of Taylor expansion), assuming no
     covariance between input variables.  This is the approach used by
     Dickson and Riley (1978). It is the default method.

\item Method of moments: The method of moments is a more general form
     of the Gaussian method. But in addition, it also accounts for
     covariance between input variables.  In this case, the 'errors'
     routine allows the user to specify a value of the correlation
     coefficient 'r', having a value between -1.0 and 1.0, to indicate
     the correlation between uncertainties of the input pair of
     carbonate system variables. That correlation is used to comput
     the covariance. But by default, it is assumed that there is no
     covariance (r=0.0).

\item Monte Carlo:
     The Monte Carlo method is a brute-force approach relying on
     repeated random sampling of input errors, adding those to each
     input variables, calculating the corresponding output variables
     for each sample, and finally assessing the standard deviation in
     each output variables.
}

This function has many input parameters that are identical to those
in the carb function. For their details, refer to the 'carb' documentation.

All parameters may be scalars or vectors except epK, eBt, method, runs, and gas.

\itemize{
  \item runs and eBt must be scalars
  \item method and gas must each consist of a character string
  \item epK may be a vector of 7 values. In that case, it must list errors for 
    pK0, pK1, pK2, pKb, pKw, pKspa and pKspc, respectively.
    That set of errors is identical for all input data.
    Alternatively, users may specify 'epK=NULL' or 'epK=0' to set all 7 values to zero
    and thus neglect errors in the equilibrium constants. 
}
In constrast, for evar1, evar2, r, eS, eT, ePt and eSit:

\itemize{
  \item  if they are vectors, they represent errors associated with each data point
  \item  if they are scalars (single real numbers), they represent one error value each associated to all data points
}

The same remark applies to parameter r (correlation coefficient).

long and lat are used as conversion parameters from absolute to practical salinity: when seawater is not of standard composition, practical salinity alone is not sufficient to compute absolute salinity and vice-versa. One needs to know the density. When long and lat are given, density is inferred from WOA silicate concentration at given location. When they are not, an arbitrary geographic point is chosen: mid equatorial Atlantic. Note that this implies an error on computed salinity up to 0.02 g/kg.
}

\section{Correlation coefficient}{
By default, 'r' is zero. However, for some pairs the user may want to specify a different value. For example, measurements of pCO2 and pH are often anti-correlated. The same goes for two other pairs: 'CO2 and CO3' and 'pCO2 and CO3'.  But even for these cases, care is needed before using non-zero values of 'r'.

When the user wishes to propagate errors for an individual measurement, 'r' should ALWAYS be zero if each member of the input pair is measured independently. In this case, we are interested in the correlation between the uncertainties in those measurements, not in the correlation between the measurements themselves. Uncertainties from those measurements are probably not correlated if they come from different instruments. Conversely, if users are interested in the error in the mean of a distribution of measurements (i.e., if they are propagating standard errors instead of standard deviations), one should then also account for the correlation between the measurements of the two variables of the input pair.

For input pairs where one member is pH (flags 1, 6, 7, 8, 9, and 21), this 'errors' function automatically inverses the sign of 'r'. The reason for that is that the associated derivatives are computed in terms of the hydrogen ion concentration (H+), not pH. Therefore for each of these 6 flags, if the user wants to compute their own 'r' that should be done by (1) using the H+ concentration instead of pH, and (2) inversing the sign of that computed 'r' before passing it as an argument to this routine. Usually though (when not calculating r for pH), the user may just use the 'r' in the expected way.  For example, to include the covariance term when there is a perfect anticorrelation of pH with pCO2, one would use 'r=-1.0'.

}

\section{Computation time}{
Computation time depends on the method chosen; the Monte Carlo method takes much longer to execute. The computational time required for the Monte Carlo method is proportional to the number of runs. More runs, implies improved accuracy: runs = 10000 appears a minimum to obtain an accuracy of less than 1\%. Accuracy is inversely proportional to the number of runs.

Computation time also depends on the chosen pair of input variables. For example, with the input pair DIC and Total alkalinity (flag=15), it is much longer than for input pair pH and Total alkalinity (flag=8)
}

\value{The function returns a 2-dimensional dataframe, with the following columns:
 \itemize{
 \item{H}{total error in [H+] concentration (mol/kg)}
 \item{pH}{total error in pH}
 \item{CO2}{total error in CO2 concentration (mol/kg)}
 \item{pCO2}{total error in "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
 \item{fCO2}{total error in "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
 \item{HCO3}{total error in HCO3 concentration (mol/kg)}
 \item{CO3}{total error in CO3 concentration (mol/kg)}
 \item{DIC}{total error in DIC concentration (mol/kg)}
 \item{ALK}{total error in ALK, total alkalinity (mol/kg)}
 \item{OmegaAragonite}{total error in Omega aragonite (aragonite saturation state)}
 \item{OmegaCalcite}{total error in Omega calcite   (calcite saturation state)}
}
 
If all input data have the same 'flag' value, the returned data frame does not show errors on input pair of carbonate system variables. For example, if all input flags are 15, the input pair is DIC and ALK; hence, errors on DIC and ALK are not returned.
}

\references{
Dickson, A. G. and Riley, J. P., 1978 The effect of analytical error on the evaluation of the components of the aquatic carbon-dioxide system, \emph{Marine  Chemistry}, \bold{6}, 77-85.

Dickson A. G. and Riley J. P., 1979 The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water. \emph{Marine Chemistry} \bold{7}, 89-99.

Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of hydrogen ion concentration in seawater from 5 to 40oC: standard potentials at salinities from 20 to 45. \emph{Analytical Chemistry} \bold{22}, vol49 29-34.

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

## 1) For the input pair ALK and DIC (var1 and var2 when flag=15),
## compute resulting uncertainty from given uncertainty on ALK and DIC (5 umol/kg)
## and default uncertainties in dissociation constants and total boron
## using the default method (Gaussian)
errors(flag=15, var1=2300e-6, var2=2000e-6, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=5e-6, evar2=5e-6, eS=0, eT=0, ePt=0, eSit=0, 
       pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")
## Typical output:
## H             pH          CO2           fCO2      pCO2      HCO3          ...
## 3.721614e-10  0.01796767  5.441869e-07  19.25338  19.31504  9.170116e-06  ...

## 2) Do the same as in one, but assign a 4% uncertainty to total boron
##    This uncertainty is the amount by which estimates from Lee et al (2010) and 
##    Uppstrom (1974) differ. The default for the latter is eBt=0.01.
errors(flag=15, var1=2300e-6, var2=2000e-6, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=5e-6, evar2=5e-6, eS=0, eT=0, ePt=0, eSit=0, eBt=0.04,
       pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 3) For the input pair pH and ALK (var1 and var2 when flag=8)
## compute standard errors in output variables from errors in input variables, i.e., 
## for pH (0.005 pH units) and in ALK (5 umol/kg), along with
## errors in total dissolved inorganic phosphorus (0.1 umol/kg) and
## total dissolved inorganic silicon (2 umol/kg) concentrations, while
## assuming no uncertainty in dissociation constants & boron, using the Gaussian method:
errors(flag=8, var1=8.25, var2=2300e-6,  S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=0.005, evar2=5e-6, eS=0, eT=0, ePt=0.1, eSit=2, epK=0, eBt=0,
       method="ga", pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 4) For the input pair pCO2 and pH (var1 and var2 when flag=21)
## compute standard errors in output variables from errors in input variables, i.e., 
## for pCO2 (2 uatm) and pH (0.005 pH units), with no uncertainties in Pt and Sit
## nor in the dissociation constants BUT a perfect anticorrelation between pCO2 and pH,
## (the input pair) using the Method of moments:
errors(flag=21, var1=400, var2=8.1,  S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0, 
       evar1=2, evar2=0.005, eS=0, eT=0, ePt=0.0, eSit=0, epK=0, eBt=0, 
       method="mo", r=-1.0, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")

## 5) Use vectors as arguments and compute errors on all output variables
## using Monte Carlo method taking into account input errors on pH, ALK, DIC
## and dissociation constants (pKx)
flag <- c(8, 15, 8)
var1 <- c(8.2, 0.002394, 8.25)
var2 <- c(0.002343955, 0.002017, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 22)
P <- 0
Pt <- 0
Sit <- 0
evar1 <- c(0.005, 2e-6, 0.005)
evar2 <- c(2e-6, 2e-6, 2e-6)
epKx <- c(0.002, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01)
eBtx = 0.01
method <- "mc"
kf <- "pf"
k1k2 <- "l"
pHscale <- "T"
b <- "u74"
## NOTE that the following is executable but enclosed in "donttest" 
## because it takes too long to run when submiting to CRAN
## and is therefore rejected
\donttest{errors(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
       evar1=evar1, evar2=evar2, eS=0, eT=0, ePt=0, eSit=0, epK=epKx, eBt=eBtx,
       method=method, runs=10000, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)}
}

\keyword{utilities}
