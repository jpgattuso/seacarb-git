\name{SIR_b}
\alias{SIR_b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameters of the seawater carbonate system with boron addition, including the substrate-inhibitor-ratio (SIR)
}
\description{
This function is a modification of \emph{carbb}. It returns parameters of the seawater carbonate system when boron is added, as well as additional parameters relevant to the SIR.
}
\usage{
SIR_b(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, 
  k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential", badd=0, 
  warn="y", eos = "eos80", long = 1e+20, lat = 1e+20)
}

\arguments{
  \item{flag}{select the couple of variables available. The flags which can be used are: 

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
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010) and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickson (1990) and "k" for using Ks from Khoo et al. (1977), default is "d"} 
	\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
	\item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74" }
  \item{gas}{used to indicate the convention for INPUT pCO2, i.e., when it is an input variable (flags 21 to 25): "insitu" indicates it is referenced to in situ pressure and in situ temperature; "potential" indicates it is referenced to 1 atm pressure and potential temperature; and "standard" indicates it is referenced to 1 atm pressure and in situ temperature. All three options should give identical results at surface pressure. This option is not used when pCO2 is not an input variable (flags 1 to 15). The default is "potential".}
  \item{badd}{Amount of boron added in mol/kg.}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for constants; "n" to supress warnings. The default is "y".}
  \item{eos}{"teos10" to specify T and S according to Thermodynamic Equation Of Seawater - 2010 (TEOS-10); "eos80" to specify T and S according to EOS-80.}
  \item{long}{longitude of data point, used when eos parameter is "teos10" as a conversion parameter from absolute to practical salinity.}
  \item{lat}{latitude of data point, used when eos parameter is "teos10".}

}
\details{Calculates the "substrate-inhibitor ratio (SIR)" (i.e. [HCO3-]/[H+]) which is used as a metric to predict carbonate chemistry dependency of biotic calcium carbonate formation (Bach, 2015).  Please note that hydrogen ion concentrations in the output [HCO3-]/[H+] are on the \emph{free scale}, regardless of the input pH scale. We have included pH and [H+] outputs on multiple scales as a teaching exercise to show how calculating the SIR on different scales changes its meaning.

The function is based on \emph{carbb} and therefore also returns parameters of the seawater carbonate system when boron is added. 

\bold{\emph{carbb} details}:

The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\emph{For K1 and K2:}
\itemize{
\item Roy et al. (1993): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Millero et al. (2006): S ranging between 0.1 and 50 and T ranging between 1 and 50oC.

\item Millero (2010): S ranging between 1 and 50 and T ranging between 0 and 50oC. Millero (2010) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale". 
}

\emph{For Kf:}
\itemize{
\item Perez and Fraga (1987): S ranging between 10 and 40 and T ranging between 9 and 33oC.

\item Dickson and Riley (1979 in Dickson and Goyet, 1994): S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

\emph{For Ks:}
\itemize{
\item Dickson (1990): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Khoo et al. (1977): S ranging between 20 and 45 and T ranging between 5 and 40oC.
}

The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It is recommended to use either vectors with the same dimension or one vector for one argument and numbers for the other arguments.

\emph{Pressure corrections and pH scale:}

\itemize{
\item For K0, the pressure correction term of Weiss (1974) is used.

\item For K1, K2, pK1, pK2, pK3, Kw, Kb, Khs and Ksi, the pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

\item For Kf, the pressure correction was applied on the free scale. The formulation of Dickson and Riley (1979 in Dickson and Goyet, 1994) provides Kf on the free scale but that of Perez and Fraga (1987) provides it on the total scale. Hence, in that case, Kf was first transformed from the total scale to the free scale. With both formulations, the pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Ks, the pressure correction was applied on the free scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Kn, The pressure correction was applied on the seawater scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).
}

long and lat are used as conversion parameters from absolute to practical salinity: when seawater is not of standard composition, practical salinity alone is not sufficient to compute absolute salinity and vice-versa. One needs to know the density. When long and lat are given, density is inferred from WOA silicate concentration at given location. When they are not, an arbitrary geographic point is chosen: mid equatorial Atlantic. Note that this implies an error on computed salinity up to 0.02 g/kg.

}


\value{
The function returns a data frame containing the following columns:
  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm}
  \item{P}{Hydrostatic pressure in bar}
  \item{pH}{pH}
  \item{CO2}{CO2 concentration (mol/kg)}
  \item{pCO2}{"standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{fCO2}{"standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{pCO2pot}{"potential" pCO2, CO2 partial pressure computed at potential temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{fCO2pot}{"potential" fCO2, CO2 fugacity computed at potential temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{pCO2insitu}{"in situ" pCO2, CO2 partial pressure computed at in situ temperature and total pressure (atm + hydrostatic) (\eqn{\mu}atm)}
  \item{fCO2insitu}{"in situ" fCO2, CO2 fugacity computed at in situ temperature and total pressure (atm + hydrostatic) (\eqn{\mu}atm)}
  \item{HCO3}{HCO3 concentration (mol/kg)}
  \item{CO3}{CO3 concentration (mol/kg)}
  \item{DIC}{DIC concentration (mol/kg)}
  \item{ALK}{ALK, total alkalinity (mol/kg)}
  \item{OmegaAragonite}{Omega aragonite, aragonite saturation state}
  \item{OmegaCalcite}{Omega calcite, calcite saturation state}
  \item{SIR}{substrate inhibitor ratio: [HCO3]/[HCO3]/[H_free]  (mol:\eqn{\mu}mol)}
  \item{H_free}{H ion concentration on free scale (\eqn{\mu}mol/kg)}
  \item{H_sws}{H ion concentration on seawater scale (\eqn{\mu}mol/kg)}
  \item{H_t}{H ion concentration on total scale (\eqn{\mu}mol/kg)}
  \item{pH_free}{pH calculated on the free scale}
  \item{pH_sws}{pH calculated on the seawater scale}
  \item{pH_t}{pH calculated on the total scale}

}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2015)}

\references{

Bach L. T., 2015 Reconsidering the role of carbonate ion concentration in calcification by marine organisms. \emph{Biogeosciences} \bold{12}, 4939-4951.


Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen ion concentration in seawater from 5 to 40oC: standard potentials at salinities from 20 to 45. \emph{Analytical Chemistry} \bold{49}, 29-34. 

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.

Millero F. J., 2010. Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006. Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2015. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences} \bold{12}, 1483-1510.

Perez F. F. and Fraga F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of
a non-ideal gas, \emph{Mar.  Chem.}, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and
seawater, \emph{Marine Chemistry}, \bold{8}, 347-359.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.

}
\author{
Heloise Lavigne, James Orr, Jean-Pierre Gattuso, Lennart Bach and Kimberlee Baldry \email{gattuso@obs-vlfr.fr}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## With a couple of variables
SIR_b(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", badd=0)

}
\keyword{utilities}

