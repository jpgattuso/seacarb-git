\encoding{latin1}
\name{psi}
\alias{psi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Molar ratio of CO2 released vs CaCO3 precipitated}
\description{Returns the molar ratio of CO2 released vs CaCO3 precipitated described by Frankignoulle et al. (1994).}
\usage{
psi(flag, var1, var2, S=35, T=20, Patm, P=0, Pt=0, Sit=0, pHscale="T", 
	kf="x", k1k2="x", ks="d", eos = "eos80", long = 1e+20, lat = 1e+20)
}
%- maybe also 'usage' for other objects documented here.
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
	\item{var1}{enter value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{enter value of the second  variable in mol/kg, except for pH}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
  \item{P}{Hydrostatic pressure in bar (surface = 0)}
  \item{Pt}{Concentration of total phosphate in mol/kg; set to 0 if NA}
  \item{Sit}{Concentration of total silicate in mol/kg; set to 0 if NA}
	\item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
  \item{k1k2}{"cw" for using K1 and K2 from Cai & Wang (1998), "l" from Lueker et al. (2000), "m02" from Millero et al. (2002), "m06" from Millero et al. (2006), "m10" from Millero (2010), "mp2" from Mojica Prieto et al. (2002), "p18" from Papadimitriou et al. (2018), "r" from Roy et al. (1993), "sb21" from Shockman & Byrne (2021), "s20" from Sulpis et al. (2020), and "w14" from Waters et al. (2014). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
  \item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickon (1990), "k" for using Ks from Khoo et al. (1977), default is "d"}
	\item{eos}{"teos10" to specify T and S according to Thermodynamic Equation Of Seawater - 2010 (TEOS-10); "eos80" to specify T and S according to EOS-80.}
  \item{long}{longitude of data point, used when eos parameter is "teos10" as a conversion parameter from absolute to practical salinity.}
  \item{lat}{latitude of data point, used when eos parameter is "teos10".}
}

\details{The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\emph{For K1 and K2:}
\itemize{
\item Cai and Wang (1998): S ranging between 0 and 40 and T ranging between 0.2 and 30oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Millero et al. (2002): S ranging from 34 to 37 and T ranging between -1.6 and 35oC.

\item Millero et al. (2006): S ranging between 0.1 and 50 and T ranging between 1 and 50oC.

\item Millero (2010): S ranging between 1 and 50 and T ranging between 0 and 50oC. Millero (2010) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale". 

\item Mojica Prieto et al. (2002): S ranging from 5 to 42 and T ranging between 0 and 45oC.

\item Papadimitriou et al. (2018): S ranging from 33 to 100 and T ranging between -6 to 25oC.

\item Roy et al. (1993): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Shockman & Byrne (2021): for K2, S ranging from 19.6 to 41 and T ranging between 15 to 35oC. For K1, formulation is that of Waters et al.

\item Sulpis et al. (2020): S ranging from 30.7 to 37.6 and T ranging between -1.7 to 31.8oC.

\item Waters et al.(2014): S ranging between 1 and 50 and T ranging between 0 and 50oC. Waters (2014) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale".
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

\value{The function returns a data frame containing the following columns:
	\item{psi}{ratio of CO2 released vs CaCO3 precipitated (mol/mol)}
}

\references{
Cai W. J., and Wang Y., 1998. The chemistry, fluxes, and sources of carbon dioxide in the estuarine waters of the Satilla and Altamaha Rivers, Georgia. \emph{Limnology and Oceanography} \bold{43}, 657-668.

Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen Ion Concentration in Seawater from 5 to 40oC: Standard Potentials at Salinities from 20 to 45. \emph{Analytical Chemistry} \bold{49}, 29-34. 

Frankignoulle M., 1994 A complete set of buffer factors for acid/base CO2 system in seawater. \emph{Journal of Marine Systems} \bold{5}, 111-118.

Frankignoulle M., Canon C. and Gattuso J.P., 1994 Marine calcification as a source of carbon dioxide- Positive feedback of increasing atmospheric CO2. \emph{Limnology and Oceanography} \bold{2}, 458-462.

Millero F. J., 2010 Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006 Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Waters, J., Millero, F. J., and Woosley, R. J., 2014. Corrigendum to ``The free proton concentration scale for seawater pH'', [MARCHE: 149 (2013) 8-22], Marine Chemistry 165, 66-67.
}

\author{
Jean-Pierre Gattuso and Heloise Lavigne \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link[seacarb]{speciation}}.
}

\examples{
## Calculation using the numerical example given in Frankignoulle et al. (1994)
psi(flag=24, var1=350, var2=2400e-6, S=35, T=25, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d")
}
\keyword{utilities}

