\encoding{latin1}
\name{oa}
\alias{oa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perturbation of the seawater carbonate system}
\description{Describes the various approaches that can be used to alter the seawater carbonate system. Its main purpose is to assist the design of ocean acidification perturbation experiments.}
\usage{oa(flag, var1, var2, pCO2f, pCO2s=1e6, S=35, T=25, P=0, 
	Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", plot=FALSE, 
	b="u74", eos = "eos80", long = 1e+20, lat = 1e+20)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{flag}{select the couple of variables available to describe the initial seawater. The flags which can be used are: 

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
	\item{var1}{Value of the first variable available to describe the initial seawater, in mol/kg except for pH and for pCO2 in uatm}
	\item{var2}{Value of the second variable available to describe the initial seawater, in mol/kg except for pH}
	\item{pCO2f}{pCO2 target value, in uatm}
	\item{pCO2s}{pCO2s is the pCO2, in uatm, of the \dQuote{high-CO2} seawater that will be mixed with \dQuote{normal seawater}. The default value is 10^6 uatm, that is seawater bubbled with pure CO2 gas and saturated with CO2.}
	\item{S}{Salinity, default is 35}
	\item{T}{Temperature in degrees Celsius, default is 25}
  	\item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  	\item{Pt}{Concentration of total phosphate in mol/kg, default is 0}
  	\item{Sit}{Concentration of total silicate in mol/kg, default is 0}
        \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010), "w14" from Waters et al. (2014), and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickon (1990), "k" for using Ks from Khoo et al. (1977), default is "d"} 
	\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
	\item{plot}{A plot of the different perturbation methods can be plotted in a DIC vs ALK field with pCO2 isoclines are drawn in the back. Default is false.}
	\item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74".}
	 \item{eos}{"teos10" to specify T and S according to Thermodynamic Equation Of Seawater - 2010 (TEOS-10); "eos80" to specify T and S according to EOS-80.}
  \item{long}{longitude of data point, used when eos parameter is "teos10" as a conversion parameter from absolute to practical salinity.}
  \item{lat}{latitude of data point, used when eos parameter is "teos10".}
}

\details{The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

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
}


\section{Warnings}{

\itemize{
\item It is recommended to use concentrated solutions of acid or base in order to add small volumes.

\item The addition of strong acid does not simulate well natural ocean acidification (higher concentration of dissolved inorganic carbon at constant total alkalinity) since it generates a decrease in total alkalinity while dissolved inorganic carbon is kept constant.

\item Other important advice is provided in Gattuso and Lavigne (2009), Schulz et al. (2009) and in the \dQuote{Guide for Best Practices on Ocean Acidification Research and Data Reporting} (\url{https://www.iaea.org/ocean-acidification/page.php?page=2194})
}
long and lat are used as conversion parameters from absolute to practical salinity: when seawater is not of standard composition, practical salinity alone is not sufficient to compute absolute salinity and vice-versa. One needs to know the density. When long and lat are given, density is inferred from WOA silicate concentration at given location. When they are not, an arbitrary geographic point is chosen: mid equatorial Atlantic. Note that this implies an error on computed salinity up to 0.02 g/kg.
}


\value{The function returns a list built as follows: 

\item{description}{A table describing in plain English the various ways to reach the target pCO2. Note that if a vector is given in argument only the first value is used.}

\item{perturbation}{Table providing key parameters for the following methods:

\emph{CO2 bubbling}: high-CO2 air is bubbled in seawater. The first parameter is the value of the pCO2 in the air required to bubble the seawater (in uatm).

\emph{SW mixing}: mixing of \dQuote{normal} and \dQuote{high-CO2} seawater. The fist parameter, \dQuote{Weight fraction high-CO2 SW} or wf, is the weight fraction of the high-CO2 seawater per kg seawater.

\emph{Addition of acid}: strong acid is added to seawater. Note that this method is not recommended because it does not closely mimic natural ocean acidification (Gattuso and Lavigne, 2009). The first parameter, H+ (mol/kg), is the amount of H+ that must be added (mol/kg). The acid must be fortified with NaCl in order to have the same salinity than seawater.

\emph{Addition of HCO3 and acid}: bicarbonate (HCO3) and a strong acid are added. The first parameter, HCO3, is the amount of HCO3 that must be added (mol/kg). The second parameter, H+, is the quantity of H+ that must be added (mol/kg). The acid must be fortified with NaCl in order to have the same salinity than seawater.

\emph{Addition of CO3 and acid}: carbonate, CO3, and a strong acid are added. The first parameter, HCO3, is the quantity of CO3 that must be added (mol/kg). The second parameter, H+, is the quantity of H+ that must be added (mol/kg).
}

\item{summary}{Table summarizing the carbonate chemistry before and after applying each perturbation: pCO2 bubbling, mixing with high-CO2 seawater, addition of strong acid, and addition of bicarbonate/carbonate and strong acid.}
}

\references{
Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen Ion Concentration in Seawater from 5 to 40oC: Standard Potentials at Salinities from 20 to 45. \emph{Analytical Chemistry} \bold{49}, 29-34. 

Gattuso J.-P. and Lavigne H., 2009 Technical note: approaches and software tools to investigate the impact of ocean acidification. \emph{Biogeosciences} \bold{21}, 6:2121-2133. 

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Millero F. J., 2010 Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006. Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Perez F. F. and Fraga F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993 The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Schulz K. G., Barcelos e Ramos J., Zeebe R. E. and Riebesell U., 2009 CO2 perturbation experiments: similarities and differences between dissolved inorganic carbon and total alkalinity manipulations. \emph{Biogeosciences} \bold{6}, 2145-2153.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.

Waters, J., Millero, F. J., and Woosley, R. J., 2014. Corrigendum to ``The free proton concentration scale for seawater pH'', [MARCHE: 149 (2013) 8-22], Marine Chemistry 165, 66-67.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}

\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{carb}}, \code{\link{pgas}}, \code{\link{pmix}}, \code{\link{ppH}}, \code{\link{pTA}}.
}

\examples{
oa(flag=24, var1=384, var2=2325e-6, pCO2s=1e6, pCO2f=793, S=34.3, T=16, 
	P=0, pHscale="T", kf="pf", k1k2="l", ks="d", plot=TRUE, b="u74")
}

\keyword{utilities}
