\encoding{latin1}
\name{Om}
\alias{Om}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Carbonate saturation state for magnesian calcites}
\description{Calculates the calcium carbonate saturation state for magnesian calcite}
\usage{Om(x, flag, var1, var2, k1k2='x', kf='x', ks="d", pHscale="T", b="l10")}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{mole fraction of magnesium ions, note that the function is only valid for x ranging between 0 and 0.25}
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
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010) and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m10".}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickon (1990) and "k" for using Ks from Khoo et al. (1977), default is "d"} 
	\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
	\item{b}{"l10" for computing boron total from the Lee et al. (2010) formulation or "u74" for using the Uppstrom (1974) formulation, default is "l10" }
}


\details{It is important to note that this function is \bold{only valid} for:
\itemize{
	\item Salinity = 35	
	\item Temperature = 25 degrees Celsius	
  	\item Hydrostatic pressure = 0 bar (surface)  	
  	\item Concentration of total phosphate = 0 mol/kg
    	\item Concentration of total silicate = 0 mol/kg
}


Note that the stoichiometric solubility products with respect to Mg-calcite minerals have not been determined experimentally. The saturation state with respect to Mg-calcite minerals is therefore calculated based on ion activities, i.e.,

\deqn{
\Omega_{x} = \frac{ \{Ca^{2+}\}^{1-x} \{Mg^{2+}\}^{x} \{CO_{3}\}^{2-} } { K_{x} }
}

The ion activity \{a\} is calculated based on the observed ion concentrations [C] multiplied by the total ion activity coefficient, \eqn{\gamma_T}, which has been determined experimentally or from theory (e.g. Millero & Pierrot 1998): \{a\}=\eqn{\gamma_T}[C]. Because a true equilibrium cannot be achieved with respect to Mg-calcite minerals, \eqn{K_x} represents a metastable equilibrium state obtained from what has been referred to as stoichiometric saturation (Thorstenson & Plummer 1977; a term not equivalent to the definition of the stoichiometric solubility product, see for example Morse et al. (2006) and references therein). In the present calculation calcium and magnesium concentrations were calculated based on salinity. Total ion activity coefficients with respect to \eqn{Ca^{2+}}, \eqn{Mg^{2+}}, and \eqn{CO_{3}^{2-}} were adopted from Millero & Pierrot (1998).  


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
}


\value{The function returns a list with
  \item{OmegaMgCa_biogenic}{Mg-calcite saturation state for minimally prepared biogenic Mg-calcite.}
  \item{OmegaMgCa_biogenic_cleaned}{Mg-calcite saturation state for cleaned and annealed biogenic Mg-calcite.}
}

\references{
Only the references related to the saturation state of magnesian calcite are listed below; the other references are listed under the carb function.

Andersson A. J., Mackenzie F. T., Nicholas R. B., 2008, Life on the margin: implications of ocean acidification on Mg-calcite, high latitude and cold-water marine calcifiers. \emph{Marine Ecology Progress Series} \bold{373}, 265-273.

Bischoff W. D., Bertram M. A., Mackenzie F. T. and Bishop F.C., 1993 Diagenetic stabilization pathways of magnesian calcites. \emph{Carbonates and Evaporites} \bold{8}, 82-89.

Millero F. J. and Pierrot D., 1998. A chemical equilibrium model for natural waters. \emph{Aquatic Geochemistry} \bold{4}, 153-199.

Morse J. W., Andersson A. J. and Mackenzie F. T., 2006. Initial responses of carbonate-rich shelf sediments to rising atmospheric pCO2 and 
ocean acidification: Role of high Mg-calcites. \emph{Geochimica et Cosmochimica Acta} \bold{70}, 5814-5830.

Thorstenson D.C. and Plummer L.N., 1977. Equilibrium criteria for two component solids reacting with fixed composition in an aqueous phase-example: the magnesian calcites. \emph{American Journal of Science} \bold{277}, 1203-1233.
}

\author{
Heloise Lavigne, Andreas J. Andersson and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
Om(x=seq(0.01, 0.252, 0.01), flag=8, var1=8.2, var2=0.00234, 
  k1k2='x', kf='x', ks="d", pHscale="T", b="l10")
}
\keyword{utilities}
