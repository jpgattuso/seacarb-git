\encoding{latin1}
\name{K2}
\alias{K2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Second dissociation constant of carbonic acid (mol/kg)}
\description{Second dissociation constant of carbonic acid (mol/kg)}

\usage{
K2(S=35, T=25, P=0, k1k2="x", pHscale="T", kSWS2scale="x", ktotal2SWS_P0="x", warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m02" from Millero et al. (2002), "m06" from Millero et al. (2006), "m10" from Millero (2010), "mp2" from Mojica Prieto et al. (2002), "p18" from Papadimitriou et al. (2018), "r" from Roy et al. (1993), "s20" from Sulpis et al. (2020), "sb21" from Shockman & Byrne (2021), and "w14" from Waters et al. (2014),. "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. It is not required for all formulations of K1 and K2, nor when the hydrostatic pressure is 0. It is advised to use default value "x", in which case it is computed when required.}
  \item{ktotal2SWS_P0}{Conversion factor from the total scale to the SWS at an hydrostatic pressure of 0. It is not required for all formulations of K1 and K2. It is advised to use default value "x", in which case it is computed when required.}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for K2; "n" to supress warnings. The default is "y".}
}

\details{The Lueker et al. (2000) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\itemize{
\item Roy et al. (1993): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Mojica Prieto et al. (2002): S ranging from 5 to 42 and T ranging between 0 and 45oC.

\item Millero et al. (2002): S ranging from 34 to 37 and T ranging between -1.6 and 35oC.

\item Millero et al. (2006): S ranging between 0.1 and 50 and T ranging between 1 and 50oC.

\item Millero (2010): S ranging between 1 and 50 and T ranging between 0 and 50oC. Millero (2010) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale". 
 
\item Waters et al.(2014): S ranging between 1 and 50 and T ranging between 0 and 50oC. Waters (2014) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale". 

\item Papadimitriou et al. (2018): S ranging from 33 to 100 and T ranging between -6 to 25oC.

\item Sulpis et al. (2020): S ranging from 30.7 to 37.6 and T ranging between -1.7 to 31.8oC.

\item Shockman & Byrne (2021): S ranging from 19.6 to 41 and T ranging between 15 to 35oC.
}

The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length. 

The pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).}


\value{
  \item{K2}{Second dissociation constant of carbonic acid (mol/kg)}

}
\references{The Lueker et al. (2000) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). The Roy et al. (1993) constants is recommended by DOE (1994).

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006 Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Millero F. J., 2010 Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Mojica Prieto F. J. and Millero F. J., 2002. The values of pK1 + pK2 for the dissociation of carbonic acid in seawater. \emph{Geochimica et Cosmochimica Acta} \bold{66}, 2529-2540.

Papadimitriou S., Loucaides S., Rérolle V. M. C., Kennedy P., Achterberg E. P., Dickson A. G., Mowlem M. and Kennedy H., 2018. The stoichiometric dissociation constants of carbonic acid in seawater brines from 298 to 267 K. \emph{Geochimica et Cosmochimica Acta} \bold{220}, 55-70.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Schockman, K.M., Byrne, R.H., 2021. Spectrophotometric Determination of the Bicarbonate Dissociation Constant in Seawater, \emph{Geochimica et Cosmochimica Acta}

Sulpis O., Lauvset S. K. and Hagens M., 2020. Current estimates of K1* and K2* appear inconsistent with measured CO2 system parameters in cold oceanic regions. \emph{Ocean Science} \bold{16}, 847-862.

Waters, J., Millero, F. J., and Woosley, R. J., 2014. Corrigendum to ``The free proton concentration scale for seawater pH'', [MARCHE: 149 (2013) 8-22], \emph{Marine Chemistry} \bold{165}, 66-67.
}

\author{
Jean-Marie Epitalon, Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K1}}.
}

\examples{
  K2(35,25,0)
}
\keyword{utilities}
