\name{sir_full}
\alias{sir_full}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameters of the seawater carbonate system including the substrate-inhibitor-ratio (SIR) - extension of carb
}
\description{
This function is a modification of carb. It returns parameters of the seawater carbonate system, including the ammonium and sulfide acid-base systems, as well as full acid-base speciation and additional parameters relevant to the SIR.
}
\usage{
sir_full(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0,
        k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
        NH4t=0, HSt=0)
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
	\item{var1}{Value of the first variable in mol/kg-soln, except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{Value of the second variable in mol/kg-soln, except for pH}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
  \item{P}{Hydrostatic pressure in bar (surface = 0)}
  \item{Pt}{Concentration of total phosphate in mol/kg-soln; set to 0 if NA}
  \item{Sit}{Concentration of total silicate in mol/kg-soln; set to 0 if NA}
  \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "mp2" from Mojica Prieto et al. (2002), "m02" from Millero et al. (2002), "m06" from Millero et al. (2006), "m10" from Millero (2010), "w14" from Waters et al. (2014), "p18" from Papadimitriou et al. (2018), "s20" from Sulpis et al. (2020), and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m10".}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickson (1990) and "k" for using Ks from Khoo et al. (1977), default is "d"} 
	\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
	\item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74" }
  \item{gas}{used to indicate the convention for INPUT pCO2, i.e., when it is an input variable (flags 21 to 25): "insitu" indicates it is referenced to in situ pressure and in situ temperature; "potential" indicates it is referenced to 1 atm pressure and potential temperature; and "standard" indicates it is referenced to 1 atm pressure and in situ temperature. All three options should give identical results at surface pressure. This option is not used when pCO2 is not an input variable (flags 1 to 15). The default is "potential".}
  \item{NH4t}{Concentration of total ammonium in mol/kg-soln; set to 0 if NA}
  \item{HSt}{Concentration of total hydrogen sulfide in mol/kg-soln; set to 0 if NA}

}
\details{Calculates the "substrate-inhibitor ratio (SIR)" (i.e. [HCO3-]/[H+]) which is used as a metric to predict carbonate chemistry dependency of biotic calcium carbonate formation (Bach, 2015).  Please note that hydrogen ion concentrations in the output [HCO3-]/[H+] are on the \emph{free scale}, regardless of the input pH scale. We have included pH and [H+] outputs on multiple scales as a teaching exercise to show how calculating the SIR on different scales changes its meaning.

The function is based on \emph{carbfull} and therefore also returns parameters of the seawater carbonate system, including the ammonium and sulfide acid- base systems, as well as full acid-base speciation.

\bold{\emph{carbfull} details}:

The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\emph{For K1 and K2:}
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

\item For K1, K2, pK1, pK2, pK3, Kw, Kb, Khs, Ksi and K2si, the pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

\item For Kf, the pressure correction was applied on the free scale. The formulation of Dickson and Riley (1979 in Dickson and Goyet, 1994) provides Kf on the free scale but that of Perez and Fraga (1987) provides it on the total scale. Hence, in that case, Kf was first transformed from the total scale to the free scale. With both formulations, the pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Ks, the pressure correction was applied on the free scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Kn, the pressure correction was applied on the seawater scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).
}
}


\value{
The function returns a data frame containing the following columns:
  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm}
  \item{P}{Hydrostatic pressure in bar}
  \item{pH}{pH}
  \item{CO2}{CO2 concentration (mol/kg-soln)}
  \item{pCO2}{"standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{fCO2}{"standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{pCO2pot}{"potential" pCO2, CO2 partial pressure computed at potential temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{fCO2pot}{"potential" fCO2, CO2 fugacity computed at potential temperature and atmospheric pressure (\eqn{\mu}atm)}
  \item{pCO2insitu}{"in situ" pCO2, CO2 partial pressure computed at in situ temperature and total pressure (atm + hydrostatic) (\eqn{\mu}atm)}
  \item{fCO2insitu}{"in situ" fCO2, CO2 fugacity computed at in situ temperature and total pressure (atm + hydrostatic) (\eqn{\mu}atm)}
  \item{HCO3}{HCO3 concentration (mol/kg-soln)}
  \item{CO3}{CO3 concentration (mol/kg-soln)}
  \item{DIC}{DIC concentration (mol/kg-soln)}
  \item{ALK}{ALK, total alkalinity (mol/kg-soln)}
  \item{OmegaAragonite}{Omega aragonite, aragonite saturation state}
  \item{OmegaCalcite}{Omega calcite, calcite saturation state}
  \item{NH4}{NH4 concentration (mol/kg-soln)} 
  \item{NH3}{NH3 concentration (mol/kg-soln)}
  \item{BOH3}{B(OH)3 concentration (mol/kg-soln)} 
  \item{BOH4}{B(OH)4 concentration (mol/kg-soln)} 
  \item{H3PO4}{H3PO4 concentration (mol/kg-soln)} 
  \item{H2PO4}{H2PO4 concentration (mol/kg-soln)} 
  \item{HPO4}{HPO4 concentration (mol/kg-soln)} 
  \item{PO4}{PO4 concentration (mol/kg-soln)} 
  \item{H2S}{H2S concentration (mol/kg-soln)} 
  \item{HS}{HS concentration (mol/kg-soln)} 
  \item{SiOH4}{SiOH4 concentration (mol/kg-soln)} 
  \item{SiOOH3}{SiOOH3 concentration (mol/kg-soln)} 
  \item{SiO2OH2}{SiO2OH2 concentration (mol/kg-soln)} 
  \item{HF}{HF concentration (mol/kg-soln)}
  \item{F}{F concentration (mol/kg-soln)} 
  \item{HSO4}{HSO4 concentration (mol/kg-soln)} 
  \item{SO4}{SO4 concentration (mol/kg-soln)}
  \item{H}{H concentration at specified pH scale (mol/kg-soln)}
  \item{OH}{OH concentration (mol/kg-soln)}
  \item{NH4t}{Supplied NHt concentration (mol/kg-soln); values are recycled if necessary} 
  \item{BOR}{Calculated total borate concentration (mol/kg-soln)} 
  \item{Pt}{Supplied Pt concentration (mol/kg-soln); values are recycled if necessary} 
  \item{HSt}{Supplied HSt concentration (mol/kg-soln); values are recycled if necessary}  
  \item{Sit}{Supplied Sit concentration (mol/kg-soln); values are recycled if necessary} 
  \item{FLUO}{Calculated total fluoride concentration (mol/kg-soln)} 
  \item{ST}{Calculated total sulfate concentration (mol/kg-soln)} 
  \item{SIR}{substrate inhibitor ratio: [HCO3]/[H_free]  (mol:\eqn{\mu}mol)}
  \item{H_free}{H ion concentration on free scale (\eqn{\mu}mol/kg)}
  \item{H_sws}{H ion concentration on seawater scale (\eqn{\mu}mol/kg)}
  \item{H_t}{H ion concentration on total scale (\eqn{\mu}mol/kg)}
  \item{pH_free}{pH calculated on the free scale}
  \item{pH_sws}{pH calculated on the seawater scale}
  \item{pH_t}{pH calculated on the total scale}

}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2014)}

\references{

Bach L. T., 2015 Reconsidering the role of carbonate ion concentration in calcification by marine organisms. \emph{Biogeosciences} \bold{12}, 4939-4951.

Dickson A. G. and Riley J. P., 1979 The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water. \emph{Marine Chemistry} \bold{7}, 89-99.

Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen Ion Concentration in Seawater from 5 to 40oC: Standard Potentials at Salinities from 20 to 45. \emph{Analytical Chemistry} \bold{22}, vol49 29-34. 

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.

Millero F. J., 2010. Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006. Dissociation constants of carbonic acid in seawater as a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2014. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences Discussions} \bold{11}, 5327-5397.

Perez F. F. and Fraga F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.

Waters, J., Millero, F. J., and Woosley, R. J., 2014. Corrigendum to ``The free proton concentration scale for seawater pH'', [MARCHE: 149 (2013) 8-22], Mar. Chem. 165, 66-67.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of a non-ideal gas, \emph{Mar.  Chem.}, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and seawater, \emph{Mar. Chem.}, \bold{8}, 347-359.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.

}
\author{
 Kimberlee Baldry, Lennart Bach and Jean-Pierre Gattuso, 
\email{gattuso@obs-vlfr.fr}
}

\examples{
## With a couple of variables
sir_full(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=0, HSt=0)

## With a couple of variables and non-zero nutrient concentrations
sir_full(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=5e-6, Sit=2e-6,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=10e-6, HSt=0.1e-6)

## Using vectors as arguments
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
b <- c("l10", "l10", "l10")
gas <- c("potential", "potential", "potential")
NH4t <- c(0, 0, 0)
HSt <- c(0, 0, 0)
sir_full(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
  kf=kf, k1k2=k1k2, pHscale=pHscale, b=b, gas=gas, NH4t=NH4t, HSt=HSt)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
sir_full(flag=flag, var1=var1, var2=var2)
}
\keyword{utilities}

