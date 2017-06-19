\encoding{latin1}
\name{buffergen}
\alias{buffergen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Buffer factors of the seawater carbonate system as defined by Hagens and 
  Middelburg (2016)}
\description{Returns the suite of buffer factors presented in Table 3 of Hagens and Middelburg 
  (2016), as well as the proton concentration buffer factor (beta.H of Hofmann et al, 2010) and 
  the classic Revelle factor. For practical purposes, this function excludes the nitrate and 
  nitrite acid-base systems presented in this paper, as well as the fully protonoted form of 
  sulfate (H2SO4) and fully deprotonated form of sulfide (S2-), as their contributions to total 
  alkalinity under natural seawater conditions are negligible. Its input arguments are 
  indentical to those in the “carbfull” function of seacarb.}
\usage{buffergen(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, k1k2="x", kf="x",
                 ks="d", pHscale="T", b="u74", gas="potential", NH4t=0, HSt=0)}
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
	\item{var1}{Value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{Value of the second  variable in mol/kg, except for pH}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
  \item{P}{Hydrostatic pressure in bar (surface = 0)}
  \item{Pt}{Concentration of total phosphate in mol/kg; set to 0 if NA}
  \item{Sit}{Concentration of total silicate in mol/kg; set to 0 if NA}
  \item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010), "w14" from Waters et al. (2014), and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m10".}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{ks}{"d" for using Ks from Dickson (1990) and "k" for using Ks from Khoo et al. (1977), default is "d"} 
	\item{pHscale}{"T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
	\item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74" }
  \item{gas}{used to indicate the convention for INPUT pCO2, i.e., when it is an input variable (flags 21 to 25): "insitu" indicates it is referenced to in situ pressure and in situ temperature; "potential" indicates it is referenced to 1 atm pressure and potential temperature; and "standard" indicates it is referenced to 1 atm pressure and in situ temperature. All three options should give identical results at surface pressure. This option is not used when pCO2 is not an input variable (flags 1 to 15). The default is "potential".}
  \item{NH4t}{Concentration of total ammonium in mol/kg; set to 0 if NA}
  \item{HSt}{Concentration of total hydrogen sulfide in mol/kg; set to 0 if NA}
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

\item For K1, K2, pK1, pK2, pK3, Kw, Kb, Khs, Ksi and K2si, the pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

\item For Kf, the pressure correction was applied on the free scale. The formulation of Dickson and Riley (1979 in Dickson and Goyet, 1994) provides Kf on the free scale but that of Perez and Fraga (1987) provides it on the total scale. Hence, in that case, Kf was first transformed from the total scale to the free scale. With both formulations, the pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Ks, the pressure correction was applied on the free scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Kn, the pressure correction was applied on the seawater scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).
}
}

\value{The function returns a list containing the following matrices:
  \item{Carbfull}{Output of the 'carbfull' function that is used within 'buffergen'}
  \item{dALK.dH}{Sensitivity of ALK to a change in proton concentration (dimensionless). Species-specific.}
  \item{dtotX.dH}{Sensitivity of an acid-base species to a change in proton concentration (dimensionless). Species-specific.}
  \item{dALK.dX}{Sensitivity of ALK to a change in an acid-base species (dimensionless). Species-specific.}
  \item{dtotX.dX}{Sensitivity of an acid-base species to a change in its total concentration (dimensionless). Species-specific.}
  \item{dALK.dpH}{Sensitivity of ALK to a change in pH (mol/kg-soln). Species-specific.}
  \item{dtotX.dpH}{Sensitivity of an acid-species to a change in pH (mol/kg-soln). Species-specific.}
  \item{dH.dALK}{Sensitivity of proton concentration to a change in ALK (dimensionless). Values are the same for all species and all acid-base systems, except for the fluoride and sulfate acid-base systems, which slightly deviate due to pH scale conversion effects.}
  \item{dH.dtotX}{Sensitivity of an acid-species to a change in its total concentration (dimensionless). Values are the same for all species of a specific acid-base system.}
  \item{dX.dALK}{Sensitivity of an acid-species to a change in its total concentration (dimensionless). Species-specific.}
  \item{dX.dtotX}{Sensitivity of an acid-species to a change in its total concentration (dimensionless). Species-specific.}
  \item{dpH.dALK}{Sensitivity of pH due to a change in ALK ((mol/kg-soln)-1). Values are the same for all species and all acid-base systems, except for the fluoride and sulfate acid-base systems, which slightly deviate due to pH scale conversion effects.}
  \item{dpH.dtotX}{Sensitivity of pH due to a change in the total concentration of an acid-base system ((mol/kg-soln)-1). Values are the same for all species of a specific acid-base system.}
  \item{beta.H}{proton concentration buffer factor (Eq.4 of Hagens and Middelburg (2016), dimensionless)}   
  \item{RF}{Revelle factor (dimensionless)}
  If the total concentration of an acid-base system is 0, the values of the buffer factors corresponding to all species of that acid-base system return NA.
}

\references{
Hagens M. and Middelburg J. J., 2016 Generalised expressions for the response of pH to changes in ocean chemistry. \emph{Geochimica et Cosmochimica Acta}.
}

\author{
Mathilde Hagens \email{m.hagens@uu.nl}
}

\examples{

## With a couple of variables
buffergen(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="potential", NH4t=0, HSt=0)

## With a couple of variables and non-zero nutrient concentrations
buffergen(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Patm=1.0, Pt=5e-6, Sit=2e-6,
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
buffergen(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, 
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
buffergen(flag=flag, var1=var1, var2=var2)
}
\keyword{utilities}
