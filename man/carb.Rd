\encoding{latin1}
\name{carb}
\alias{carb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameters of the seawater carbonate system}
\description{Returns parameters of the seawater carbonate system.}
\usage{carb(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0,
        k1k2="x", kf="x", ks="d", pHscale="T", b="l10")}
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
        \item{Patm}{Surface atmospheric pressure in atm, default 1 atm}
  	\item{P}{Hydrostatic pressure in bar (surface = 0)}
  	\item{Pt}{Concentration of total phosphate in mol/kg; set to 0 if NA}
  	\item{Sit}{Concentration of total silicate in mol/kg; set to 0 if NA}
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m06" from Millero et al. (2006), "m10" from Millero (2010) and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m10".}
 F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I}, \bold{21} 161-162.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of
a non-ideal gas, \emph{Mar.  Chem.}, \bold{2}, 203–215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and
seawater, \emph{Mar. Chem.}, \bold{8}, 347–359.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}

\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{

## With a couple of variables
carb(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Pt=0, Sit=0,
	pHscale="T", kf="pf", k1k2="l", ks="d", b="l10")

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
carb(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P,
	Pt=Pt, Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

## Test with all flags 
flag <- c((1:15), (21:25))
var1 <- c(8.200000, 7.308171e-06, 7.308171e-06, 7.308171e-06, 7.308171e-06, 
	8.2, 8.2, 8.2, 8.2, 0.001646857, 0.001646857, 0.001646857, 0.0002822957, 
	0.0002822957, 0.00234, 258.2164, 258.2164, 258.2164, 258.2164, 258.2164 )
var2 <- c(7.308171e-06, 0.001646857, 0.0002822957, 0.00234, 0.001936461, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461, 0.0002822957, 
	0.00234, 0.001936461,  0.00234, 0.001936461, 0.001936461, 8.2, 
	0.001646857, 0.0002822957, 0.00234, 0.001936461)
carb(flag=flag, var1=var1, var2=var2)

## Test using a data frame 
data(seacarb_test_P0)	#test data set for P=0 (surface)
tab <- seacarb_test_P0

## method 1 using the column numbers
carb(flag=tab[[1]], var1=tab[[2]], var2=tab[[3]], S=tab[[4]], T=tab[[5]], 
	P=tab[[6]], Sit=tab[[8]], Pt=tab[[7]])

## method 2 using the column names
carb(flag=tab$flag, var1=tab$var1, var2=tab$var2, S=tab$S, T=tab$T, P=tab$P, 
	Sit=tab$Sit, Pt=tab$Pt)

}
\keyword{utilities}
