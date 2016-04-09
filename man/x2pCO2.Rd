\encoding{latin1}
\name{x2pCO2}
\alias{x2pCO2}
\title{Converts mole fraction to partial pressure of CO2}
\description{Converts xCO2 (mole fraction of CO2) into pCO2 (partial pressure of CO2)}
\usage{
x2pCO2(S=35, T=25, Patm=1.0, xCO2=400)
}

\arguments{
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
\item{Patm}{Atmospheric pressure in atmospheres, default is 1.0}
  \item{xCO2}{Mole fraction of CO2 in ppm, default is 400}
  }

\value{
  \item{pCO2}{Partial pressure of CO2 in \eqn{\mu}atm.}
}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2014)}

\details{
The mole fraction xCO2 (ppm) is computed from pCO2 (\eqn{\mu}atm) using the following equation: \eqn{ pCO2 = xCO2 (Patm - pH2O)}, where pH20 is the vapor pressure of seawater computed following best practices (Dickson et al., 2007). That computed pH20 is identical, when rounded to the 4th decimal place, with that computed by the equation from Weiss and Price (1980).}

\references{ 
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2014. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences Discussions} \bold{11}, 5327-5397.

Weiss, R. F. (1974) Carbon dioxide in water and seawater: the
solubility of a non-ideal gas, Mar. Chem., \bold{2}, 203-215.

Weiss, R. F. and Price, B. A. (1980) Nitrous oxide solubility in water
and seawater, Mar. Chem., \bold{8}, 347-359.
}
             
\author{
James Orr \email{james.orr@lsce.ipsl.fr}
}

\seealso{
	\code{\link{p2xCO2}} and \code{\link{vapress}}
}

\examples{
## Atmospheric pressure is rarely equal to 1 atm exactly 
## Over the Southern Ocean Patm=0.97 is more realistic
   pCO2_socn <- x2pCO2(S=35, T=0, Patm=0.97, xCO2=400.0)
   print(pCO2_socn)
## The result (385.6322 uatm) is 12 uatm less than if it was wrongly assumed that Patm=1.0

## Show effect of temperature on pCO2 computed from xCO2, and on resulting variables from "carb"
S <- 35
ALK <- 2300e-6
T <- seq(0,30,5)
xCO2 <- 400
pCO2 <- x2pCO2(S=35, T=T, Patm=1, xCO2=400)
results <- carb(flag=24, var1=pCO2, var2=ALK, S=S, T=T, P=0, Pt=0, Sit=0, 
  pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")
print(results)

}

\keyword{utilities}
