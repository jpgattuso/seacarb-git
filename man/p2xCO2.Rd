\encoding{latin1}
\name{p2xCO2}
\alias{p2xCO2}
\title{Converts partial pressure of CO2 to mole fraction of CO2}
\description{Converts pCO2 (partial pressure of CO2) into xCO2 (mole fraction of CO2)}
\usage{
p2xCO2(S=35, T=25, Patm=1, pCO2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{Patm}{Atmospheric pressure in atmospheres, default is 1.0 (this is not the hydrostatic pressure)}
  \item{pCO2}{Partial pressure of CO2 in \eqn{\mu}atm}
  }

\value{
  \item{xCO2}{Mole fraction of CO2 in ppm.}
}

\details{
The xCO2 (ppm) is computed from pCO2 (\eqn{\mu}atm) using the
following equation: \eqn{ xCO2 = pCO2 / (Patm - pH2O)}, where pH20 is the vapor pressure of seawater computed following best practices (Dickson et al., 2007). That computed pH20 is identical, when rounded to the 4th decimal place, with that computed by the equation from Weiss and Price (1980).}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2014)}

\references{
Dickson A. G., Sabine C. L. and Christian J. R. (2007) Guide to best
practices for ocean CO2 measurements. \emph{PICES Special Publication}
\bold{3}, 1-191.

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
	\code{\link{x2pCO2}} and \code{\link{vapress}}
}

\examples{
## Convert atmospheric pressure from mbar to atm
   Patm_mbar = 1052                 # in millibar
   Patm      = Patm_mbar / 1013.25  # in atm
## Compute xCO2 from pCO2
   pCO2 = 380
   xCO2 = p2xCO2(T=25, S=35, Patm=Patm, pCO2=pCO2)
   print(xCO2)
## The result is 377.1546 ppm 
}

\keyword{utilities}
