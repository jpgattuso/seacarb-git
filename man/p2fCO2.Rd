\encoding{latin1}
\name{p2fCO2}
\alias{p2fCO2}
\title{Converts pCO2 (partial pressure in CO2) into fCO2 (fugacity of CO2)}
\description{Converts pCO2 (partial pressure in CO2) into fCO2 (fugacity of CO2)}
\usage{
p2fCO2(T = 25, Patm=1, P=0, pCO2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
  \item{P}{Hydrostatic pressure in bar, default is 0 bar (surface)}
  \item{pCO2}{Partial pressure in CO2 in \eqn{\mu}atm, the same units as that for the fugacity output}
  }

\value{
  \item{fCO2}{Fugacity of CO2 in \eqn{\mu}atm, the same units as that for the pCO2 input.}
}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2015)}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2015. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences} \bold{12}, 1483-1510.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of a non-ideal gas, \emph{Mar.  Chem.}, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and seawater, \emph{Marine Chemistry}, \bold{8}, 347-359.
}
             
\author{
Heloise Lavigne, Jean-Pierre Gattuso, and James Orr \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{f2pCO2}}.
}

\examples{
p2fCO2(T=25, Patm=0.97, P=0, pCO2=380)  
}

\keyword{utilities}
