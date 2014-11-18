\encoding{latin1}
\name{f2pCO2}
\alias{f2pCO2}
\title{Converts the CO2 fugacity to CO2 partial pressure}
\description{Converts fCO2 (fugacity of CO2) into pCO2 (partial pressure in CO2)}
\usage{
f2pCO2(T = 25, Patm=1, P=0, fCO2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{Patm}{Surface atmospheric pressure in atm, default 1 atm}
  \item{P}{Hydrostatic pressure in bar, default is 0 bar (surface)}
  \item{fCO2}{Fugacity of CO2 in \eqn{\mu}atm, the same units as that for the pCO2 output}
  }

\value{
  \item{pCO2}{Partial pressure of CO2 in \eqn{\mu}atm, the same units as that for the fCO2 input}
}

\note{\bold{Warning:} pCO2 estimates below 100 m are subject to considerable uncertainty. See Weiss (1974) and Orr et al. (2014)}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007. Guide to best
practices for ocean CO2 measurements. \emph{PICES Special Publication}
\bold{3}, 1-191.

Orr J. C., Epitalon J.-M. and Gattuso J.-P., 2014. Comparison of seven packages that compute ocean carbonate chemistry. \emph{Biogeosciences Discussions} \bold{11}, 5327-5397.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and
seawater, \emph{Mar. Chem.}, \bold{8}, 347-359, 1980.

Weiss, R. F., 1974. Carbon dioxide in water and seawater: the solubility of
a non-ideal gas, \emph{Mar.  Chem.}, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A., 1980. Nitrous oxide solubility in water and
seawater, \emph{Mar. Chem.}, \bold{8}, 347-359.
}
             
\author{
Heloise Lavigne, Jean-Pierre Gattuso, and James Orr \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{p2fCO2}}.
}

\examples{
f2pCO2(T=25, Patm=1.0, P=0, fCO2=380)  
}

\keyword{utilities}
