\encoding{latin1}
\name{f2pCO2}
\alias{f2pCO2}
\title{Converts the CO2 fugacity to CO2 partial pressure}
\description{Converts fCO2 (fugacity of CO2) into pCO2 (partial pressure in CO2)}
\usage{
f2pCO2(T = 25, fCO2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{fCO2}{Fugacity of CO2 in \eqn{\mu}atm}
  }

\value{
  \item{pCO2}{Partial pressure of CO2 in \eqn{\mu}atm.}
}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}
             
\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{p2fCO2}}.
}

\examples{
f2pCO2(T=25, fCO2=380)  
}

\keyword{utilities}
