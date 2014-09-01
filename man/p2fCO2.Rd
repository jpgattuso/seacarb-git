\encoding{latin1}
\name{p2fCO2}
\alias{p2fCO2}
\title{Converts pCO2 (partial pressure in CO2) into fCO2 (fugacity of CO2)}
\description{Converts pCO2 (partial pressure in CO2) into fCO2 (fugacity of CO2)}
\usage{
p2fCO2(T = 25, pCO2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{pCO2}{Partial pressure in CO2, the unit given here define the unit of fugacity in output}
  }

\value{
  \item{fCO2}{Fugacity of CO2, given in the same unit than pCO2, entered in argument.}
}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}
             
\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{f2pCO2}}.
}

\examples{
p2fCO2(T=25, pCO2=380)  
}

\keyword{utilities}
