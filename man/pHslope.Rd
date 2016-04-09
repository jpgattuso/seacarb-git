\encoding{latin1}
\name{pHslope}
\alias{pHslope}

\title{Slope of the calibration curve of a pH electrode}
\description{Slope of the calibration curve of a pH electrode (percent of theoretical slope)}
\usage{
pHslope(Etris=-72.4,Eamp=4.9,S=35,T=25)
%- maybe also 'usage' for other objects documented here.
}

\arguments{
  \item{Etris}{e.m.f. of the TRIS buffer in mV, default is -72.4}
  \item{Eamp}{e.m.f. of the AMP buffer in mV, default is 4.9}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
  \item{pHslope}{Slope of the calibration curve (in percent of theoretical slope)}

}
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}

\author{
Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{tris}}, \code{\link{amp}}, \code{\link{pH}}.
}

\examples{
	##Example from Dickson et al. (2007)
	pHslope(Etris=-72.4,Eamp=4.9,S=35,T=25)
}
\keyword{utilities}
