\encoding{latin1}
\name{pH}
\alias{pH}

\title{Potentiometric pH}
\description{Calculation of potentiometric pH}
\usage{
	pH(Ex=-67,Etris=-72.4,S=35,T=25)
	%- maybe also 'usage' for other objects documented here.
}

\arguments{
	\item{Ex}{e.m.f. of the seawater sample in mV, default is 67}
	\item{Etris}{e.m.f. of the TRIS buffer in mV, default is -72.4}
	\item{S}{Salinity, default is 35}
	\item{T}{Temperature in degrees Celsius, default is 25oC}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
  \item{pH}{Potentiometric pH (in mol/kg on the total scale)}

}
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}

\author{
Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{tris}}, \code{\link{amp}}, \code{\link{pHslope}}.
}

\examples{
	##Example from Dickson et al. (2007)
	pH(Ex=-67,Etris=-72.4,S=35,T=25)
}
\keyword{utilities}
