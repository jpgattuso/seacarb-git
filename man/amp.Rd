\encoding{latin1}
\name{amp}
\alias{amp}

\title{pH value of the AMP buffer}
\description{pH value of the AMP buffer (on the total scale in mol/kg)}
\usage{
amp(S=35,T=25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
  \item{AMP}{pH value of the AMP buffer (on the total scale in mol/kg)}

}
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}

\author{
Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{tris}}, \code{\link{pHslope}}, \code{\link{pH}}.
}

\examples{
	##Example from Dickson et al. (2007)
	amp(S=35,T=25)
}
\keyword{utilities}
