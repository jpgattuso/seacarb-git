\encoding{latin1}
\name{p2d}
\alias{p2d}
\title{Converts pressure in dbar to depth in meters}
\description{Converts pressure in dbar to depth in meters}
\usage{
p2d(pressure, lat=40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pressure}{Pressure in dbar}
  \item{lat}{Latitude in degrees, N and S is irrelevant, default is 40o}
  }

\value{
  \item{depth}{Depth corresponding to the pressure given, in meters}
}

\references{
Saunders P. M., 1981. Practical conversion of pressure to depth. \emph{J. Phys. Oceanogr.} \bold{11}: 573-574.
}
             
\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{d2p}}
}

\examples{
p2d(pressure=7686, lat=30)  
}

\keyword{utilities}

