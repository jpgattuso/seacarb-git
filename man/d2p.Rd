\encoding{latin1}
\name{d2p}
\alias{d2p}
\title{Converts depth in meters to pressure in dbar}
\description{Converts depth in meters to pressure in dbar}
\usage{
d2p(depth, lat=40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depth}{Depth in meters}
  \item{lat}{Latitude in degrees, N and S is irrelevant, default is 40o}
  }

\value{
  \item{pressure}{Pressure corresponding to the depth given, in dbar}
}

\references{
Saunders P. M., 1981. Practical conversion of pressure to depth. \emph{J. Phys. Oceanogr.} \bold{11}: 573-574.
}
             
\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{p2d}}
}

\examples{
d2p(depth=7500, lat=30)  
}

\keyword{utilities}
