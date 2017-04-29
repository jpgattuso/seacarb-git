\encoding{latin1}
\name{sa2sp_geo}
\alias{sa2sp_geo}

\title{From Absolute to Practical Salinity}
\description{Converts from absolute (SA) to practical salinity (SP) based on depth and geographic location.}
\usage{
sa2sp_geo(SA, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute salinity in g/kg}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{Latitude in decimal degrees [-90 ... 90]}
  }

\value{
  \item{SP}{Practical salinity (psu)}
}

\details{
This function is almost an alias of subroutine gsw_SP_from_SA from gsw package on which it relies. The only difference is in that depth and location are optional. If location is not given, or incomplete (either longitude or latitude missing), an arbitrary location is chosen: the mid equatorial atlantic ocean. Note that this implies an error on computed SA ranging from 0 up to 0.02 g/kg.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
sp2sa_geo does the reverse, sa2sp_chem
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate the practical salinity of a sample whose absolute Salinity is 35,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   SP <- sa2sp_geo(35, 10, 188, 4)
}

\keyword{utilities}
