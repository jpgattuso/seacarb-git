\encoding{latin1}
\name{sp2sa_geo}
\alias{sp2sa_geo}

\title{From practical to absolute salinity}
\description{Converts from practical (SP) to absolute salinity (SA) based on depth and geographic location.}

\usage{
sp2sa_geo(SP, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{SP}{Practical salinity on the practical salinity scale}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{Latitude in decimal degrees [-90 ... 90]}
  }

\value{
  \item{SA}{Absolute salinity (g/kg)}
}

\details{
This function is almost an alias of function gsw_SA_from_SP of the gsw package on which it relies. The only difference is in that depth and location are optional.  If location is not given, or incomplete (either longitude or latitude missing), an arbitrary location is chosen: the mid equatorial atlantic ocean. Note that this implies an error on computed SA ranging from 0 up to 0.02 g/kg.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
             
\seealso{
sa2sp_geo does the reverse, sp2sa_chem
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate the absolute salinity of a sample whose practical Salinity is 35,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   SA <- sp2sa_geo(35, 10, 188, 4)     # 34.711778344814114
}

\keyword{utilities}
