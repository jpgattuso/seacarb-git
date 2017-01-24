\encoding{latin1}
\name{sp2sa_geo}
\alias{sp2sa_geo}

\title{From Practical to Absolute Salinity}
\description{Convert from Practical (SP) to Absolute Salinity (SA)
based on depth and geographic location.}
\usage{
sp2sa_geo(SP, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical Salinity on the practical salinity scale}
  \item{P}{Sea water pressure in dbar}
  \item{lon}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{latitude in decimal degrees [-90 ... 90]}
  }

\value{
  \item{SA}{Absolute Salinity (g/kg)}
}

\details{
This subroutine is almost an alias of subroutine gsw_SA_from_SP of gsw package on which it relies.
The only difference is in that depth and location are optional. 
If location is not given, or incomplete (either longitude or latitude missing),
an arbitrary location is chosen: the mid equatorial atlantic ocean.
Note that this implies an error on computed SA ranging from 0 up to 0.02 g/kg
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

T. J. McDougall 1 , D. R. Jackett 2 , F. J. Millero 3 , R. Pawlowicz 4 , and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
sa2sp_geo does the reverse

sp2sa_chem
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
