\encoding{latin1}
\name{eos2teos_geo}
\alias{eos2teos_geo}

\title{Convert temperature and salinity from EOS-80 to TEOS-10}
\description{Convert in-situ to Conservative temperature and Practical (SP) to Absolute Salinity (SA)
Salinity conversion depends on depth and geographic location. 
}
\usage{
eos2teos_geo(SP, T, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical Salinity on the practical salinity scale}
  \item{T}{in-situ temperature in deg. C}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{latitude in decimal degrees [-90 ... 90]}
  }

\value{The function returns a data frame containing the following columns:
  \item{CT}{Conservative Temperature (deg C)}
  \item{SA}{Absolute Salinity (g/kg)}
}

\details{
Conversion from Practical (SP) to Absolute Salinity (SA) depends on water density anomaly which is correlated with Silicate concentration.
This function relies on silicate concentration taken from WOA (World Ocean Atlas) to evaluate density anomaly.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

T. J. McDougall 1 , D. R. Jackett 2 , F. J. Millero 3 , R. Pawlowicz 4 , and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
teos2eos_geo does the reverse, eos2teos_chem, sp2sa_geo

package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate Conservative Temperature and Absolute Salinity of a sample with 
   # Practical Salinity of 35 psu, in-situ Temperature of 18 deg C,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   f <- eos2teos_geo(SP=35, T=18, P=10, long=188, lat=4)
   CT <- f$CT     # Conservative Temperature
   SA <- f$SA     # Absolute Salinity
}

\keyword{utilities}
