\encoding{latin1}
\name{teos2eos_geo}
\alias{teos2eos_geo}

\title{Convert temperature and salinity from TEOS-10 to EOS-80}
\description{Convert conservative temperature to in-situ temperature and Absolute Salinity (SA) to Practical Salinity (SP).
Salinity conversion depends on depth and geographic location. 
}
\usage{
teos2eos_geo(SA, CT, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute Salinity in g/kg}
  \item{CT}{Conservative Temperature in degrees C}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{latitude in decimal degrees [-90 ... 90]}
  }

\value{The function returns a data frame containing the following columns:
  \item{T}{in situ Temperature (deg C)}
  \item{SP}{Practical Salinity (psu)}
}

\details{
Conversion from Absolute (SA) to Practical Salinity (SP) depends on water density anomaly which is correlated with Silicate concentration.
This function relies on silicate concentration taken from WOA (World Ocean Atlas) to evaluate density anomaly.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

T. J. McDougall 1 , D. R. Jackett 2 , F. J. Millero 3 , R. Pawlowicz 4 , and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}

\seealso{
eos2teos_geo does the reverse, teos2eos_chem, sa2sp_geo

package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate insitu Temperature and practical salinity of a sample with 
   # Absolute Salinity of 35 g/kg, Conservative Temperature of 18 deg C,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   f <- teos2eos_geo(SA=35, CT=18, P=10, long=188, lat=4)
   T <- f$T     # insitu temperature
   SP <- f$SP     # Practical Salinity
}

\keyword{utilities}
