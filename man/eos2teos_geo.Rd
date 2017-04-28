\encoding{latin1}
\name{eos2teos_geo}
\alias{eos2teos_geo}

\title{Convert temperature and salinity from EOS-80 to TEOS-10}
\description{Converts in situ to conservative temperature and practical (SP) to absolute salinity (SA). Salinity conversion depends on depth and geographic location. 
}
\usage{
eos2teos_geo(SP, T, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical salinity on the practical salinity scale}
  \item{T}{In situ temperature in deg. C}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{Latitude in decimal degrees [-90 ... 90]}
  }

\value{The function returns a data frame containing the following columns:
  \item{CT}{Conservative temperature (deg C)}
  \item{SA}{Absolute salinity (g/kg)}
}

\details{
Conversion from Practical (SP) to Absolute Salinity (SA) depends on water density anomaly which is correlated with silicate concentration. This function relies on silicate concentration taken from WOA (World Ocean Atlas) to evaluate density anomaly.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
teos2eos_geo does the reverse, eos2teos_chem, sp2sa_geo, package gsw
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
