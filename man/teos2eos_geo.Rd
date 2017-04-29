\encoding{latin1}
\name{teos2eos_geo}
\alias{teos2eos_geo}

\title{Convert temperature and salinity from TEOS-10 to EOS-80}
\description{Converts conservative temperature to in situ temperature and absolute salinity (SA) to practical salinity (SP). Salinity conversion depends on depth and geographic location.
}
\usage{
teos2eos_geo(SA, CT, P=0, long=1.e20, lat=1.e20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute salinity in g/kg}
  \item{CT}{Conservative temperature in degrees C}
  \item{P}{Sea water pressure in dbar}
  \item{long}{Longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]}
  \item{lat}{Latitude in decimal degrees [-90 ... 90]}
  }

\value{The function returns a data frame containing the following columns:
  \item{T}{In situ temperature (deg C)}
  \item{SP}{Practical salinity (psu)}
}

\details{
Conversion from absolute (SA) to practical salinity (SP) depends on water density anomaly which is correlated with silicate concentration. This function relies on silicate concentration taken from WOA (World Ocean Atlas) to evaluate the density anomaly.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., Wright D. G. and Millero F. J., 2011. The effects of biogeochemical processes on oceanic conductivity/salinity/density relationships and the characterization of real seawater. \emph{Ocean Science} \bold{7}, 363-387.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}

\seealso{
eos2teos_geo does the reverse, teos2eos_chem, sa2sp_geo, package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate in situ Temperature and practical salinity of a sample with 
   # Absolute Salinity of 35 g/kg, Conservative Temperature of 18 deg C,
   # depth is 10 dbar and location is 188 degrees East and 4 degrees North.
   f <- teos2eos_geo(SA=35, CT=18, P=10, long=188, lat=4)
   T <- f$T     # insitu temperature
   SP <- f$SP     # Practical Salinity
}

\keyword{utilities}
