\encoding{latin1}
\name{teos2eos_chem}
\alias{teos2eos_chem}

\title{Convert temperature and salinity from TEOS-10 to EOS-80}
\description{Convert conservative temperature to in-situ temperature and Absolute Salinity (SA) to Practical Salinity (SP).
Salinity conversion depends on Total Alkalinity, Dissolved Inorganic Carbon 
and Nitrate and Silicate concentrations.}
\usage{
teos2eos_chem(SA, CT, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute Salinity in g/kg}
  \item{CT}{Conservative Temperature in degrees C}
  \item{P}{Sea water pressure in dbar}
  \item{TA}{Total Alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved Inorganic Carbone concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total Nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total Silicate concentration in mol/kg, default is 0}
  }

\value{The function returns a data frame containing the following columns:
  \item{T}{in situ Temperature (deg C)}
  \item{SP}{Practical Salinity (psu)}
}

\details{
Conversion from Absolute (SA) to Practical Salinity (SP)
depends on carbon system parameters and ion concentration which most affect water density anomalies.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

R. Pawlowicz 1 , D. G. Wright 2,† , and F. J. Millero, 2011: The effects of biogeochemical processes on oceanic conductivity/
salinity/density relationships and the characterization of real seawater

T. J. McDougall 1 , D. R. Jackett 2 , F. J. Millero 3 , R. Pawlowicz 4 , and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
eos2teos_chem does the reverse, teos2eos_geo, sa2sp_cham

package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate insitu Temperature and practical salinity of a sample with 
   # Absolute Salinity of 35 g/kg, Conservative Temperature of 18 deg C,
   # at 0 dbar and Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- teos2eos_chem(SA=35, CT=18, P=0, TA=0.00234, DIC=0.00202)
   T <- f$T     # insitu temperature
   SP <- f$SP     # Practical Salinity
}

\keyword{utilities}
