\encoding{latin1}
\name{eos2teos_chem}
\alias{eos2teos_chem}

\title{Convert temperature and salinity from EOS-80 to TEOS-10}
\description{Convert in-situ to Conservative temperature and Practical (SP) to Absolute Salinity (SA)
Salinity conversion depends on Total Alkalinity, Dissolved Inorganic Carbon 
and Nitrate and Silicate concentrations.}
\usage{
eos2teos_chem(SP, T, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical Salinity on the practical salinity scale}
  \item{T}{in-situ temperature in deg. C}
  \item{P}{Sea water pressure in dbar}
  \item{TA}{Total Alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved Inorganic Carbone concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total Nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total Silicate concentration in mol/kg, default is 0}
  }

\value{The function returns a data frame containing the following columns:
  \item{CT}{Conservative Temperature (deg C)}
  \item{SA}{Absolute Salinity (g/kg)}
}

\details{
Conversion from Practical (SP) to Absolute Salinity (SA)
depends on carbon system parameters and ion concentration which most affect water density anomalies.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

R. Pawlowicz 1 , D. G. Wright 2,† , and F. J. Millero, 2011: The effects of biogeochemical processes on oceanic conductivity/
salinity/density relationships and the characterization of real seawater

T. J. McDougall 1 , D. R. Jackett 2 , F. J. Millero 3 , R. Pawlowicz 4 , and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
teos2eos_chem does the reverse, eos2teos_geo, sp2sa_chem

package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate Conservative Temperature and Absolute Salinity of a sample with 
   # Practical Salinity of 35 psu, in-situ Temperature of 18 deg C,
   # at 0 dbar and Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- eos2teos_chem(SP=35, T=18, P=0, TA=0.00234, DIC=0.00202)
   CT <- f$CT     # Conservative Temperature
   SA <- f$SA     # Absolute Salinity
}

\keyword{utilities}
