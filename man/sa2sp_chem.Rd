\encoding{latin1}
\name{sa2sp_chem}
\alias{sa2sp_chem}

\title{From Absolute to Practical Salinity}
\description{Convert from Absolute (SA) to Practical Salinity (SP)
based on Total Alkalinity, Dissolved Inorganic Carbon 
and Nitrate and Silicate concentrations.}
\usage{
sa2sp_chem(SA, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute Salinity in g/kg}
  \item{TA}{Total Alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved Inorganic Carbone concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total Nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total Silicate concentration in mol/kg, default is 0}
  }

\value{
  \item{SP}{Practical Salinity (in psu)}
}

\details{
Convert from Absolute (SA) to Practical Salinity (SP)
from carbon system parameters and ion concentration which most affect water density anomalies.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The effects of biogeochemical processes on oceanic conductivity/
salinity/density relationships and the characterization of real seawater

T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
sp2sa_chem does the reverse

sa2sp_geo
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate the practical salinity of a sample with Absolute Salinity of 35 g/kg,
   # Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   SP <- sa2sp_chem(SA=35, TA=0.00234, DIC=0.00202)
}

\keyword{utilities}
