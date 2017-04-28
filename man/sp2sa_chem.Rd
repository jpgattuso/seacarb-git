\encoding{latin1}
\name{sp2sa_chem}
\alias{sp2sa_chem}

\title{From Practical to Absolute Salinity}
\description{Convert from Practical (SP) to Absolute Salinity (SA)
based on Total Alkalinity, Dissolved Inorganic Carbon 
and Nitrate and Silicate concentrations.}
\usage{
sp2sa_chem(SP, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical Salinity on the practical salinity scale}
  \item{TA}{Total Alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved Inorganic Carbone concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total Nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total Silicate concentration in mol/kg, default is 0}
  }

\value{
  \item{SA}{Absolute Salinity (g/kg)}
}

\details{
Convert from Practical (SP) to Absolute Salinity (SA)
from carbon system parameters and ion concentration which most affect water density anomalies.
}

\references{
TEOS-10 web site: http://www.teos-10.org/

What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer) by Rich Pawlowicz (on TEOS-10 web site)

R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The effects of biogeochemical processes on oceanic conductivity/
salinity/density relationships and the characterization of real seawater

T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, and P. M. Barker, 2012: Algorithm for estimating Absolute Salinity}
             
\seealso{
sa2sp_chem does the reverse

sp2sa_geo
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate the absolute salinity of a sample with practical Salinity of 35,
   # Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   SA <- sp2sa_chem(SP=35, TA=0.00234, DIC=0.00202)
}

\keyword{utilities}
