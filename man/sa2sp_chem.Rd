\encoding{latin1}
\name{sa2sp_chem}
\alias{sa2sp_chem}

\title{From Absolute to Practical Salinity}
\description{Converts from absolute (SA) to practical salinity (SP). Salinity conversion depends on total alkalinity as well as the concentrations of dissolved inorganic carbon, nitrate and silicate.}
\usage{
sa2sp_chem(SA, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{Absolute salinity in g/kg}
  \item{TA}{Total alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved inorganic carbon concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total silicate concentration in mol/kg, default is 0}
  }

\value{
  \item{SP}{Practical salinity (in psu)}
}

\details{
Convert from Absolute (SA) to Practical Salinity (SP)
from carbon system parameters and ion concentration which most affect water density anomalies.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., Wright D. G. and Millero F. J., 2011. The effects of biogeochemical processes on oceanic conductivity/salinity/density relationships and the characterization of real seawater. \emph{Ocean Science} \bold{7}, 363-387.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
sp2sa_chem does the reverse, sa2sp_geo
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
