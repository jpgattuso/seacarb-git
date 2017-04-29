\encoding{latin1}
\name{sp2sa_chem}
\alias{sp2sa_chem}

\title{From Practical to absolute salinity}
\description{Converts from practical (SP) to absolute salinity (SA) based on total alkalinity as well as on the concentrations of dissolved inorganic carbon, nitrate and silicate.}
\usage{
sp2sa_chem(SP, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{SP}{Practical salinity on the practical salinity scale}
  \item{TA}{Total alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved inorganic carbon concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total silicate concentration in mol/kg, default is 0}
  }

\value{
  \item{SA}{Absolute salinity (g/kg)}
}

\details{
Converts from practical (SP) to absolute salinity (SA) from carbonate system parameters and ion concentration which mostly affect water density anomalies.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., Wright D. G. and Millero F. J., 2011. The effects of biogeochemical processes on oceanic conductivity/salinity/density relationships and the characterization of real seawater. \emph{Ocean Science} \bold{7}, 363-387.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
sa2sp_chem does the reverse, sp2sa_geo
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
