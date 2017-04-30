\encoding{latin1}
\name{eos2teos_chem}
\alias{eos2teos_chem}

\title{Convert temperature and salinity from EOS-80 to TEOS-10}
\description{Converts in situ temperature to conservative temperature and practical to absolute salinity (SA). Salinity conversion depends on total alkalinity as well as the concentrations of dissolved inorganic carbon, nitrate and silicate.}
\usage{
eos2teos_chem(SP, T, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{Practical salinity on the practical salinity scale}
  \item{T}{In situ temperature in deg. C}
  \item{P}{Sea water pressure in dbar}
  \item{TA}{Total alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved inorganic carbon concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total silicate concentration in mol/kg, default is 0}
  }

\value{The function returns a data frame containing the following columns:
  \item{CT}{Conservative temperature (deg C)}
  \item{SA}{Absolute salinity (g/kg)}
}

\details{
Conversion from practical to absolute salinity depends on carbonate system parameters and ion concentration which mostly affect water density anomalies.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., Wright D. G. and Millero F. J., 2011. The effects of biogeochemical processes on oceanic conductivity/salinity/density relationships and the characterization of real seawater. \emph{Ocean Science} \bold{7}, 363-387.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
teos2eos_chem does the reverse, eos2teos_geo, sp2sa_chem

package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate Conservative Temperature and Absolute Salinity of a sample with 
   # Practical salinity of 35 psu, in-situ temperature of 18 deg C,
   # at 0 dbar and total alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- eos2teos_chem(SP=35, T=18, P=0, TA=0.00234, DIC=0.00202)
   CT <- f$CT     # Conservative Temperature
   SA <- f$SA     # Absolute Salinity
}

\keyword{utilities}
