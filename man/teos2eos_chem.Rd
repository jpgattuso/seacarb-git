\encoding{latin1}
\name{teos2eos_chem}
\alias{teos2eos_chem}

\title{Convert temperature and salinity from TEOS-10 to EOS-80}
\description{Converts conservative temperature to in situ temperature and absolute salinity to practical salinity (SP). Salinity conversion depends on total alkalinity as well as on the concentrations of dissolved inorganic carbon, nitrate and silicate.}

\usage{
teos2eos_chem(SA, CT, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{SA}{Absolute salinity in g/kg}
  \item{CT}{Conservative temperature in degrees C}
  \item{P}{Sea water pressure in dbar}
  \item{TA}{Total alkalinity, in mol/kg, default is 2300 µmol/kg}
  \item{DIC}{Dissolved inorganic carbon concentration in mol/kg, default is 2000 µmol/kg}
  \item{NO3}{Total nitrate concentration in mol/kg, default is 0}
  \item{SIOH4}{Total silicate concentration in mol/kg, default is 0}
  }

\value{The function returns a data frame containing the following columns:
  \item{T}{In situ temperature (deg C)}
  \item{SP}{Practical salinity (psu)}
}

\details{
Conversion from absolute to practical salinity depends on carbonate system parameters and ion concentration which mostly affect water density anomalies.
}

\references{
McDougall T. J., Jackett D. R., Millero F. J., Pawlowicz R. and Barker P. M., 2012. A global algorithm for estimating Absolute Salinity. \emph{Ocean Science} \bold{8}, 1123-1134.

Pawlowicz R., Wright D. G. and Millero F. J., 2011. The effects of biogeochemical processes on oceanic conductivity/salinity/density relationships and the characterization of real seawater. \emph{Ocean Science} \bold{7}, 363-387.

Pawlowicz R., 2013. What every oceanographer needs to know about TEOS-10 (The TEOS-10 Primer). http://www.teos-10.org/
}
             
\seealso{
eos2teos_chem, teos2eos_geo, sa2sp_cham, package gsw
}

\author{
Jean-Marie Epitalon
}

\examples{
   # Calculate in situ temperature and practical salinity of a sample with 
   # Absolute salinity of 35 g/kg, Conservative temperature of 18 deg C,
   # at 0 dbar and Total alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
   f <- teos2eos_chem(SA=35, CT=18, P=0, TA=0.00234, DIC=0.00202)
   T <- f$T     # insitu temperature
   SP <- f$SP     # Practical salinity
}

\keyword{utilities}
