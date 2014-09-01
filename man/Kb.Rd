\encoding{latin1}
\name{Kb}
\alias{Kb}
\title{Dissociation constant of boric acid (mol/kg)}

\description{Dissociation constant of boric acid (mol/kg)}
\usage{
Kb(S = 35, T = 25, P = 0, pHscale="T", kSWS2scale=0, ktotal2SWS_P0=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from SWS pH scale to chosen scale at chosen pressure;
if needed and not given, it is computed; if given, computation speed is increased}
  \item{ktotal2SWS_P0}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. It is not required when the selected pH scale is SWS and the hydrostatic pressure is 0. If it is required and not given, it is computed, which slows down computations.}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 5 and 45 and T ranging between 0 and 45oC.
}

The pressure correction is applied on the seawater scale. Hence, values are first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{Kb}{Dissociation constant of boric acid (mol/kg)}
}

\references{
Dickson A. G., 1990 Thermodynamics of the dissociation of boric acid in synthetic seawater from 273.15 to 318.15 K. \emph{Deep-Sea Research} \bold{375}, 755-766.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.}

\author{
Jean-Marie Epitalon and Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
  Kb(S=35,T=25,P=0,pHscale="T")
}
\keyword{utilities}
