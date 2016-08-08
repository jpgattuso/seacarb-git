\encoding{latin1}
\name{Kw}
\alias{Kw}
\title{Ion product of water (mol2/kg2)}
\description{Ion product of water (mol2/kg2)}
\usage{
Kw(S=35, T=25, P=0, pHscale="T", kSWS2scale=0, warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. If it is required and not given, it is computed, which slows down computations.}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for Kw; "n" to supress warnings. The default is "y".}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

The pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{Kw}{Ion product of water (mol2/kg2)}
}

\references{Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica et Cosmochimica Acta} \bold{59} 661-677.
}

\author{
Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
  Kw(S=35,T=25,P=0,pHscale="T")
}
\keyword{utilities}
