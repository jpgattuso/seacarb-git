\encoding{latin1}
\name{Ksi}
\alias{Ksi}
\title{Dissociation constant of Si(OH)4}
\description{Dissociation constant of Si(OH)4 on total scale (mol/kg)}
\usage{
Ksi(S=35, T=25, P=0, pHscale="T", kSWS2scale=0, warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. If it is required and not given, it is computed, which slows down computations.}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for Ksi; "n" to supress warnings. The default is "y".}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

The pressure correction is applied on the seawater scale. Hence, values are first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value is transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
	\item{Ksi}{Dissociation constant of Si(OH)4 (mol/kg)}
}

\references{
DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl} and Heloise Lavigne
}

%\seealso{}

\examples{
  Ksi(S=35, T=25, P=0, pHscale="T")
}

\keyword{utilities}
