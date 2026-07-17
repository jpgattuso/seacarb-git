\encoding{latin1}
\name{K2p}
\alias{K2p}
\title{Second dissociation constant of phosphoric acid (mol/kg)}

\description{Second dissociation constant of phosphoric acid (mol/kg)}
\usage{
K2p(S=35, T=25, P=0, pHscale="T", kSWS2scale="x", warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. It is advised to use default value "x", in which case it is computed when required.}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for K2p; "n" to supress warnings. The default is "y".}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

The pressure correction was applied on the seawater scale. Hence, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{K2p}{Second dissociation constant of phosphoric acid (mol/kg)}

}
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.}

\author{
Jean-Marie Epitalon, Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{K1p}}, \code{\link{K3p}}.
}


\examples{
  K2p(35,25,0)
}
\keyword{utilities}
