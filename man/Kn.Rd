\encoding{latin1}
\name{Kn}
\alias{Kn}
\title{Dissociation constant of ammonium (mol/kg)}
\description{Dissociation constant of ammonium on the total scale (mol/kg)}
\usage{
Kn(S=35, T=25, P=0, pHscale="T", kSWS2scale="x", warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Optional conversion factor from the seawater pH scale to the chosen pH scale. Supply the value already computed for the other constants in the same calculation (as \code{calculate_carb} does for \code{Kw}, \code{Ksi}, \code{K1p}, \code{K2p} and \code{K3p}), so that every constant is placed on the same pH scale by the same route. If omitted, the factor is recomputed here with \code{\link{kconv}} using its default \code{kf="x"}; a caller who requested \code{kf="dg"} would then silently receive a constant built on a different \code{Kf}, an inconsistency of up to 0.9 percent.}
\item{warn}{"y" to show warnings when T or S go beyond the valid range for Kn; "n" to supress warnings. The default is "y".}
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
	\item{Kn}{Dissociation constant of ammonium (mol/kg)}
}

\references{
Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl} and Heloise Lavigne
}

%\seealso{}

\examples{
  Kn(S=35,T=25,P=0, pHscale="T")
}

\keyword{utilities}
