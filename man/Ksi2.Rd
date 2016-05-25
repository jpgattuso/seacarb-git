\encoding{latin1}
\name{Ksi2}
\alias{Ksi2}
\title{Second dissociation constant of Si(OH)4}
\description{Second dissociation constant of Si(OH)4 on total scale (mol/kg)}
\usage{
Ksi2(S=35, T=25, P=0, pHscale="T", kSWS2scale=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. If it is required and not given, it is computed, which slows down computations.}
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
	\item{Ksi2}{Second dissociation constant of Si(OH)4 (mol/kg)}
}

\references{
Hofmann A. F., Soetaert K., Middelburg J. J. and Meysman F. J. R., 2010. AquaEnv: an aquatic acid£base modelling environment in R. \emph{Aquatic Geochemistry} \bold{16}: 507-546.
}
             
\author{
Mathilde Hagens \email{M.Hagens@uu.nl}
}

%\seealso{}

\examples{
  Ksi2(S=35, T=25, P=0, pHscale="T")
}

\keyword{utilities}
