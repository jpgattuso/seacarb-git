\encoding{latin1}
\name{K2si}
\alias{K2si}
\title{Second dissociation constant of Si(OH)4}
\description{Second dissociation constant of Si(OH)4 (mol/kg)}
\usage{
K2si(S=35, T=25, P=0, pHscale="T", kSWS2scale=0, ktotal2SWS_P0=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
  \item{kSWS2scale}{Conversion factor from the seawater scale (SWS) to the pH scale selected at the hydrostatic pressure value indicated. If it is required and not given, it is computed, which slows down computations.}
  \item{ktotal2SWS_P0}{Conversion factor from the total scale to the SWS at an hydrostatic pressure of 0. If it is not given, it is computed, which slows down computations.}
  
}

\details{This equation is modified from Wischmeyer et al. (2003), who fitted the temperature-dependent K2si from Nordstrom et al. (1990) for freshwater to a value of 12.56 for T=25 and an ionic strength of 0.5 mol/kg. The temperature and salinity ranges in which it is valid are not well constrained.

The pressure correction is applied on the seawater scale. Hence, values are first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value is transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
	\item{K2si}{Second dissociation constant of Si(OH)4 (mol/kg)}
}

\references{
Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.

Nordstrom D. K., L. N. Plummer, D. Langmuir, E. Busenberg, H. M. May, B. F. Jones, D. L. Parkhurst, 1990 Revised chemical equilibrium data from major mineral reactions and their limitations. In: Melchior, D.C., R. L. Bassett (Eds.) Chemical Modeling of Aqueous Systems. \emph{IIACS Series} \bold{416}. American Chemical Society, Washington, DC.

Wischmeyer A. G., Y. Del Amo, M. Brzezinski, D. A. Wolf-Gladrow, 1995 Theoretical constraints on the uptake of silicic acid species by marine diatoms. \emph{Marine Chemistry} \bold{82}: 13-29.
}
             
\author{
Mathilde Hagens (\email{M.Hagens@uu.nl})
}

%\seealso{}

\examples{
  K2si(S=35, T=25, P=0, pHscale="T")
}

\keyword{utilities}
