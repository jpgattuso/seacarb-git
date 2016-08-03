\encoding{latin1}
\name{Kspc}
\alias{Kspc}

\title{Solubility product of calcite (mol/kg)}
\description{Solubility product of calcite (mol/kg)}
\usage{
Kspc(S=35, T=25, P=0, warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for Kspc; "n" to supress warnings. The default is "y".}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 5 and 44 and T ranging between 5 and 40oC.
}

The pressure coorection was performed as described by Millero (1995).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{Kspc}{Solubility product of calcite (mol2/kg)}

}
\references{
Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica et Cosmochimica Acta} \bold{59} 661-677.

Mucci A., 1983 The solubility of calcite and aragonite in seawater at various salinities, temperature, and one atmosphere total pressure. \emph{American Journal of Science} \bold{283}: 780-799.
}

\author{
Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{Kspa}}.
}

\examples{
	Kspc(S=35,T=25,P=0)
}
\keyword{utilities}
