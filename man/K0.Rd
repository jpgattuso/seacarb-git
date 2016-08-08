\encoding{latin1}
\name{K0}
\alias{K0}
\title{Henry's constant mol/(kg/atm)}
\description{Henry's constant mol/(kg/atm)}
\usage{
K0(S=35, T=25, P=0, Patm=1, warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{Patm}{Surface atmospheric pressure in atm, default is 1 atm}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for K0; "n" to supress warnings. The default is "y".}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 0 and 45 and T ranging between -1 and 45oC.
}

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.

For pressure corrections: the pressure correction term of Weiss (1974) is used.
} 

\value{
  \item{K0}{Henry's constant mol/(kg/atm)}
 
}

\references{Weiss R. F., 1974 Carbon dioxide in water and seawater: the solubility of a non-ideal gas. \emph{Marine Chemistry} \bold{2}, 203-215.}

\author{
Jean-Marie Epitalon, Aurelien Proye, and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
  K0(S=35,T=25,P=0)
}
\keyword{utilities}
