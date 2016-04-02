\encoding{latin1}
\name{Ks}
\alias{Ks}
\title{Stability constant of hydrogen sulfate (mol/kg)}
\description{Stability constant of hydrogen sulfate (mol/kg)}
\usage{
Ks(S = 35, T = 25, P = 0, ks="d")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{ks}{"d" for using Ks from Dickson (1990), "k" for using Ks from Khoo et al. (1977), default is "d"} 
}

\details{The Dickson (1990) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\itemize{
\item Dickson (1990): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Khoo et al. (1977): S ranging between 20 and 45 and T ranging between 5 and 40oC.
}

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length. 

The pressure correction is applied on the free scale as described by Millero (1995), and the value transformed back to the required scale (T, F or SWS).
} 

\value{
  \item{Ks}{Stability constant of hydrogen sulfate (mol/kg), pHscale = free scale}
}

\references{
Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of Hydrogen Ion Concentration in Seawater from 5 to 40oC: Standard Potentials at Salinities from 20 to 45. \emph{Analytical Chemistry} \bold{22}, vol49 29-34. }

\author{
Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
  Ks(S=35,T=25,P=0, ks="d")
}
\keyword{utilities}
