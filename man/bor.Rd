\encoding{latin1}
\name{bor}
\alias{bor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Total boron concentration (mol/kg)}
\description{total boron concentration (\eqn{mol\ kg^{-1}})}
\usage{
bor(S, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{b}{"l10" for using the Lee et al. (2010) formulation or "u74" for using the Uppstrom (1974) formulation, default is "u74" }
  }


\value{
  \item{bor}{total boron concentration (\eqn{mol\ kg^{-1}}))}

}
\references{
DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.
}

\author{
Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
bor(35, "l10")
}

\keyword{utilities}
