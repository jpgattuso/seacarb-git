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
  \item{b}{"l10" for using the formulation of Lee et al. (2010), "u74" for using the Uppstrom (1974), or "k18" for using the Kulinski et al. (2018), default is "u74" }
  }

\details{Note that the formulation of Kulinski et al. (2018) is specifically designed for the Baltic Sea. Three formulations are described in their paper:
\itemize{
  \item based on their measurements: TB = [umol/kg] = 10.838 * S + 13.821
  \item based on Kremling (1970 and 1972): TB [umol/kg] = 11.44 * S + 12.6; R2 = 0.95
  \item consensus regression (Kremling + their data): TB [umol/kg] = 11.405 * S + 11.869; R2 = 0.97s
}
The latter formulation is used here.
}

\value{
  \item{bor}{total boron concentration (\eqn{mol\ kg^{-1}}))}

}
\references{
DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Kulinski K., Szymczycha B., Koziorowska K., Hammer K. & Schneider B., 2018. Anomaly of total boron concentration in the brackish waters of the Baltic Sea and its consequence for the CO2 system calculations. \emph{Marine Chemistry}. doi:s10.1016/j.marchem.2018.05.007.

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
