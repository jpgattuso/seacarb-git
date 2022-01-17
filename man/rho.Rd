\encoding{latin1}
\name{rho}
\alias{rho}

\title{Density of seawater  (kg/m3)}
\description{Calculates the density of seawater (\eqn{kg\ m^{-3}})}
\usage{
rho(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Practical Salinity (PSS-78), default is 35}
  \item{T}{Temperature in degrees Celsius (ITS-90), default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}


\value{
  \item{rho}{Density of seawater (kg/m3)}

}
\references{Millero F. J. and Poisson A., 1981 International one-atmosphere equation of state of seawater. \emph{Deep-Sea Research} \bold{28A}, 625-629.}

\author{
Aurelien Proye and Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}


\examples{
	rho(35,25,0)
}
\keyword{utilities}
