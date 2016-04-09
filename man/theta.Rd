\encoding{latin1}
\name{theta}
\alias{theta}

\title{Potential temperature of seawater}
\description{Computes theta, the potential temperature of seawater given
original temperature, salinity, pressure, and reference pressure}
\usage{
theta(S=35, T=25, P=0, Pref=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0; 1000 db = 100 bar), default is 0}
  \item{Pref}{Reference hydrostatic pressure in bar, default is 0}
  }

\value{
  \item{theta}{potential temperature of seawater (C)}
}

\details{
Computes the potential temperature of seawater relative to a
chosen reference pressure following Fofonoff and Millard (1983).
The potential temperature \eqn{\theta}{theta} is the temperature that a
water parcel would have if were moved adiabatically to another pressure
level \code{Pref}.  Typically, the potential temperature is referenced
to the surface (\eqn{Pref=0}).  The potential teperature depends on
the original salinity \code{S}, \emph{in-situ} temperature \code{T}
and pressure \code{P}.

This routine is essentially a wrapper for the \code{\link[oce]{swTheta}}
routine of the 'oce' package. Unlike the latter, pressure units here
are given in bars and method="unesco" is prescribed.
}

\references{
Fofonoff, P. and R. C. Millard Jr, 1983. Algorithms for computation of
fundamental properties of seawater. \emph{Unesco Technical Papers in
Marine Science}, \bold{44}, 53 pp.}
             
\author{
James Orr \email{james.orr@lsce.ipsl.fr}
}

\seealso{
	\code{\link[oce]{swTheta}}.
}

\examples{
   #Calculate the potential temperature for a sample at 1000 db referenced to the surface
   theta <- theta(S=35, T=25, P=100, Pref=0)
}

\keyword{utilities}
