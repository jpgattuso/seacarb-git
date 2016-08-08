\encoding{latin1}
\name{vapress}
\alias{vapress}

\title{Computes vapor pressure of seawater}
\description{Computes vapor pressure of seawater (atm) from temperature and salinity}
\usage{
vapress(S=35, T=25, form="d2007")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{form}{choose either "d2007" for the best-practices formulation to compute vapor pressure of seawater from Dickson et al. (2007) or "wp1980" for the formulation from weiss and Price (1980).}
  }

\value{
  \item{vapress}{Vapor pressure of seawater in atm}
}

\details{
Computes the vapor pressure of seawater pH20 following best practices
(Dickson et al., 2007). That computed pH20 is identical, when rounded
to the 4th decimal place, with that computed by the equation from
Weiss and Price (1980).}

\references{
Dickson A. G., Sabine C. L. and Christian J. R. (2007) Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Weiss, R. F. (1974) Carbon dioxide in water and seawater: the solubility of a non-ideal gas, Marine Chemistry, \bold{2}, 203-215.

Weiss, R. F. and Price, B. A. (1980) Nitrous oxide solubility in water and seawater, Marine Chemistry, \bold{8}, 347-359.
}
             
\author{
James Orr \email{james.orr@lsce.ipsl.fr}
}

\seealso{
	\code{\link{x2pCO2}}, and \code{\link{p2xCO2}}.
}

\examples{
   pH20 <- vapress(S=35, T=25, form="d2007")
}

\keyword{utilities}
