\encoding{latin1}
\name{fH}
\alias{fH}
\title{Total activity coefficient for H+}
\description{Compute total hydrogen ion activity coefficient.
The activity coefficient (fH) is used to convert from H+ conccentration on SWS scale 
to H+ activity (ah), as used for NBS scale and vice versa. Here, fH is taken from Takahashi 
et al (1982, GEOSECS Pacific Expedition, Chap 3, p. 80) who say:
fH is the total activity coeff., which includes contributions from HSO4- and HF [as well as H+].

Takahashi et al. (1982) computed a relationship for fH based on experimental data
from Culberson & Pytkowicz (1973), who determined it experimentally as a function
of temperature and salinity.  The approach is old and full of uncertainty.
Newer approaches are more complicated (Pitzer equations) 
and big uncertainties remain (Marion et al., 2011; Pilson, 2013).
}
\usage{
fH(S=35, T=25)
}
\arguments{
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  }

\value{
  \item{fH}{Total activity coefficient for H+}
}

\details{
This total activity coefficient appears in
the following basic chemistry equation: \eqn{ ah = fH * hsws}, 
where  \eqn{ah}  is the activity of hydrogen ion,
\eqn{fH}  is the total activity coefficient, and
\eqn{hsws = [H+] + [HSO4-] + [HF]}. In other words, 
\eqn{hsws} is the total hydrogen ion conccentration on the seawater scale.

The two pH scales of concern are defined as \eqn{pHNBS = -log10(ah)} and 
\eqn{pHSWS = -log10(hsws)}.  
}

\references{
Culberson, C.H., & Pytkowicz, R.M. (1973). Ionization of water in seawater. 
\emph{Marine Chemistry}, \bold{1(4)}, 309-316.

Marion G.M., Millero F.J., Camoes M.F., Spitzer P., Feistel R., Chen C.T.A. 2011. pH of seawater. 
\emph{Marine Chemistry} \bold{126} 89-96.

Pilson M.E.Q. (2013) An introduction to the chemistry of the sea, 
2 edn. Cambridge, UK: Cambridge University Press.

Takahashi T., Williams R.T., and Ros D.L. (1982) Carbonate chemistry. 
GEOSECS Pacific Expedition, Volume 3, Hydrographic Data 1973-1974, 77-83.
}
             
\author{
James Orr \email{james.orr@lsce.ipsl.fr}
}

\seealso{
	\code{\link{pHnbs2sws}} and \code{\link{pHsws2nbs}}
}

\examples{
## Compute fH
   f = fH(T=25, S=35)
   print(f)
## The result is (give check value here)
}

\keyword{utilities}
