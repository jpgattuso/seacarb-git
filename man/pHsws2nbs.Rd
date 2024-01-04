\encoding{latin1}
\name{pHsws2nbs}
\alias{pHsws2nbs}
\title{Converts pH from SWS to NBS scale}
\description{Converts pH on SWS scale to pH on NBS scale pCO2. 
The SWS-to-NBS conversion is done with the total activity coefficient fH (combined activity coeff for H+, HSO4-, and HF)
from Takahashi et al. (1982) based on data from Culberson and Pytkowicz (1973). The approach is old and full of uncertainty.
Newer approaches are more complicated (Pitzer equations) and big uncertainties remain (Marion et al., 2011; Pilson, 2013).
}
\usage{
pHsws2nbs(pHSWS, S=35, T=25)
}
\arguments{
  \item{pHSWS}{pH on SWS scale}
  \item{S}{Salinity on the practical salinity scale, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  }

\value{
  \item{pHNBS}{pH on NBS scale}
}

\details{
The pHNBS is computed from pHSWS using the total activity coefficient and 
relying on the following basic chemistry equation: \eqn{ ah = fH * hsws}, 
where  \eqn{ah}  is the activity of hydrogen ion,
\eqn{fH}  is the total activity coefficient, and
\eqn{hsws = [H+] + [HSO4-] + [HF]}. In other words, 
\eqn{hsws} is the total hydrogen ion conccentration on the seawater scale.

The two pH scales are defined as \eqn{pHNBS = -log10(ah)} and 
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
	\code{\link{fH}} and \code{\link{pHnbs2sws}}
}

\examples{
## Convert pHSWS to pHNBS
   pHSWS = 8.0                          # pH on the SWS scale
   pHNBS = pHnbs2sws(pHNBS, T=25, S=35) # pH on the NBS scale
   print(pHNBS)
## The result is (give check value here)
}

\keyword{utilities}
