\encoding{latin1}
\name{fCO2insi}
\alias{fCO2insi}
\title{fCO2 at in situ temperature}
\description{Correction to compensate for the difference in temperature between the temperature of measurement and in situ temperature.}
\usage{
fCO2insi(fCO2lab = 400, Tlab = 20, SST = 19)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fCO2lab}{Fugacity of CO2 measured in the lab in \eqn{\mu}atm}
  \item{Tlab}{Temperaure of measurment in the lab in in degrees Celsius, default is 20oC}
  \item{SST}{Temperature in degrees Celsius, default is 19oC}
  }

\value{
  \item{fCO2insi}{Fugacity of CO2 at in situ temperature in \eqn{\mu}atm, the same units as that for the fCO2 input.}
}

\note{The empirical correction applied is from Takahashi (1993) as recommended by Pierrot et al. (2009)}

\references{
Pierrot D., Neill C., Sullivan K., Castle R., Wanninkhof R., Lüger H., Johannessen T., Olsen A., Feely R. A. and Cosca C. E., 2009. Recommendations for autonomous underway pCO2 measuring systems and data-reduction routines. \emph{Deep-Sea Res. II} \bold{56}, 512-522.

Takahashi T., Olafsson J., Goddard J. G., Chipman D. W. and Sutherland S. C., 1993. Seasonal variation of CO2 and nutrients in the high-latitude surface oceans: a comparative study. \emph{Glob. Biogeochem. Cycles} \bold{7}, 843-878.
}
             
\author{
Jean-Pierre Gattuso \email{jean-pierre.gattuso@imev-mer.fr}
}

\seealso{
	\code{\link{f2pCO2}}, \code{\link{pCO2insi}}.
	}

\examples{
fCO2insi(fCO2lab = 400, SST = 15, Tlab = 16)
}

\keyword{utilities}
