\encoding{latin1}
\name{pHCalib}
\alias{pHCalib}

\title{Calibration coefficients for seaFET sensor}
\description{Calibration coefficients E0INT,25, E0EXT,25}
\usage{
pHCalib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calsal=35.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calEint}{EINT (V), default is 0.0865}
  \item{calEext}{EEXT (V), default is -0.93}
  \item{calpH}{spectrophotometric pH in Total scale, default is 8.132}
  \item{calT}{Temperature in degrees Celsius, default is 16.2}
  \item{calsal}{Salinity, default is 35.6}
}

\details{

This function is adapted from Bresnahan et al........

Samples for spectrophotometric measurements must be taken at the vicinity of the seaFET sensor and at the same time of a seaFET measurement itself. The accepted offset in the time between the seaFET measurement and the discrete sampling depends on the hydrology of the location. Strong currents could bring new seawater masses with different carbonate chemistry parameters.

Talk about the different periodes of calibration and the fact that we use the mean of E0int25 and E0extt25 per periode of deployement for the next function.... 

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
This function returns a dataframe comprising 2 variables: 
  \item{E0int25}{Calibration coefficients of the internal sensor at 25oC.}
  \item{E0ext25}{Calibration coefficients of the external sensor at 25oC.}

}
\references{
Martz, T.R., Connery, J.G., Johnson, K.S., 2010. Testing the Honeywell Durafet for seawater pH applications. Limnol. Oceanogr. Methods 8, 172–184. doi:10.4319/lom.2010.8.172
 
Bresnahan, P. J. J., Martz, T. R., Takeshita, Y., Johnson, K. S. & LaShomb, M.,2014. Best practices for autonomous measurement of seawater pH with the Honeywell Durafet. Methods Oceangr. 9, 44–60.
}

\author{
Samir Alliouane, Lydia Kapsenberg, Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{sf_calc}}.
}

\examples{
pHCalib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calsal=35.6)
# There is a test file "seaFET".

pHCalib(calEint=seaFET$Eint, calEext=seaFET$Eext, calpH=seaFET$pHspectro, calT=seaFET$Temperature, calsal=seaFET$Salinity)
}
\keyword{utilities}
