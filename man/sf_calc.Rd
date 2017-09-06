\encoding{latin1}
\name{sf_calc}
\alias{sf_calc}

\title{Calculation of calibrated pH for seaFET sensor}

\description{
The function sf_calc() calculates pH time series (pHint_tot and pHext_tot) for SeaFET pH sensors, using calibration coefficients E0int25 and E0ext25 from the function sf_calib(). Both functions are R-adaptations from MATLAB scripts published by Bresnahan et al. (2014).
}

\usage{
sf_calc(calEint=0.0865, calEext= -0.93, E0int25 =-0.39, 
  E0ext25=-1.46, calT=16.2, calSal=35.6)
sf_calc(calEint=0.0865, calEext= -0.93, E0int25 =-0.39, 
  E0ext25=-1.46, calT=16.2, calSal=35.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calEint}{EINT (V), default is 0.0865}
  \item{calEext}{EEXT (V), default is -0.93}
  \item{E0int25}{Coeficient of calibration related to the internal sensor and obtained via sf_calib function, default is -0.39. If time-serie, we use the mean per periode of deployement}
  \item{E0ext25}{Coeficient of calibration related to the external sensor and obtained via sf_calib function, default is -1.46. If time-serie, we use the mean per periode of deployement}
  \item{calT}{Temperature in degrees Celsius, default is 16.2}
  \item{calSal}{Salinity, default is 35.6}
}

\details{
Input values should be vectors of equal length. E0int25 and E0ext25 should be constant throughout the time series. When multiple reference samples are available for one SeaFET deployment, mean E0int25 and mean E0ext25 should be calculated and used in sf_calc(). Each unique SeaFET deployment requires a new calculation of mean E0int25 and mean E0ext25 based on reference pH samples (total hydrogen ion scale). For detailed SeaFET calibration instructions and recommendations see Bresnahan et al. (2014) and Rivest et al. (2016).
} 

\value{
This function returns a dataframe comprising 2 variables: 
  \item{pHint_tot}{Calibrated pH of the internal sensor at in situ temperature.}
  \item{pHext_tot}{Calibrated pH of the external sensor at in situ temperature.}

}
\references{
Bresnahan, P. J. J., Martz, T. R., Takeshita, Y., Johnson, K. S., and LaShomb, M., 2014. Best practices for autonomous measurement of seawater pH with the Honeywell Durafet. \emph{Methods Oceanogr.} \bold{9}, 44–60.

Martz, T.R., Connery, J.G., and Johnson, K.S., 2010. Testing the Honeywell Durafet for seawater pH applications. \emph{Limnol. Oceanogr. Meth.} \bold{8}, 172–184.

Rivest, E.B., O'Brien, M., Kapsenberg, L., Gotschalk, C.C., Blanchette, C.A., Hoshijima, U., and Hofmann, G.E., 2016. Beyond the benchtop and the benthos: dataset management planning and design for time series of ocean carbonate chemistry associated with Durafet®-based pH sensors. \emph{Ecological Informatics} \bold{36}, 209–220.
}

\author{
Samir Alliouane, Lydia Kapsenberg, Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{sf_calib}}.
}

\examples{
sf_calc(calEint=0.0865, calEext= -0.93, E0int25 =-0.39, 
  E0ext25=-1.46, calT=16.2, calSal=35.6)

## Using the test file seaFET
sf_calc(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  E0int25=seaFET$E0int25, E0ext25=seaFET$E0ext25, 
  calT=seaFET$Temperature, calSal=seaFET$Salinity)

# There is a test file "seaFET".
sf_calc(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  E0int25=seaFET$E0int25, E0ext25=seaFET$E0ext25, 
  calT=seaFET$Temperature, calSal=seaFET$Salinity)
}
\keyword{utilities}
