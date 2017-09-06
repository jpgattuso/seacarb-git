\encoding{latin1}
\name{sf_calib}
\alias{sf_calib}

\title{Calibration coefficients for seaFET sensor}
\description{Calibration coefficients E0INT,25, E0EXT,25}
\usage{
sf_calib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calSal=35.6)

sf_calib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calSal=35.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calEint}{EINT (V), default is 0.0865}
  \item{calEext}{EEXT (V), default is -0.93}
  \item{calpH}{spectrophotometric pH in Total scale, default is 8.132}
  \item{calT}{Temperature in degrees Celsius, default is 16.2}
  \item{calSal}{Salinity, default is 35.6}
}

\details{
Outputs E0INT25 and E0EXT2 must be calculated for each reference sample collected during a SeaFET deployment. Multiple E0INT25 and E0EXT2 may be calculated if there is more than one reference sample for a given deployment. As such, arguments can be given as unique numbers or as vectors (vectors should be of the same length).

It is critical that Eint (calEint) and Eext (calEext) recorded by the SeaFET match reference sample measurements of temperature (calT), salinity (calSal), and spectrophotometric pH (calpH, total hydrogen ion scale) taken at the same time. Note that SeaFET temperature measurements may require calibration via an applied offset. When possible, calibrated CTD temperature and salinity measurements may be used, while spectrophotometric pH measurements always require discrete 'reference' water samples (unless in situ, certified, seawater-based, Tris pH buffer is used). The accepted time offset between collection of reference samples and SeaFET measurements depends on the hydrology and pH variability of the location. For detailed SeaFET calibration instructions and recommendations see Bresnahan et al. (2014) and Rivest et al. (2016).
} 

\value{
This function returns a dataframe comprising 2 variables: 
  \item{E0int25}{Calibration coefficients of the internal sensor at 25oC.}
  \item{E0ext25}{Calibration coefficients of the external sensor at 25oC.}

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
	\code{\link{sf_calc}}.
}

\examples{
sf_calib(calEint=0.0865, calEext=-0.93, calpH=8.132, calT=16.2, calSal=35.6)

## Using the test file seaFET
sf_calib(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  calpH=seaFET$pHspectro, calT=seaFET$Temperature, 
  calSal=seaFET$Salinity)

# There is a test file "seaFET".

sf_calib(calEint=seaFET$Eint, calEext=seaFET$Eext, 
  calpH=seaFET$pHspectro, calT=seaFET$Temperature, 
  calSal=seaFET$Salinity)

}
\keyword{utilities}
