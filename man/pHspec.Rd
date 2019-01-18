\encoding{latin1}
\name{pHspec}
\alias{pHspec}

\title{Calculates pHT from results of spectrophotometric measurements}
\description{Calculates pH of a water sample from absorbance ratios R, obtained from spectrophotometric measurements with pH indicator dyes (on the total scale in mol/kg-soln)}
\usage{
pHspec(S=35, T=25, R=1, d="mCP", k="m18", warn="y")
%- maybe also 'usage' for other objects documented here.
}

\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{R}{Absorbance ratio, default is 1}
  \item{d}{Dye used for spectrophotometric measurement, default is "mCP"}
  \item{k}{"m18" for using mCP characterization by Mueller and Rehder (2018)} 
  \item{warn}{"y" to show warnings when S and/or T go beyond the valid range for the chosen d and k; "n" to supress warnings. The default is "y".}
}

\details{The model used to calculate the return value of this function is based on experimental data. It is critical to consider that the formulation refers to the conditions studied during the characterization experiment and is only valid for the studied ranges of temperature and salinity:
\itemize{
\item Mueller and Rehder (2018): S ranging between 0 and 40, T ranging between 5 and 35oC, and the dye used being m-Cresol purple (mCP) with R referring to the ratio of absorbances at wavelengths 578 and 434 nm.

}
Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{pHspec}{The function returns the pH value of a water sample from absorbance ratios R, obtained from spectrophotometric measurements with pH indicator dyes (on the total scale in mol/kg-soln)}

}
\references{
Mueller, J. D. and Rehder, G., 2018 Metrology for pH Measurements in Brackish Waters - Part 2: Experimental Characterization of Purified meta-Cresol Purple for Spectrophotometric pHT Measurements. \emph{Frontiers in Marine Science} \bold{5:177}, 1-9. https://doi.org/10.3389/fmars.2018.00177


}

\author{
Jens Daniel Mueller \email{jens.mueller@io-warnemuende.de}

Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{amp}}, \code{\link{pHslope}}, \code{\link{pH}}. \code{\link{tris}},
}

\examples{
	##Example should give test value pHT = 7.6713
	 pHspec(S=35, T=25, R=1, d="mCP", k="m18", warn="y")
}
\keyword{utilities}
