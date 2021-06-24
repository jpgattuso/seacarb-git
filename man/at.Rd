\encoding{latin1}
\name{at}
\alias{at}
\title{Calculates total alkalinity from potentiometric titration data}
\description{Calculates total alkalinity from potentiometric titration data}
\usage{
at(S=35, T=25, CHCl=0.1, HClunit="mol/L", HCltemp=20, CNaCl=0, pHtot="yes", pHTris=NULL, ETris=NULL, weightSample, E, volume)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35. S must be a single value, not a vector.}
  \item{T}{Temperature in degrees Celsius, default is 25oC, can be given as a vector or as a single value.}
  \item{CHCl}{Normality of the acid in mol/kg-H2O, default is 0.1. CHCl must be a single value, not a vector.}
  \item{HClunit}{Unit of the normality of the acid, default is "mol/L". HClunit must be a single value, not a vector.}
  \item{HCltemp}{Temperature of the acid in degrees Celsius, default is 20oC. HCltemp must be a single value, not a vector.}
  \item{CNaCl}{Concentration of NaCl in mol/kg-H2O, default is 0. CNaCl must be a single value, not a vector.}
  \item{pHtot}{Is pH on the total scale?, default is "yes". }
  \item{pHTris}{pH used for the calibration of the electrode with the TRIS buffer. pHTris must be a single value, not a vector.}
  \item{ETris}{Potential used for the calibration of the electrode in mV. ETris must be a single value, not a vector.}
  \item{weightSample}{Weight of the sample in g. weight must be a single value, not a vector.}
  \item{E}{Potential measured during the titration in mV. E must be a vector.}
  \item{volume}{Volume of acid added to the sample in ml. volume must be a vector.}
  }

\details{Total alkalinity is estimated using the non-linear least-square procedure described by Dickson et al. (2007).
} 

\value{A list including;
	\item{AT}{Total alkalinity in mol/kg)}
	\item{E0}{Equivalence point in mV)}
}

\references{Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}
             
\author{
May-Lin Paulsen and Jean-Pierre Gattuso
}

\seealso{alkalinity}

\examples{
data(alkalinity)
d <- alkalinity
at(S=38, T=d$T, CHCl=0.1, HClunit="mol/L", CNaCl=0, HCltemp=20, pHtot="yes", pHTris=8.13, ETris=-71.544, weightSample = 25.753, E = d$E, volume = d$volume)
}

\keyword{utilities}
