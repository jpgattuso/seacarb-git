\encoding{latin1}
\name{at}
\alias{at}
\title{Calculates total alkalinity from potentiometric titration data}
\description{Calculates total alkalinity from potentiometric titration data}
\usage{
at(S=35, T=25, C=0.1, d=1, pHTris=NULL, ETris=NULL, weight, E, volume)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35. S must be a single value, not a vector.}
  \item{T}{Temperature in degrees Celsius, default is 25oC, can be given as a vector or as a single value.}
  \item{C}{Normality of the acid, default is 0.1. C must be a single value, not a vector.}
  \item{d}{Density of the acid, default is 1. d must be a single value, not a vector.}
  \item{pHTris}{pH used for the calibration of the electrode with the TRIS buffer. pHTris must be a single value, not a vector.}
  \item{ETris}{Potential used for the calibration of the electrode in mV. ETris must be a single value, not a vector.}
  \item{weight}{Weight of the sample in g. weight must be a single value, not a vector.}
  \item{E}{Potential measured during the titration in mV. E must be a vector.}
  \item{volume}{Volume of acid added to the sample in ml. volume must be a vector.}
  }

\details{Total alkalinity is estimated using the non-linear least-square procedure described by Dickson et al. (2007).
} 

\value{
	\item{AT}{Total alkalinity in mol/kg)}
}

\references{Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.
}
             
\author{
Steeve Comeau, Heloise Lavigne and Jean-Pierre Gattuso
}

\seealso{
  \code{\link{alkalinity}}.
  }

\examples{
data(alkalinity)
data <- alkalinity
AT <- at(S=data$S[1], T=data$temperature, C=data$normality[1], pHTris=data$pHTris[1],
      ETris=data$ETris[1], E=data$E, weight=data$weight[1], volume=data$volume)
}

\keyword{utilities}
