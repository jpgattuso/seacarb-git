\encoding{latin1}
\name{tris}
\alias{tris}

\title{pH of TRIS buffer}
\description{Calculates the pH value of TRIS buffered artificial seawater solutions (on the total scale in mol/kg-soln)}
\usage{
tris(S=35,T=25,b=0.04,k="d98",warn="y")
%- maybe also 'usage' for other objects documented here.
}

\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{b}{Molality if TRIS/TRISH+ in moles per kg of water, default is 0.04 mol/kg-H20}
  \item{k}{"d98" for DelValls and Dickson 1998, "m18" for using tris characterization by Mueller et al (2018), default is "d98"} 
  \item{warn}{"y" to show warnings when S,T and/or b go beyond the valid range for the chosen k; "n" to supress warnings. The default is "y".}
}

\details{The models used to calculate the return value of this function are based on experimental data. It is critical to consider that each formulation refers to the artificial seawater solution applied during the characterization experiment and is only valid for the studied ranges of temperature and salinity:
\itemize{
  \item DelValls and Dickson (1998): S ranging between 20 and 40, T ranging between 0 and 45oC, and b being 0.04 mol/kg-H20.
  \item Mueller et al. (2018): S ranging between 5 and 40, T ranging between 5 and 45oC, and b ranging between 0.01 and 0.04 mol/kg-H20.
}
Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
  \item{tris}{The function returns the pH value of TRIS buffered artificial seawater solutions (on the total scale in mol/kg-soln)}

}
\references{
DelValls, T. A., and Dickson, A. G., 1998 The pH of buffers based on 2-amino-2-hydroxymethyl-1,3-propanediol ('tris') in synthetic sea water. \emph{Deep Sea Research Part I: Oceanographic Research Papers} \bold{45(9)}, 1541-1554. https://doi.org/10.1016/S0967-0637(98)00019-3

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Mueller, J. D., Bastkowski, F., Sander, B., Seitz, S., Turner, D. R., Dickson, A. G., and Rehder, G., 2018 Metrology for pH measurements in brackish waters- Part 1: Extending electrochemical pHT measurements of TRIS buffers to salinities 5-20. \emph{Frontiers in Marine Science} \bold{5:176}, 1-12. https://doi.org/10.3389/fmars.2018.00176
}

\author{
Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr} and Jens Daniel Mueller \email{jens.mueller@io-warnemuende.de} 
}

\seealso{
	\code{\link{amp}}, \code{\link{pHslope}}, \code{\link{pH}}.
}

\examples{
	##Example from Mueller et al. (2018), should give test value pHT = 8.0703
	 tris(S=20,T=25,b=0.04,k="m18")
}
\keyword{utilities}
