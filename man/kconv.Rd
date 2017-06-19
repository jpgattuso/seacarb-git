\encoding{latin1}
\name{kconv}
\alias{kconv}
\title{Conversion factors to change the pH scale of dissociation constants}
\description{Conversion factors from the total scale to the free and seawater scales}
\usage{
kconv(S=35, T=25, P=0, kf, Ks, Kff, warn="y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994);
if the fonction Kf was used previously, the default value is the value given for the argument kf in the fonction Kf. If the Kf function was not used previously, the default value is "pf", except if T is outside the range 9 to 33oC or of S is outside the range 10 to 40. In these cases, the default is "dg".}	
  \item{Ks}{Stability constant of hydrogen sulfate (mol/kg) at given S, T and P, optional;
if not given, it will be computed, if given, it allows for speed optimisation}
  \item{Kff}{Stability constant of hydrogen fluoride (mol/kg) on free pH scale at given S, T and P, optional;
if not given, it will be computed, if given, it allows for speed optimisation and kf parameter is then ignored}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for constants; "n" to supress warnings. The default is "y".}
}

\value{
The function returns a list with 6 conversion factors :
  \item{ktotal2SWS}{to convert from the total scale to seawater scale}
  \item{ktotal2free}{to convert from the total scale to the free scale}
  \item{kfree2SWS}{to convert from the free scale to the seawater scale}  
  \item{kfree2total}{to convert from the free scale to total scale}
  \item{kSWS2total}{to convert from the seawater scale to the total scale}
  \item{kSWS2free}{to convert from the seawater scale to the free scale}  
}

\details{It is critical to consider that each formulation is valid in specific ranges of temperature and salinity:

\itemize{
  \item Perez and Fraga (1987): S ranging between 10 and 40 and T ranging between 9 and 33oC.
  \item Dickson and Riley (1979 in Dickson and Goyet, 1994): S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

Note that kconv may be called in many functions (i.e. K1, K2, K1p, K2p, K3p, Kw, Ksi, K2si, etc...) without user controls it. 
To force a particular formulation for Kf, it is recommended to call kconv() first then pass the resulting conversion factors
to these functions.
}

\references{Dickson A.G. and F.J. Millero, 1987 A comparison of the equilibrium constants for the dissociation of carbonic acid in seawater media. \emph{Deep-Sea Research} \bold{34}:1733-1743.}

             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

\seealso{
	\code{\link{pHconv}}.
}


\examples{
##To convert dissociation constants from the total scale to the free scale
## (at salinity=35, temperature=25oC and atmospheric pressure):
kconv(35,25,0)
conv <- kconv()
 c(K1_total=K1(),K1_SWS=K1()*conv$ktotal2SWS,K1_free=K1()*conv$ktotal2free)  
}

\keyword{utilities}
