\encoding{latin1}
\name{Pcorrect}
\alias{Pcorrect}
\title{Pressure correction of equilibrium constants}
\description{Computes the pressure correction of the equilibrium constants}
\usage{
Pcorrect(Kvalue, Ktype, T=25, S=35, P=0, pHscale="T", 
  kconv2ScaleP0=0, kconv2Scale=, warn="y")
}

\arguments{
  \item{Kvalue}{Value of the constant at P=0 (hydrostatic pressure in bar, surface = 0)}
  \item{Ktype}{Name of the constant, 
	\itemize{
	\item K1	First dissociation constant of carbonic acid (mol/kg)
	\item K2	Second dissociation constant of carbonic acid (mol/kg)
	\item Kb	Dissociation constant of boric acid (mol/kg)
	\item Kw	Ion product of water (mol2/kg2)
	\item Ks	Stability constant of hydrogen sulfate (mol/kg)
	\item Kf	Stability constant of hydrogen fluoride (mol/kg)
	\item Kspc	Solubility product of calcite (mol/kg)
	\item Kspa	Solubility product of aragonite (mol/kg)
	\item K1p	First dissociation constant of phosphoric acid (mol/kg)
	\item K2p	Second dissociation constant of phosphoric acid (mol/kg)
	\item K3p	Third dissociation constant of phosphoric acid (mol/kg)
	\item Khs	Dissociation constant of hydrogen sufide(mol/kg)
	\item Kn	Dissociation constant of ammonium (mol/kg)
	\item Ksi	Dissociation constant of Si(OH)4 (mol/kg)
}}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{S}{Salinity, default is 35}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{pHscale of the constant given in Kvalue}
  \item{kconv2ScaleP0}{Conversion factor from the pH scale selected to the SWS (or free for Kf) scale at pressure zero. It is computed if it is not given, which significantly slows down the computation}
  \item{kconv2Scale}{Conversion factor from the pH scale selected to the SWS (or free for Kf) scale at the hydrostatic pressure value indicated. It is computed if it is not given, which significantly slows down the computation}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for constants; "n" to supress warnings. The default is "y".}
}    

\details{
  \itemize{
  \item The pressure correction is applied on the seawater scale for K1, K2, K1p, K2p, K3p, Kb, Khs, Kn, Ksi and Kw. Hence the K value is first converted on the seawater scale if needed. After pressure correction, the constant is converted back to the initial pH scale.
  \item The pressure correction is applied on the free scale for Kf.
  \item There is no issue of pH scale for Ks, Kspa and Kspc.
}
}

\value{
  The equilibrium constant given in argument but after pressure correction
}

\references{
Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica et Cosmochimica Acta} \bold{59} 661-677.
}

\author{
  Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
Pcoeffs
}

\examples{
k10 <- K1(T=25, P=0, S=35)
Pcorrect(Kvalue=k10, Ktype="K1", P=300, T=25, S=35, pHscale="T")
}

\keyword{utilities}
