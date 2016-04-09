\name{Pcoeffs}

\docType{data}

\alias{Pcoeffs}

\title{Coefficients used for pressure-correcting the equilibrium constants}

\description{
Pressure corrections are based on the following equations: 

\deqn{
  \ln{\frac{K_{i}^{P}}{K_{i}^{0}}} = -\frac{\Delta V_{i}}{RT}.P + 0.5\frac{\Delta K_{i}}{RT}.P^{2} 
 }

with 
\deqn{
   \Delta V_{i} = a_{0} + a_{1}T + a_{2}T^{2}
 }
and 
\deqn{
   \Delta K_{i} = b_{0} + b_{1}T + b_{2}T^{2}
}


The variables are:

\itemize{
	\item K indicating the type of equilibrium constant
	\item coefficient \eqn{a_0}
	\item coefficient \eqn{a_1}
	\item coefficient \eqn{a_2}
	\item coefficient \eqn{b_0}
	\item coefficient \eqn{b_1}
	\item coefficient \eqn{b_2}

	}
}
\usage{Pcoeffs}

\format{A data frame with 14 rows and 7 variables}

\source{
Millero F. J., 1979 The thermodynamics of the carbonate system in seawater. \emph{Geochemica et Cosmochemica Acta} \bold{43}: 1651-1661.

Millero F. J., 1983 Influence of pressure on chemical processes in the sea. pp. 1-88. In J. P. Riley and R. Chester (eds.), Chemical Oceanography. Academic Press, New York.

Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} \bold{59}: 661-677.
}

\seealso{
Pcorrect
}

\details{
For Kb, to be consistent with Millero (1979) a2 was changed to -2.608e-3 instead of 2.608e-3 (value given in Millero, 1995)
For Kw, coefficients are from Millero (1983).
}

\keyword{datasets}
