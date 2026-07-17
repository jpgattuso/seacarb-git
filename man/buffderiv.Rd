\encoding{latin1}
\name{buffderiv}
\alias{buffderiv}
\title{Analytic partial derivatives of the Bolin sensitivity Be}
\description{
Computes, in closed form, the four partial derivatives of the Bolin sensitivity
\eqn{B_e = \partial DIC / \partial [CO_2^*]}{Be = d(DIC)/d[CO2*]} (at constant total
alkalinity, temperature and salinity) with respect to temperature, salinity, DIC and
total alkalinity.
}
\usage{
buffderiv(flag, var1, var2, S = 35, T = 25, Patm = 1, Pt = 0, Sit = 0,
          NH4t = 0, HSt = 0, k1k2 = "l", kf = "dg", ks = "d",
          pHscale = "T", b = "u74", warn = "y", npolish = 3)
}
\arguments{
  \item{flag}{select the pair of variables available, as in \code{\link{carb}}}
  \item{var1}{value of the first variable in mol/kg, except for pH, and for pCO2 in \eqn{\mu}atm}
  \item{var2}{value of the second variable in mol/kg, except for pH}
  \item{S}{salinity}
  \item{T}{temperature in degrees Celsius}
  \item{Patm}{surface atmospheric pressure in atm, default is 1 atm}
  \item{Pt}{concentration of total dissolved inorganic phosphorus (mol/kg)}
  \item{Sit}{concentration of total dissolved inorganic silicon (mol/kg)}
  \item{NH4t}{concentration of total dissolved inorganic ammonia (mol/kg)}
  \item{HSt}{concentration of total dissolved inorganic hydrogen sulfide (mol/kg)}
  \item{k1k2}{formulation of K1 and K2. \strong{Must be "l"} (Lueker et al., 2000). See Details.}
  \item{kf}{formulation of Kf. \strong{Must be "dg"} (Dickson and Riley, 1979). See Details.}
  \item{ks}{formulation of Ks. Must be "d" (Dickson, 1990).}
  \item{pHscale}{choice of pH scale. Must be "T" (total scale).}
  \item{b}{formulation of total boron. Must be "u74" (Uppstrom, 1974).}
  \item{warn}{"y" to show warnings when T or S go beyond the valid range}
  \item{npolish}{number of Newton refinement steps applied to [H+] (default 3). See Details.}
}
\details{
Surface only (P = 0).

\strong{Method.} The proton concentration \eqn{h} is fixed implicitly by the total
alkalinity constraint
\deqn{F(h; T,S,C_T,A_T) = A_c(h,C_T) + A_{nc}(h) - A_T = 0}{F(h) = Ac(h,DIC) + Anc(h) - ALK = 0}
and \eqn{B_e}{Be} is then an explicit algebraic function \eqn{G(h,C_T,T,S)}{G(h,DIC,T,S)},
the same one evaluated by \code{\link{buffsun}}. For any \eqn{Y} in
\eqn{\{T,S,C_T,A_T\}}{{T, S, DIC, ALK}}, the implicit function theorem gives
\deqn{\frac{dB_e}{dY} = G_Y + G_h \left(-\frac{F_Y}{F_h}\right)}{dBe/dY = G_Y + G_h * (-F_Y / F_h)}
Every term is closed form. The only nonlinear solve is the single call to \code{\link{carb}}
(or \code{\link{carbfull}}) that supplies \eqn{h}. This makes the routine roughly two orders
of magnitude cheaper than differencing \code{carb()} numerically, and free of step-size error.

\strong{Acid-base systems.} Carbonate, borate, water, phosphate, silicate, fluoride, and
(on the \code{carbfull} path) ammonia and sulfide. Sulfate and the free proton are carried
exactly as in SolveSAPHE. The alkalinity is precisely the one that seacarb itself inverts;
\code{alk_residual} is returned so this can be checked at every point.

\strong{Silicate follows the solver.} When \code{NH4t = HSt = 0} the routine calls
\code{\link{carb}}, whose alkalinity treats silicate as \emph{monoprotic} (K1si only).
Otherwise it calls \code{\link{carbfull}}, which is \emph{diprotic}. K2si is set to zero on
the first path, which reduces the diprotic expressions to the monoprotic ones exactly.

\strong{Constants are not optional.} The analytic dK/dT and dK/dS (in \code{dlnK.R}) are the
exact derivatives of specific formulations, so \code{buffderiv} stops rather than silently
returning derivatives of the wrong function. Two of the guards matter in practice:
\code{k1k2 = "x"} silently switches to Waters et al. (2014) when T < 2 or T > 35 or S < 19 or
S > 43, which is exactly where polar surface waters lie; and \code{kf = "x"} silently switches
between Perez and Fraga and Dickson and Riley by T and S range, changing the SWS-to-total
conversion and hence Kw, Ksi, K1p, K2p and K3p by 0.5 to 0.9 percent. Neither switch warns.

\strong{Newton refinement.} \code{carb()} returns \eqn{h} converged to roughly 1e-10
relative. \eqn{B_e}{Be} is steep in \eqn{h} and its derivatives steeper still, so
\code{buffderiv} applies \code{npolish} Newton steps to \eqn{F(h) = 0} using \eqn{F_h}, which
it needs anyway. Speciation, \eqn{B_e}{Be} and all four derivatives then follow from the same
\eqn{h} to machine precision. Set \code{npolish = 0} to use \code{carb()}'s \eqn{h} untouched.
}
\value{
A data frame containing:
  \item{Be}{Bolin sensitivity, dimensionless}
  \item{dBe_dT}{derivative of Be with respect to temperature, per degree Celsius}
  \item{dBe_dS}{derivative of Be with respect to salinity, per practical salinity unit}
  \item{dBe_dCT}{derivative of Be with respect to DIC, per (mol/kg)}
  \item{dBe_dAT}{derivative of Be with respect to total alkalinity, per (mol/kg)}
  \item{dh_dT}{derivative of [H+] with respect to temperature, (mol/kg) per degree Celsius}
  \item{dh_dS}{derivative of [H+] with respect to salinity, (mol/kg) per psu}
  \item{dh_dCT}{derivative of [H+] with respect to DIC, dimensionless}
  \item{dh_dAT}{derivative of [H+] with respect to total alkalinity, dimensionless}
  \item{alk_residual}{Ac(h) + Anc(h) - ALK (mol/kg). A diagnostic. It must vanish to machine
    precision; if it does not, the alkalinity being differentiated is not the one seacarb
    solved, and every derivative above is suspect.}
}
\author{
James Orr \email{james.orr@lsce.ipsl.fr}
}
\references{
Bolin, B. and Eriksson, E., 1959 Changes in the carbon dioxide content of the atmosphere and
sea due to fossil fuel combustion. \emph{The Atmosphere and the Sea in Motion}, Rossby
Memorial Volume, Rockefeller Institute Press, 130-142.

Munhoven, G., 2013 Mathematics of the total alkalinity-pH equation - pathway to robust and
universal solution algorithms: the SolveSAPHE package v1.0.1. \emph{Geoscientific Model
Development} \bold{6}, 1367-1388.

Orr, J. C., Epitalon, J.-M. and Gattuso, J.-P., 2015 Comparison of ten packages that compute
ocean carbonate chemistry. \emph{Biogeosciences} \bold{12}, 1483-1510.
}
\seealso{
\code{\link{buffsun}}, \code{\link{buffer}}, \code{\link{buffesm}}, \code{\link{buffzwg}},
\code{\link{carb}}, \code{\link{carbfull}}, \code{\link{derivnum}}.
}
\examples{
##  Bolin sensitivity and its four partial derivatives,
##  from total alkalinity (2300 umol/kg) and DIC (2000 umol/kg)
buffderiv(flag = 15, var1 = 2300e-6, var2 = 2000e-6, S = 35, T = 20,
          Pt = 0.2e-6, Sit = 3e-6,
          k1k2 = "l", kf = "dg", ks = "d", pHscale = "T", b = "u74")

##  Vectorised over a field. The single carb() call is the only nonlinear solve.
n   <- 1000
S   <- runif(n, 30, 37)
T   <- runif(n, -1.8, 30)
AT  <- runif(n, 2100, 2400) * 1e-6
CT  <- AT * runif(n, 0.85, 0.97)
d   <- buffderiv(15, AT, CT, S = S, T = T, k1k2 = "l", kf = "dg", ks = "d",
                 pHscale = "T", b = "u74", warn = "n")
summary(d[, c("Be", "dBe_dT", "dBe_dS", "dBe_dCT", "dBe_dAT")])

##  The alkalinity residual must be at machine precision
max(abs(d$alk_residual) / AT)
}
\keyword{utilities}
