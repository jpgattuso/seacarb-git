\encoding{latin1}
\name{buffzwg}
\alias{buffzwg}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Bolin sensitivity and Revelle factor of the seawater carbonate
  system from Zeebe and Wolf-Gladrow (2001), corrected and enhanced}

\description{Returns the Bolin sensitivity \eqn{B_e = \partial
  DIC/\partial [CO_2^*]} at constant total alkalinity, temperature, and
  salinity, following the closed-form expression given by Zeebe and
  Wolf-Gladrow (2001, Eq. 1.5.88). Egleston et al., (2010) mentioned
  that their Eq. 1.5.88 contains a typographical error but did not
  elaborate. To identify the error, we rederived the equation, finding
  it to be a missing factor of \eqn{h/s} on one of the two numerator
  terms, as detailed in comments of the \code{buffzwg} function. It was
  later found that this same corrected equation had been given in the errata
  of the Zeebe and Wolf-Gladow book (in an errata document from 2010,
  where it was mentioned that the correction had not yet been
  implemented in more recent editions of the book).  In any case, the
  corrected formula, which in Zeebe and Wolf-Gladrow (2001) already accounted for
  carbonate, borate, and water alkalinity, has been extended here
  to also account for phosphoric and silicic acid systems as well as
  NH3 and H2S acid-base systems, following the
  same approach used to extend seacarb's \code{buffesm} function. The
  classic Revelle factor is also returned, computed as \eqn{RF =
  (DIC/[CO_2^*])/B_e}. This routine was validated against
  \code{buffesm}: when Pt = Sit = 0, and separately when either is
  nonzero, buffzwg's \code{RF} and buffesm's \code{R}
  agree to within floating-point precision. This cross-validation
  also led to the correction of an unrelated bug in
  buffesm's own phosphate term (an R line-continuation parsing
  error in \code{numPt}, fixed in this package version). The
  \code{buffzwg} input arguments are identical to those in the
  \code{buffesm} and \code{carb} functions of seacarb.}

\usage{
buffzwg(flag, var1, var2, S = 35, T = 25, Patm = 1, P = 0, Pt = 0, Sit = 0, 
  NH4t = 0, HSt = 0, k1k2 = "x", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y", 
  eos = "eos80", long = 1e+20, lat = 1e+20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{flag}{select the pair of variables available. The flags which can be used are: 

flag = 1      pH and CO2 given

flag = 2      CO2 and HCO3 given

flag = 3      CO2 and CO3 given

flag = 4      CO2 and ALK given

flag = 5      CO2 and DIC given

flag = 6      pH and HCO3 given

flag = 7      pH and CO3 given

flag = 8      pH and ALK given

flag = 9      pH and DIC given

flag = 10     HCO3 and CO3 given

flag = 11     HCO3 and ALK given

flag = 12     HCO3 and DIC given

flag = 13     CO3 and ALK given

flag = 14     CO3 and DIC given

flag = 15     ALK and DIC given

flag = 21     pCO2 and pH given

flag = 22     pCO2 and HCO3 given

flag = 23     pCO2 and CO3 given

flag = 24     pCO2 and ALK given

flag = 25     pCO2 and DIC given
}
  \item{var1}{enter value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
  \item{var2}{enter value of the second  variable in mol/kg, except for pH}
  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{Patm}{Surface atmospheric pressure in atm, default 1 atm}
  \item{P}{Hydrostatic pressure in bar (surface = 0)}
  \item{Pt}{Concentration of total phosphate in mol/kg; set to 0 if NA; when nonzero, accounts for the phosphoric acid system, unlike in Zeebe and Wolf-Gladrow (2001).}
  \item{Sit}{Concentration of total silicate in mol/kg; set to 0 if NA; when nonzero, accounts for the silicic acid system, unlike in Zeebe and Wolf-Gladrow (2001).}
  \item{NH4t}{Concentration of total ammonium in mol/kg; set to 0 if NA; when nonzero (or when \code{HSt} is nonzero), accounts for ammonia alkalinity and \code{carbfull()} is used internally instead of \code{carb()}}
  \item{HSt}{Concentration of total hydrogen sulfide in mol/kg; set to 0 if NA; when nonzero (or when \code{NH4t} is nonzero), accounts for sulfide alkalinity and \code{carbfull()} is used internally instead of \code{carb()}}
  \item{k1k2}{"cw" for using K1 and K2 from Cai & Wang (1998), "l" from Lueker et al. (2000), "m02" from Millero et al. (2002), "m06" from Millero et al. (2006), "m10" from Millero (2010), "mp2" from Mojica Prieto et al. (2002), "p18" from Papadimitriou et al. (2018), "r" from Roy et al. (1993), "sb21" from Shockman & Byrne (2021), "s20" from Sulpis et al. (2020), and "w14" from Waters et al. (2014). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "w14".}
  \item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979 in Dickson and Goyet, 1994). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
  \item{ks}{"d" for using Ks from Dickon (1990), "k" for using Ks from Khoo et al. (1977), default is "d"} 
  \item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}
  \item{b}{Concentration of total boron. "l10" for the Lee et al. (2010) formulation or "u74" for the Uppstrom (1974) formulation, default is "u74" }
  \item{warn}{"y" to show warnings when T or S go beyond the valid range for constants; "n" to supress warnings. The default is "y".}
  \item{eos}{"teos10" to specify T and S according to Thermodynamic Equation Of Seawater - 2010 (TEOS-10); "eos80" to specify T and S according to EOS-80.}
  \item{long}{longitude of data point, used when eos parameter is "teos10" as a conversion parameter from absolute to practical salinity.}
  \item{lat}{latitude of data point, used when eos parameter is "teos10".}
}


\details{The Lueker et al. (2000) constants for K1 and K2, the Perez and Fraga (1987) constant for Kf and the Dickson (1990) constant for Ks are recommended by Dickson et al. (2007). It is, however, critical to consider that each formulation is only valid for specific ranges of temperature and salinity:

\emph{For K1 and K2:}
\itemize{
\item Cai and Wang (1998): S ranging between 0 and 40 and T ranging between 0.2 and 30oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Millero et al. (2002): S ranging from 34 to 37 and T ranging between -1.6 and 35oC.

\item Millero et al. (2006): S ranging between 0.1 and 50 and T ranging between 1 and 50oC.

\item Millero (2010): S ranging between 1 and 50 and T ranging between 0 and 50oC. Millero (2010) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale". 

\item Mojica Prieto et al. (2002): S ranging from 5 to 42 and T ranging between 0 and 45oC.

\item Papadimitriou et al. (2018): S ranging from 33 to 100 and T ranging between -6 to 25oC.

\item Roy et al. (1993): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Shockman & Byrne (2021): for K2, S ranging from 19.6 to 41 and T ranging between 15 to 35oC. For K1, formulation is that of Waters et al.

\item Sulpis et al. (2020): S ranging from 30.7 to 37.6 and T ranging between -1.7 to 31.8oC.

\item Waters et al.(2014): S ranging between 1 and 50 and T ranging between 0 and 50oC. Waters (2014) provides a K1 and K2 formulation for the seawater, total and free pH scales. Therefore, when this method is used and if P=0, K1 and K2 are computed with the formulation corresponding to the pH scale given in the flag "pHscale".
}

\emph{For Kf:}
\itemize{
\item Perez and Fraga (1987): S ranging between 10 and 40 and T ranging between 9 and 33oC.

\item Dickson and Riley (1979 in Dickson and Goyet, 1994): S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

\emph{For Ks:}
\itemize{
\item Dickson (1990): S ranging between 5 and 45 and T ranging between 0 and 45oC.

\item Khoo et al. (1977): S ranging between 20 and 45 and T ranging between 5 and 40oC.
}

The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It is recommended to use either vectors with the same dimension or one vector for one argument and numbers for the other arguments.

\emph{On the underlying formula:}

Writing \eqn{s := [CO_2^*]} and \eqn{h := [H^+]}, the Bolin sensitivity is computed as

\eqn{B_e = \left( K_1/h + 4 K_1 K_2/h^2 + K_1^2 K_2/h^3 + (h/s) D_s w \right) \Big/ \left( K_1/h + 4 K_1 K_2/h^2 + (h/s) w \right)}{B_e = ( K1/h + 4 K1 K2/h^2 + K1^2 K2/h^3 + (h/s) Ds w ) / ( K1/h + 4 K1 K2/h^2 + (h/s) w )}

where \eqn{D_s = 1 + K_1/h + K_1 K_2/h^2} is the ionization fraction (\eqn{DIC/[CO_2^*]}, following the terminology of Jones et al. (2014) and Egleston et al. (2010)), and \eqn{w} is a non-carbonate buffering term built additively, one acid-base system at a time:

\eqn{w = w_{bw} + w_{Si} + w_P + w_N + w_S}

with \eqn{w_{bw} = 1 + K_w/h^2 + B_T K_B/(K_B+h)^2} (borate and water; this term alone reproduces the original Zeebe and Wolf-Gladrow (2001) formula), \eqn{w_{Si} = Si_T K_{Si}(h^2+4K_{2Si}h+K_{Si}K_{2Si})/(h^2+K_{Si}h+K_{Si}K_{2Si})^2} (silicic acid, treated as the full diprotic system, i.e. \eqn{Si(OH)_4 \rightleftharpoons SiO(OH)_3^- \rightleftharpoons SiO_2(OH)_2^{2-}}, using both \eqn{K_{Si}} and \eqn{K_{2Si}}; this matches seacarb's own \code{carbfull()} treatment of silicate alkalinity), \eqn{w_P} the analogous, algebraically more involved contribution from the three-step phosphoric acid system, \eqn{w_N = NH_{4t} K_N/(K_N+h)^2} (ammonia), and \eqn{w_S = HS_t K_{HS}/(K_{HS}+h)^2} (hydrogen sulfide). Each additional term was verified to equal \eqn{-\partial A_j/\partial h}, where \eqn{A_j(h)} is that system's own contribution to total alkalinity at fixed total concentration, and to be strictly positive over the seawater pH range, as required for a genuine buffering contribution. \eqn{w_N} and \eqn{w_S} are automatically and exactly zero unless \code{NH4t} or \code{HSt}, respectively, are nonzero.

\emph{Pressure corrections and pH scale:}

\itemize{
\item For K0, the pressure correction term of Weiss (1974) is used.

\item For K1, K2, pK1, pK2, pK3, Kw, Kb, Khs and Ksi, the pressure correction was applied on the seawater scale. Hence, if needed, values were first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS).

\item For Kf, the pressure correction was applied on the free scale. The formulation of Dickson and Riley (1979 in Dickson and Goyet, 1994) provides Kf on the free scale but that of Perez and Fraga (1987) provides it on the total scale. Hence, in that case, Kf was first transformed from the total scale to the free scale. With both formulations, the pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 

\item For Ks, the pressure correction was applied on the free scale. The pressure correction was applied as described by Millero (1995), and the value was transformed back to the required scale (T, F or SWS). 
}
long and lat are used as conversion parameters from absolute to practical salinity: when seawater is not of standard composition, practical salinity alone is not sufficient to compute absolute salinity and vice-versa. One needs to know the density. When long and lat are given, density is inferred from WOA silicate concentration at given location. When they are not, an arbitrary geographic point is chosen: mid equatorial Atlantic. Note that this implies an error on computed salinity up to 0.02 g/kg.
}

\value{The function returns a data frame containing the following columns:

  \item{Be}{\eqn{B_e}, the Bolin sensitivity, \eqn{\partial DIC/\partial [CO_2^*]} at constant total alkalinity, temperature, and salinity (unitless, since both DIC and \eqn{[CO_2^*]} are in mol/kg)}
  \item{RF}{the classic Revelle factor, \eqn{(DIC/[CO_2^*])/B_e}, equivalently \eqn{(\partial ln[CO_2]/\partial ln DIC)^{-1}} (unitless); agrees with \code{buffesm}'s \code{R} and \code{buffer}'s \code{BetaD} to floating-point precision}
}

\note{When \code{Sit > 0} and \code{NH4t = HSt = 0}, this function's diprotic treatment of silicate alkalinity (using both \code{Ksi} and \code{K2si}) is very slightly inconsistent with the underlying carbonate system state, because that state comes from \code{carb()}, which internally treats silicate as monoprotic (\code{K2si = 0}) regardless of \code{Sit} (see \code{carb.Rd}). \code{carbfull()} is only invoked here when \code{NH4t} or \code{HSt} are nonzero, so this inconsistency is not corrected for when only \code{Pt}/\code{Sit} are supplied; it is expected to affect results at only the ~1e-6 relative level. See \code{carb.Rd} and \code{carbfull.Rd} for further discussion.}

\references{
Cai W. J., and Wang Y., 1998. The chemistry, fluxes, and sources of carbon dioxide in the estuarine waters of the Satilla and Altamaha Rivers, Georgia. \emph{Limnology and Oceanography} \bold{43}, 657-668.

Dickson A. G., 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Egleston, E. S., Sabine, C. L. and Morel, F. M. M.,  2010 Revelle revisited: Buffer factors that quantify the response of ocean chemistry to changes in DIC and alkalinity, Global Biogeochem. Cycles \bold{24}, GB1002, \doi{10.1029/2008GB003407}. 

Khoo H. K., Ramette R. W., Culberson C. H. and Bates R. G., 1977 Determination of hydrogen ion concentration in seawater from 5 to 40oC: standard potentials at salinities from 20 to 45. \emph{Analytical Chemistry} \bold{49}, 29-34. 

Lee K., Tae-Wook K., Byrne R.H., Millero F.J., Feely R.A. and Liu Y-M, 2010 The universal ratio of the boron to chlorinity for the North Pacific and North Atlantoc oceans. \emph{Geochimica et Cosmochimica Acta} \bold{74} 1801-1811.

Lueker T. J., Dickson A. G. and Keeling C. D., 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero F. J., 2010 Carbonate constant for estuarine waters. \emph{Marine and Freshwater Research} \bold{61}: 139-142.

Millero F. J., Graham T. B., Huang F., Bustos-Serrano H. and Pierrot D., 2006 Dissociation constants of carbonic acid in seawateras a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.  

Perez F. F. and Fraga F., 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy R. N., Roy L. N., Vogel K. M., Porter-Moore C., Pearson T., Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45oC. \emph{Marine Chemistry} \bold{44}, 249-267.

Schockman, K.M., Byrne, R.H., 2021. Spectrophotometric determination of the bicarbonate dissociation constant in seawater, \emph{Geochimica et Cosmochimica Acta}.

Uppstrom L.R., 1974 The boron/chlorinity ratio of the deep-sea water from the Pacific Ocean. \emph{Deep-Sea Research I} \bold{21} 161-162.

Zeebe R. E. and Wolf-Gladrow D. A., 2001 CO2 in seawater: equilibrium, kinetics, isotopes. \emph{Elsevier Oceanography Series} \bold{65}, Amsterdam, 346 pp.
}

\author{
James Orr \email{James.Orr@lsce.ipsl.fr}
}

\examples{

## Computation with a pair of variables
buffzwg(flag=15, var1=0.00230, var2=0.00200, S=35, T=18, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l", b="u74")

## Using vectors as arguments
flag <- c(15, 15, 15)
var1 <- c(0.00230, 0.00230, 0.00230)
var2 <- c(0.00200, 0.00200, 0.00200)
S <- c(35, 35, 30)
T <- c(18, 25, 18)
P <- c(0, 0, 0)
Pt <- c(0, 5e-6, 5e-6)
Sit <- c(0, 0, 20e-6)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
b <- c("l10", "l10", "l10")
buffzwg(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, 
	Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale, b=b)

## Compare against buffesm's classic Revelle factor (should agree to
## floating-point precision, for any Pt and Sit)
e <- buffesm(flag=15, var1=0.00230, var2=0.00200, S=35, T=18, Pt=5e-6, Sit=20e-6)
z <- buffzwg(flag=15, var1=0.00230, var2=0.00200, S=35, T=18, Pt=5e-6, Sit=20e-6)
e$R - z$RF   # should be ~0 (floating-point precision)

}
\keyword{utilities}
