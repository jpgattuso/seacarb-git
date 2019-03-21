# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
"bor" <-
function(S=35, b="u74"){  #three formulations possible Lee et al. (2010) : "l10"; Uppstrom (1974) : "u74"; "Kulinski et al. (2018)" : "k18"

n <- max(c(length(S), length(b)))
if(length(S)!=n){ S <- rep(S[1],n)}
if(length(b)!=n){ b <- rep(b[1],n)}
bor <- rep(NA, n)
method <- rep(NA, n)

# Lee et al. (2010)
  bor_l10 <- 0.1336*S #boron in mg/kg
  method <- "Lee et al. (2010)"

# Uppstrom (1974)
  bor_u74 <- 0.1284*S #boron in mg/kg
  method <- "Uppstrom (1974)"

# Kulinski et al. (2018)
  bor_k18 <- (11.405*S + 11.869)*10.811*10^-3 #boron in mg/kg
  method <- "Kulinski et al. (2018)"

#--------------- choice between the formulations -------------------
is_l10 <- (b=='l10')
is_u74 <- (b=='u74')
is_k18 <- (b=='k18')
bor[is_l10] <- bor_l10[is_l10]
bor[is_u74] <- bor_u74[is_u74]
bor[is_k18] <- bor_k18[is_k18]

method <- rep(NA, n)
method[is_l10] <- "Lee et al. (2010)"
method[is_u74] <- "Uppstrom (1974)"
method[is_k18] <- "Kulinski et al. (2018)"

# conversion from mg/kg to mol/kg
bor <- bor*10^(-3)/10.811
attr(bor,"unit") <- "mol/kg"
attr(bor, "method") <- method

return(bor)
}
