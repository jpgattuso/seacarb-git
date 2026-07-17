## ---------------------------------------------------------------------------
## Common setup for every script and notebook in dev-validation/.
##
## These scripts test the buffer functions BEFORE they are installed, so they
## source them straight from the package's R/ directory rather than relying on
## library(seacarb).  They also need dlnK(), which is internal and unexported.
##
## Default assumes this directory sits inside the seacarb source tree, i.e.
##   <seacarb>/dev-validation/setup.R   and   <seacarb>/R/*.R
## Override with:  Sys.setenv(SEACARB_RDIR = "/path/to/seacarb/R")
## ---------------------------------------------------------------------------
library(seacarb)

RDIR <- Sys.getenv("SEACARB_RDIR", unset = "../R")
if (!dir.exists(RDIR))
  stop("Cannot find seacarb's R/ directory at '", RDIR,
       "'. Set SEACARB_RDIR, e.g. Sys.setenv(SEACARB_RDIR = '~/seacarb-git/R')")

for (f in c("dlnK.R", "buffderiv.R", "buffsun.R", "buffzwg.R", "buffesm.R", "buffer.R"))
  source(file.path(RDIR, f))
