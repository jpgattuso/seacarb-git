#create Pcoeffs data
Pcoeffs <- read.table("data-raw/Pcoeffs.txt", sep="\t", header=TRUE)
devtools::use_data(Pcoeffs, internal = TRUE, overwrite = TRUE)
