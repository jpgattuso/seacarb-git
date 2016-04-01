#create Pcoeffs data
Pcoeffs <- read.table("data-raw/Pcoeffs.txt", sep="\t", header=TRUE, row.names = 1)
devtools::use_data(Pcoeffs, internal = TRUE)
