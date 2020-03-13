panels.dt <- data.table::fread(
  "Identifier_Plus_Panels_v1X.txt", sep="\t", skip=4)
some.names <- c("Marker", "color", "min.bases", "max.bases")
names(panels.dt)[1:length(some.names)] <- some.names

panels.dt[2][allele.dt, on="Marker", nomatch=0L]
## column 5 seems to be the number of alleles

## column 6 is maybe the spacing between alleles in bases?

## column 7 not sure

## column 8 are the names of the non-virtual alleles.
