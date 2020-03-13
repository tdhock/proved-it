marker.dt <- nc::capture_all_str(
  "Identifier_Plus_Bins_v1X.txt",
  nc::field("Marker", " Name\t", ".*"),
  "\n",
  csv=nc::quantifier("[0-9].*\n", "+"))
allele.dt <- marker.dt[, data.table::fread(
  csv, colClasses=list(character=1),
  col.names=c("allele", "bases", "V3", "V4", "virtual")
), by=Marker]
table(allele.dt$V3)
table(allele.dt$V4)
