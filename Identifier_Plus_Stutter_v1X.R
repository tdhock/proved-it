num.pattern <- list("[0-9.]+", as.numeric)
stutter.dt <- nc::capture_all_str(
  "Identifier_Plus_Stutter_v1X.txt",
  nc::field("Marker", " Name\t", ".*"),
  "\n",
  ratio=num.pattern,
  "\t",
  from.distance=num.pattern,
  "\t",
  to.distance=num.pattern)
