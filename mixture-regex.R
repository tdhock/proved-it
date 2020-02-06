subject.vec <- c(
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051316LEA_5sec)/H06_RD14-0003-49_50_29-1;4;1-M4I35-0.09IP-QLAND_004.5sec.fsa",
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051616LEA_5sec)/H07_RD14-0003-35_36_37_38_39-1;4;4;4;1-M2I35-0.75IP-QLAND_004.5sec.fsa")
## Comments explaining the regex copied from
## https://lftdi.camden.rutgers.edu/wp-content/uploads/2019/12/PROVEDIt-Database-Naming-Convention-Laboratory-Methodsv1.pdf
file.dt <- nc::capture_first_vec(
  subject.vec,
  "/",
  letter="[A-Z]",
  numbers="[0-9]+",
  "_",
  project.number=list(
    ".*?",
    "-",
    "[0-9]+"),
  "-",
  ## It is to be noted that each sample is designated by the
  ## combination of project number and sample identifier. For example,
  ## RD14-0003-21 will have a different known genotype than
  ## RD12-0002-21.
  identifiers.str="[^-]+", # within that project.
  "-",
  ratios.str="[^-]+",
  "-",
  dilution="[^-]+",
  "-",
  template.nanograms="[0-9.]+", as.numeric,
  ## IP is the amplification kit type (i.e., IP=Identifiler ® Plus,
  ## GF=Globalfiler ® , PP16=PowerPlex ® 16 HS)
  kit="(?:IP|GF|PP16)",
  "-",
  ## If the extracts were quantified with both a small and large
  ## autosomal fragment, the Quality Index is presented in the form of
  ## a Q value (Q0.8 in this example). In some instances, the Q
  ## designator is followed by “LAND,” which stands for “large
  ## autosomal not detected.” This term is used for samples in which
  ## the large autosomal fragment was not detected during qPCR and so
  ## a numerical Q value was not obtained.
  nc::quantifier(nc::field("Q", "", "[^_]+"), "_", "?"),
  function(x)ifelse(x=="LAND", -Inf, as.numeric(x)),
  ## The designator 002 is the capillary number, where capillaries are
  ## numbered 001-004 on the 3130 electrophoresis platform and 01-08
  ## on the 3500 platform.
  capillary="[0-9]+",
  "[.]",
  injection.seconds="[0-9]+", as.integer)
file.dt[, identifiers.list := strsplit(identifiers.str, "_")]
file.dt[, ratios.list := lapply(strsplit(ratios.str, ";"), as.numeric)]
