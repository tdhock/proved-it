subject.vec <- c(
  "RD14-0003-21d1x-0.5IP-Q0.8_002.20sec",
  "RD12-0002-21d1-0.5IP-002.20sec")
## Comments explaining the regex copied from
## https://lftdi.camden.rutgers.edu/wp-content/uploads/2019/12/PROVEDIt-Database-Naming-Convention-Laboratory-Methodsv1.pdf
match.dt <- nc::capture_first_vec(
  subject.vec,
  project.number=list(
    ".*?",
    "-",
    "[0-9]+"),
  "-",
  ## It is to be noted that each sample is designated by the
  ## combination of project number and sample identifier. For example,
  ## RD14-0003-21 will have a different known genotype than
  ## RD12-0002-21.
  sample.identifier="[0-9]+", # within that project.
  ## d_ is the dilution number which was used by laboratory personnel
  ## to distinguish between extracts.
  "d",
  dilution="[0-9]+", as.integer,
  ## the ‘x’ designator indicates DNA condition (Supplementary Methods
  ## Table 1)
  DNA.condition="x?",
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
str(match.dt)
