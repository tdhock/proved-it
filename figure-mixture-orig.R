source("packages.R")
fsa.vec <- Sys.glob("data/*/*/*/*.fsa")
match.dt <- nc::capture_first_vec(
  fsa.vec,
  "data/PROVEDIt_",
  persons=".*?",
  "-Person Profiles")
match.dt[, table(persons)]

prefix.pattern <- list(
  "/",
  letter="[A-Z]",
  numbers="[0-9]+",
  "[-_]")
mixture.pattern <- list(
  identifiers.str="[^-]+", # within that project.
  "-",
  ratios.str="[^-]+",
  nc::quantifier("-", mix.dilution="M[^-]+", "?"))
project.pattern <- list(
  "RD[0-9]+",
  "-",
  "[0-9]+")
single.pattern <- list(
  ## It is to be noted that each sample is designated by the
  ## combination of project number and sample identifier. For example,
  ## RD14-0003-21 will have a different known genotype than
  ## RD12-0002-21.
  sample.identifier="[0-9]+", as.integer, # within that project.
  ## d_ is the dilution number which was used by laboratory personnel
  ## to distinguish between extracts.
  nc::quantifier(
    "d",
    single.dilution=".*?",
    ## the ‘x’ designator indicates DNA condition (Supplementary Methods
    ## Table 1)
    DNA.condition="x?",
    "?"))
nanograms.pattern <- nc::quantifier(
  template.nanograms="[0-9.]+", as.numeric, "?")
suffix.pattern <- list(
  ## IP is the amplification kit type (i.e., IP=Identifiler ® Plus,
  ## GF=Globalfiler ® , PP16=PowerPlex ® 16 HS)
  nc::quantifier("-", nanograms.pattern, kit="[^-_]+", "?"),
  ## If the extracts were quantified with both a small and large
  ## autosomal fragment, the Quality Index is presented in the form of
  ## a Q value (Q0.8 in this example). In some instances, the Q
  ## designator is followed by “LAND,” which stands for “large
  ## autosomal not detected.” This term is used for samples in which
  ## the large autosomal fragment was not detected during qPCR and so
  ## a numerical Q value was not obtained.
  nc::quantifier("-", nc::field("Q", "", "[^_]+"), "?"),
  "[-_]",
  ## The designator 002 is the capillary number, where capillaries are
  ## numbered 001-004 on the 3130 electrophoresis platform and 01-08
  ## on the 3500 platform.
  capillary="[0-9]+",
  "[.]",
  injection.seconds="[0-9]+", as.integer,
  "sec")
fsa.pattern <- list(
  prefix.pattern,
  lane=list(
    control="[^-_]+",
    "|",
    project=project.pattern,
    "-",
    extra=list(
      single.pattern,
      "|",
      mixture.pattern)),
  suffix.pattern)
subject.vec <- c(
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(080114LEA_10sec)/A01-Ladder-PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(120715CMG_10sec)/A12_RB121514ADG_001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(011816CMG_10sec)/A10_RB102191515LEA-IP_001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(080114LEA_10sec)/A02-RD12-0002-35-0.5PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(032516LEA_10sec)/G01_RD14-0003-35d3S30-0.01563P-Q10.0_003.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(011816CMG_10sec)/A06_RD14-0003-24d3a-0.0625IP-Q0.8_001.10sec.fsa",
  ## 1 above 2-5 below.
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051616LEA_5sec)/H07_RD14-0003-35_36_37_38_39-1;4;4;4;1-M2I35-0.75IP-QLAND_004.5sec.fsa",
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051616LEA_5sec)/D07_Ladder-IP_004.5sec.fsa")
nc::capture_first_vec(subject.vec, fsa.pattern)

match.dt <- nc::capture_first_vec(fsa.vec, fsa.pattern)

match.dt[, type := data.table::fcase(
  control != "", "control",
  sample.identifier != "", "single",
  identifiers.str != "", "mixture",
  default="unknown")]
table(match.dt$type)

## What is the P kit type? typo?
table(match.dt$kit)
match.dt[, which(kit=="P")]
i <- 9264
fsa.vec[i]
match.dt[i]

## There is supposed to be a GF kit type but I don't see any.
grep("GF", fsa.vec, value=TRUE)

## Is there a single-source fsa file for each sample referenced in the
## mixture file names?
match.dt[, file.i := 1:.N]
mix.dt <- match.dt[type=="mixture"]
mix.dt[, identifiers.list := lapply(strsplit(identifiers.str, "_"), as.integer)]

## How many sources are in each mixture?
(mix.tab <- table(sapply(mix.dt$identifiers.list, length)))

mix.dt[, mix.i := 1:.N ]
mix.samples <- mix.dt[, data.table(
  project,
  sample.identifier=identifiers.list[[1]]
), by=.(mix.i, file.i)]

single.dt <- match.dt[type=="single"]

count.dt <- single.dt[mix.samples, .(
  count=.N,
  mix.i
  ##mix.file=fsa.vec[mix.file.i]
), on=.(project, sample.identifier), by=.EACHI][order(count)]

## all samples are present.
count.dt[count==0]

mix.fsa <- fsa.vec[mix.dt[1, file.i]]
mix.fsa.data <- seqinr::read.abif(mix.fsa)

n.channels <- 5
png("figure-plotabif-%d.png", w=10, h=3, units="in", res=100)
for(channel.i in 1:n.channels){
  seqinr::plotabif(mix.fsa.data, chanel=channel.i)
}
dev.off()
