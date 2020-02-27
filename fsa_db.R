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
  well.letter="[A-Z]",
  well.number="[0-9]+", as.integer,
  "[-_]")
project.pattern <- list(
  "RD[0-9]+",
  "-",
  "[0-9]+")
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
  capillary="[0-9]+", as.integer,
  "[.]",
  injection.seconds="[0-9]+", as.integer,
  "sec")
single.pattern <- list(
  ## It is to be noted that each sample is designated by the
  ## combination of project number and sample identifier. For example,
  ## RD14-0003-21 will have a different known genotype than
  ## RD12-0002-21.
  sample.identifier="[0-9]+", as.integer)  #within that project.
mixture.pattern <- list(
  identifiers.str="[0-9_]+",
  "-",
  ratios.str="[0-9;]+",
  "-?")
alternatives <- function(...){
  in.list <- list(...)
  stopifnot(1 < length(in.list))
  out.list <- in.list[1]
  for(i in 2:length(in.list)){
    out.list[[length(out.list)+1L]] <- "|"
    out.list[[length(out.list)+1L]] <- in.list[[i]]
  }
  out.list
}
single.or.mixture <- alternatives(
  single.pattern,
  mixture.pattern)
treatment.symbol.names <- c(
  S="sonication",
  I="humic.acid",
  U="ultraviolet",
  "-"="fragmentase",
  a="none")
DNase.units <- c(b=3, c=6, d=12, e=24)
treatment.symbol.names[names(DNase.units)] <- "DNase"
treatment.name.units <- c(
  DNase="mU",
  fragmentase="minutes",
  ultraviolet="minutes",
  sonication="cycles",
  humic.acid="μL")
treatment.pattern <- list(
  treatment.symbol="[-a-eUSI]",
  treatment.amount=nc::quantifier("[0-9]+", "?"), as.integer)
maybe.number.treatment <- nc::quantifier(
  ## d_ is the dilution number which was used by laboratory personnel
  ## to distinguish between extracts.
  dilution.number="[0-9]+", as.integer,
  treatment=nc::quantifier(treatment.pattern, "?"),
  "?")
maybe.dilution <- nc::quantifier(
  dilution.letter="[dM]",
  maybe.number.treatment,
  "?")
extra.pattern <- list(
  single.or.mixture,
  maybe.dilution)
fsa.pattern <- list(
  prefix.pattern,
  lane=list(
    control="[^-_]+",
    "|",
    project=project.pattern,
    "-",
    extra=extra.pattern),
  suffix.pattern)
capture_convert <- function(fsa){
  m.dt <- nc::capture_first_vec(fsa, fsa.pattern)
  m.dt[, treatment.name := treatment.symbol.names[treatment.symbol] ]
  m.dt[, treatment.units := treatment.name.units[treatment.name] ]
  m.dt[
    treatment.symbol %in% names(DNase.units),
    treatment.amount := DNase.units[treatment.symbol] ]
  m.dt[, type := data.table::fcase(
    control != "", "control",
    sample.identifier != "", "single",
    identifiers.str != "", "mixture",
    default="unknown")]
  m.dt[, Q.num := as.numeric(ifelse(Q=="LAND", 0, Q))]
  data.table(fsa, m.dt[, .(
    well.letter, well.number, control, project, sample.identifier,
    identifiers.str, ratios.str, dilution.letter, dilution.number,
    treatment.name, treatment.amount, treatment.units,
    template.nanograms, kit, Q.num, capillary, injection.seconds, type)])
}
subject.vec <- c(
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(080114LEA_10sec)/A01-Ladder-PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(120715CMG_10sec)/A12_RB121514ADG_001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(011816CMG_10sec)/A10_RB102191515LEA-IP_001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(080114LEA_10sec)/A02-RD12-0002-35-0.5PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(032516LEA_10sec)/G01_RD14-0003-35d3S30-0.01563P-Q10.0_003.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_IDPlus28cycles/10 sec/RD14-0003(011816CMG_10sec)/A06_RD14-0003-24d3a-0.0625IP-Q0.8_001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(090514LEA_10sec)/A08-RD12-0002-01d-0.125PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(090514LEA_10sec)/A10-RD12-0002-04d1-0.0625PP16-001.10sec.fsa",
  "data/PROVEDIt_1-Person Profiles_3130 5sec_IDPlus28cycles/5 sec/RD14-0003(121715CMG_5sec)/C02_RD14-0003-15d2b-0.25IP-Q0.5_003.5sec.fsa",
  ## 1 above 2-5 below.
  "data/PROVEDIt_2-5-Person Profiles_3130 10sec_PP16HS32cycles/10sec/RD12-0002(102414LEA_10sec)/A02-RD12-0002-1_2-1;9-0.125PP16-001.10sec.fsa",
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051616LEA_5sec)/H07_RD14-0003-35_36_37_38_39-1;4;4;4;1-M2I35-0.75IP-QLAND_004.5sec.fsa",
  "data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(051616LEA_5sec)/D07_Ladder-IP_004.5sec.fsa")
(test.dt <- capture_convert(subject.vec))

(match.dt <- capture_convert(fsa.vec))

fsa_db <- data.table(fsa=fsa.vec, match.dt)
fwrite(fsa_db, "fsa_db.csv")
