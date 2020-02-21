source("packages.R")
fsa_db <- fread("fsa_db.csv")
fsa_db[, type := data.table::fcase(
  control != "", "control",
  sample.identifier != "", "single",
  identifiers.str != "", "mixture",
  default="unknown")]
fsa_db[, file.i := 1:.N]

## Untreated RD14-0003-19d2a-0.5IP-Q0.8_002.5sec a Letter a indicates
## that the sample was not treated with any protocol intended to
## induce PCR inefficiencies.

## DNase I Degradation RD14-0003-21d2b-0.5IP-Q0.8_002.20sec b to e
## Letters b to e indicate units of DNase I in degradation
## reaction. (b) 3, (c) 6, (d) 12, and (e) 24 mU enzyme.
c(b=3, c=6, d=12, e=24)

## Fragmentase ® Degradation RD14-0003-36d1-15-0.5IP-Q1.4_003.10sec
## -15 to -45, -15 indicates enzyme digestion/incubation time in
## minutes. 15, 30, and 45 minute digestion times were utilized.

## UV Damage RD14-0003-03d2U60-0.5IP-Q8.8_003.10sec, U15 to U105, U60
## indicates 60 minutes of UV exposure. Times ranged from 15-105
## minutes.

## Sonication RD14-0003-12d3S30-0.0078IP-Q17.1_002.20sec, S2 to S30,
## S30 indicates DNA was damaged with 30 cycles of sonication. 2, 10,
## and 30 cycles were utilized.

## Humic Acid Inhibition RD14-0003-49d2I22-0.5IP-Q2.4_002.10sec, I15
## to I35, I22 indicates volume of 2 mg/mL humic acid (in μL) added to
## whole blood lysate. Three volumes were utilized (15, 22, and 35
## μL).
