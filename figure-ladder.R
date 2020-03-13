source("packages.R")
file.vec <- c(
  mixture="data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(040416LEA_5sec)/A01_Ladder-IP_001.5sec.fsa",
  one.part="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(022916LEA_20sec)/A01_Ladder-IP_001.20sec.fsa",
  four.parts="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(030216LEA_20sec)/E09_Ladder-IP_001.20sec.fsa")
ladder.vec <- Sys.glob(file.path(dirname(file.vec), "*Ladder*"))

fsa2dt <- function(mix.fsa){
  L <- seqinr::read.abif(mix.fsa)
  D <- L$Data
  l.vec <- sapply(D, length)
  channel.names <- names(l.vec)[l.vec==max(l.vec)]
  dye.vec <- do.call(c, D[paste0("DyeN.", 1:5)])
  dye.colors <- c(
    "6-FAM"="blue",
    "VIC"="green",
    "NED"="yellow",
    "PET"="red",
    "LIZ"="orange")
  channel.dt.list <- list()
  for(channel.i in seq_along(channel.names)){
    channel <- channel.names[[channel.i]]
    dye <- dye.vec[[channel.i]]
    RFU <- D[[channel]]
    channel.dt.list[[channel]] <- data.table(
      dye, color=dye.colors[[dye]],
      channel=as.integer(sub(".*[.]", "", channel)),
      Time=seq_along(RFU), RFU)
  }
  do.call(rbind, channel.dt.list)
}

ladder.channels.list <- list()
for(file.i in seq_along(file.vec)){
  fsa <- file.vec[[file.i]]
  channels <- fsa2dt(fsa)
  ladder.channels.list[[file.i]] <- data.table(
    sample=names(file.vec)[[file.i]], channels)
}
ladder.channels <- do.call(rbind, ladder.channels.list)

## TODO confront Identifier_Plus_* data with the peaks in these ladder
## files.

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap(
    channel ~ sample,
    scales="free",
    labeller=label_both, ncol=3)+
  geom_line(aes(
    Time, RFU),
    data=ladder.channels)
png("figure-ladder.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()


