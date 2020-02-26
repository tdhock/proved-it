source("packages.R")
fsa_db <- fread("fsa_db.csv")
fsa_present <- fsa_db[file.exists(fsa)]
fsa_present[, type := data.table::fcase(
  control != "", "control",
  sample.identifier != "", "single",
  identifiers.str != "", "mixture",
  default="unknown")]
fsa_present[, file.i := 1:.N]
fsa_present[, table(type, kit)]
mix.dt <- fsa_present[type=="mixture" & grepl("a$", mix.dilution)]
mix.dt[, identifiers.list := lapply(strsplit(identifiers.str, "_"), as.integer)]
mix.dt[, n.sources := sapply(identifiers.list, length)]
mix.dt[, mix.i := 1:.N ]
mix.samples <- mix.dt[, data.table(
  project,
  sample.identifier=identifiers.list[[1]]
), by=.(mix.i)]
##one mixture.
mix.dt[, table(kit, n.sources)]
two.sources <- unique(mix.dt[n.sources==2, .(project, identifiers.str)])
one.pair <- two.sources[1]
one.pair.mix <- mix.dt[one.pair, on=.(project, identifiers.str)]
one.pair.mix[, table(template.nanograms, injection.seconds)]

fsa2dt <- function(mix.fsa){
  mix.fsa.data <- seqinr::read.abif(mix.fsa)
  l.vec <- sapply(mix.fsa.data$Data, length)
  channel.names <- names(l.vec)[l.vec==max(l.vec)]
  channel.dt.list <- list()
  for(channel in channel.names){
    RFU <- mix.fsa.data$Data[[channel]]
    channel.dt.list[[channel]] <- data.table(
      channel=as.integer(sub(".*[.]", "", channel)),
      Time=seq_along(RFU), RFU)
  }
  do.call(rbind, channel.dt.list)
}

one.pair.mix.channels.list <- list()
for(file.i in 1:nrow(one.pair.mix)){
  one.pair.row <- one.pair.mix[file.i]
  fsa <- fsa_present[file.i, fsa]
  channels <- fsa2dt(fsa)
  one.pair.mix.channels.list[[file.i]] <- data.table(
    one.pair.row, channels)
}
one.pair.mix.channels <- do.call(rbind, one.pair.mix.channels.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU,
    group=paste(template.nanograms, injection.seconds),
    color=template.nanograms,
    linetype=paste(injection.seconds)),
    data=one.pair.mix.channels[Time>2500])
print(gg)

png("figure-alignment-nanograms-%d.png", w=15, h=8, units="in", res=100)
for(ng in unique(one.pair.mix.channels$template.nanograms)){
  some.nanograms <- one.pair.mix.channels[template.nanograms==ng]
  gg <- ggplot()+
    ggtitle(paste0(
      some.nanograms$lane[1],
      ", two component mixture, ",
      ng, "ng"))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(channel ~ injection.seconds, labeller=label_both)+
    geom_line(aes(
      Time, RFU),
      data=some.nanograms)+
    coord_cartesian(ylim=c(0, 1000))
  print(gg)
}
dev.off()

png("figure-alignment-nanograms-wrap-%d.png", w=15, h=8, units="in", res=100)
for(ng in unique(one.pair.mix.channels$template.nanograms)){
  some.nanograms <- one.pair.mix.channels[template.nanograms==ng]
  gg <- ggplot()+
    ggtitle(paste0(
      some.nanograms$lane[1],
      ", two component mixture, ",
      ng, "ng"))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    ##facet_grid(channel ~ injection.seconds, labeller=label_both)+
    facet_wrap(
      channel ~ injection.seconds,
      labeller=label_both, ncol=3, scales="free")+
    geom_line(aes(
      Time, RFU),
      data=some.nanograms[Time>2500])
  print(gg)
}
dev.off()

is <- 5
ng <- 0.075
one.f <- one.pair.mix[
  injection.seconds == is & template.nanograms==ng]
one.mix.samples <- one.f[, .(
  project, sample.identifier=identifiers.list[[1]])]
one.ok <- one.pair.mix.channels[
  injection.seconds == is & template.nanograms==ng]
gg <- ggplot()+
  ggtitle(one.f$fsa)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both, scales="free")+
  geom_line(aes(
    Time, RFU),
    data=one.ok[2500<Time & Time<6000])
png("figure-alignment-nanograms-mixture.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()

single.dt <- fsa_present[type=="single" & grepl("a$", single.dilution)]
(one.singles <- single.dt[
  one.mix.samples, on=.(project, sample.identifier)
][order(sample.identifier, template.nanograms, injection.seconds), .(
  fsa, sample.identifier, template.nanograms, injection.seconds)])
one.singles[, table(template.nanograms, injection.seconds, sample.identifier)]

show.singles <- one.singles[template.nanograms==0.0625]
one.pair.single.channels.list <- list()
for(file.i in 1:nrow(show.singles)){
  one.single.row <- show.singles[file.i]
  fsa <- fsa_present[file.i, fsa]
  channels <- fsa2dt(fsa)
  one.pair.single.channels.list[[file.i]] <- data.table(
    one.single.row, channels)
}
one.pair.single.channels <- do.call(rbind, one.pair.single.channels.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel ~ injection.seconds + sample.identifier,
    labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=one.pair.single.channels[2500<Time])
png("figure-alignment-nanograms-single.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()

one.pair.single.good <- one.pair.single.channels[injection.seconds==20]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel ~ sample.identifier,
    labeller=label_both, scales="free")+
  geom_line(aes(
    Time, RFU),
    data=one.pair.single.good[2500<Time])
png("figure-alignment-nanograms-single-good.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()


png("figure-alignment-nanograms-%d.png", w=15, h=8, units="in", res=100)
for(ng in unique(one.pair.mix.channels$template.nanograms)){
  some.nanograms <- one.pair.mix.channels[template.nanograms==ng]
  gg <- ggplot()+
    ggtitle(paste0(
      some.nanograms$lane[1],
      ", two component mixture, ",
      ng, "ng"))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(channel ~ injection.seconds, labeller=label_both)+
    geom_line(aes(
      Time, RFU),
      data=some.nanograms)+
    coord_cartesian(ylim=c(0, 1000))
  print(gg)
}
dev.off()

png("figure-alignment-nanograms-wrap-%d.png", w=15, h=8, units="in", res=100)
for(ng in unique(one.pair.mix.channels$template.nanograms)){
  some.nanograms <- one.pair.mix.channels[template.nanograms==ng]
  gg <- ggplot()+
    ggtitle(paste0(
      some.nanograms$lane[1],
      ", two component mixture, ",
      ng, "ng"))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    ##facet_grid(channel ~ injection.seconds, labeller=label_both)+
    facet_wrap(
      channel ~ injection.seconds,
      labeller=label_both, ncol=3, scales="free")+
    geom_line(aes(
      Time, RFU),
      data=some.nanograms[Time>2500])
  print(gg)
}
dev.off()
