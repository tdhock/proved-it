source("packages.R")
fsa_db <- fread("fsa_db.csv")
fsa_db[, type := data.table::fcase(
  control != "", "control",
  sample.identifier != "", "single",
  identifiers.str != "", "mixture",
  default="unknown")]
fsa_db[, file.i := 1:.N]
mix.dt <- fsa_db[type=="mixture"]
mix.dt[, identifiers.list := lapply(strsplit(identifiers.str, "_"), as.integer)]
## How many sources are in each mixture?
(mix.tab <- table(sapply(mix.dt$identifiers.list, length)))
mix.dt[, mix.i := 1:.N ]
mix.samples <- mix.dt[, data.table(
  project,
  sample.identifier=identifiers.list[[1]]
), by=.(mix.i, file.i)]
single.dt <- fsa_db[type=="single"]

##one mixture.
one.mix <- mix.dt[1]

## Are there other samples with the same mixture components?
other.mix <- mix.dt[one.mix$identifiers.str, on="identifiers.str"]

## Do all of the other samples have the same ratios?
table(other.mix$ratios.str)

other.mix[, Q.num := as.numeric(ifelse(Q=="LAND", 0, Q))]
most.template.injection <- other.mix[
  template.nanograms==max(template.nanograms) &
    injection.seconds == max(injection.seconds)]
other.mix[, selected := file.i %in% most.template.injection$file.i]
ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(injection.seconds ~ ., labeller=label_both)+
  geom_point(aes(
    template.nanograms, Q.num, color=selected),
    data=other.mix)+
  scale_y_log10()

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

mix.raw.dt.list <- list()
for(other.i in 1:nrow(most.template.injection)){
  other.row <- most.template.injection[other.i]
  other.row.dt <- fsa2dt(other.row$fsa)
  mix.raw.dt.list[[other.i]] <- data.table(
    other.row[, .(Q.num, file.i)],
    other.row.dt)
}
mix.raw.dt <- do.call(rbind, mix.raw.dt.list)

ladder.dir <- dirname(other.row$fsa)
dir(ladder.dir)
ladder.fsa.vec <- Sys.glob(file.path(ladder.dir, "*Ladder*"))
ladder.dt.list <- list()
for(ladder.fsa in ladder.fsa.vec){
  ladder.row <- fsa_db[ladder.fsa, on="fsa"]
  one.ladder <- fsa2dt(ladder.fsa)
  ladder.dt.list[[ladder.fsa]] <- data.table(
    ladder.row[, .(letter, numbers, capillary)],
    one.ladder)
}
ladder.dt <- do.call(rbind, ladder.dt.list)

gg <- ggplot()+
  ggtitle(paste0(
    "All ladder samples in ", ladder.dir))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ letter + numbers + capillary, labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=ladder.dt)
png("figure-one-mixture-ladder-zoom-out.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()

gg+coord_cartesian(
  xlim=c(3000, 4000))

gg <- ggplot()+
  ggtitle(paste0(
    "All ladder samples in ", ladder.dir))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ letter + numbers + capillary, labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=ladder.dt[2500<Time & Time < 4000])
png("figure-one-mixture-ladder-zoom-in.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()

gg <- ggplot()+
  ggtitle(paste0(
    "One mixture, max injection seconds, ",
    "max template nanograms, different Q.num values\n",
    other.row$fsa))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ Q.num + file.i, labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=mix.raw.dt)
png("figure-one-mixture.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()
