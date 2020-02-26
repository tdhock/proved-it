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
two.sources <- mix.dt[n.sources==2]
one.mix <- two.sources[1]
one.mix.samples <- mix.samples[one.mix$mix.i == mix.i]
single.dt <- fsa_present[type=="single" & grepl("a$", single.dilution)]
(one.singles <- single.dt[
  one.mix.samples, on=.(project, sample.identifier)
][
   ##template.nanograms==one.mix$template.nanograms
   template.nanograms==max(template.nanograms)
][
  injection.seconds==one.mix$injection.seconds
 ])

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
both.cols <- c("lane", "file.i")
both.dt <- rbind(
  data.table(sample="mixture", one.mix[, both.cols, with=FALSE]),
  data.table(sample="single", one.singles[, both.cols, with=FALSE]))
both.channels.dt.list <- list()
for(both.i in 1:nrow(both.dt)){
  both.row <- both.dt[both.i]
  fsa <- fsa_present[both.i, fsa]
  channels <- fsa2dt(fsa)
  both.channels.dt.list[[both.i]] <- data.table(
    both.row, channels)
}
both.channels.dt <- do.call(rbind, both.channels.dt.list)
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(lane + channel ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=both.channels.dt)
print(gg)

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(lane + channel ~ ., labeller=label_both, scales="free")+
  geom_line(aes(
    Time, RFU),
    data=both.channels.dt[Time>2500])
print(gg)


both.105 <- both.channels.dt[channel==105]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(lane ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=both.105)
print(gg)

