source("packages.R")
file.vec <- c(
  mixture="data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(040416LEA_5sec)/A03_RD14-0003-40_41-1;4-M4a-0.075IP-Q0.6_001.5sec.fsa",
  one.part="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(022916LEA_20sec)/A06_RD14-0003-40d3a-0.0625IP-Q1.1_001.20sec.fsa",
  four.parts="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(030216LEA_20sec)/E01_RD14-0003-41d3a-0.0625IP-Q1.0_001.20sec.fsa")

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

one.pair.channels.list <- list()
for(file.i in seq_along(file.vec)){
  fsa <- file.vec[[file.i]]
  channels <- fsa2dt(fsa)
  one.pair.channels.list[[file.i]] <- data.table(
    sample=names(file.vec)[[file.i]], channels)
}
one.pair.channels <- do.call(rbind, one.pair.channels.list)

ip.sizes <- c(
  ##20, 40,
  ##60,
  80, 100,
  114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314,
  320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540,
  560, 580, 600)
(n.peaks <- length(ip.sizes))

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap(
    channel ~ sample,
    scales="free",
    labeller=label_both, ncol=3)+
  geom_line(aes(
    Time, RFU),
    data=one.pair.channels[Time>2500])
print(gg)

pdir <- "figure-one-pair-data"
if(FALSE){
  unlink(pdir, recursive=TRUE)
}
standard.after <- one.pair.channels[Time>2500 & channel==105]
max.dt <- standard.after[, {
  non.neg <- RFU-min(RFU)
  rle.list <- rle(non.neg)
  edge.dt <- data.table(
    Time=cumsum(c(Time[1]-1L, rle.list$lengths)))
  coverage.dt <- edge.dt[, data.table(
    chrom = "chrUnknown",
    chromStart = Time[-.N],
    chromEnd = Time[-1],
    count = rle.list$values)]
  sdir <- file.path(pdir, sample)
  dir.create(sdir, showWarnings=FALSE, recursive=TRUE)
  coverage.bedGraph <- file.path(sdir, "coverage.bedGraph")
  PeakSegDisk::writeBedGraph(coverage.dt, coverage.bedGraph)
  fit <- PeakSegDisk::sequentialSearch_dir(sdir, n.peaks)
  peak.dt <- fit$segments[status=="peak"][order(chromStart), .(
    mean,
    peakStart=chromStart+0.5,
    peakEnd=chromEnd+0.5,
    bases=ip.sizes)]
  setkey(peak.dt, peakStart, peakEnd)
  rfu.dt <- data.table(Time, RFU, Time0=Time)
  setkey(rfu.dt, Time, Time0)
  rfu.in.peaks <- foverlaps(rfu.dt, peak.dt, nomatch=0L)
  rfu.in.peaks[, .SD[which.max(RFU)], by=bases]
}, by=sample]

coef.dt <- max.dt[, {
  data.table(t(coef(lm(bases ~ Time))))
}, by=sample]

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_abline(aes(
    slope=Time, intercept=`(Intercept)`),
    data=coef.dt)+
  geom_point(aes(
    Time, bases),
    data=max.dt)+
  facet_grid(. ~ sample)
print(gg)

one.pair.channels[, bases := {
  samp.max <- max.dt[sample, on="sample"]
  approx(samp.max$Time, samp.max$bases, Time)$y
}, by=sample]

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap(
    channel ~ sample,
    scales="free",
    labeller=label_both, ncol=3)+
  geom_line(aes(
    bases, RFU),
    data=one.pair.channels[!is.na(bases)])
print(gg)

dput(RColorBrewer::brewer.pal(Inf, "Paired"))
sample.colors <- c(
  mixture="black",
  ## one.part="#A6CEE3", #blues
  ## four.parts="#1F78B4",
  ## one.part="#B2DF8A", #green
  ## four.parts="#33A02C",
  one.part="#FB9A99", #red
  four.parts="#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
some.data <- one.pair.channels[!is.na(bases) & bases<400]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., scales="free", labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  geom_line(aes(
    bases, RFU, color=sample),
    data=some.data)
print(gg)

some.data[, Sample := factor(sample, names(file.vec))]
some.data[, Channel := factor(channel, 1:105)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel + Sample ~ ., labeller=label_both, scales="free")+
  ##scale_color_manual(values=sample.colors)+
  geom_line(aes(
    bases, RFU, color=Channel),
    size=1, alpha=0.5,
    data=some.data)
png("figure-one-pair-panels.png", w=15, h=18, units="in", res=100)
print(gg)
dev.off()

some.data[, relative.height := (RFU-mean(RFU[bases>350]))/max(RFU), by=.(channel, sample)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  geom_line(aes(
    bases, relative.height, color=sample),
    size=1, alpha=0.5,
    data=some.data)+
  scale_y_continuous(breaks=seq(0,1, by=0.2))
png("figure-one-pair.png", w=15, h=8, units="in", res=100)
print(gg)
dev.off()
