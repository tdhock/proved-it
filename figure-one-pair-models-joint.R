source("packages.R")
file.vec <- c(
  mixture="data/PROVEDIt_2-5-Person Profiles_3130 5sec_IDPlus28cycles/5sec/RD14-0003(040416LEA_5sec)/A03_RD14-0003-40_41-1;4-M4a-0.075IP-Q0.6_001.5sec.fsa",
  one.part="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(022916LEA_20sec)/A06_RD14-0003-40d3a-0.0625IP-Q1.1_001.20sec.fsa",
  four.parts="data/PROVEDIt_1-Person Profiles_3130 20sec_IDPlus28cycles/20 sec/RD14-0003(030216LEA_20sec)/E01_RD14-0003-41d3a-0.0625IP-Q1.0_001.20sec.fsa")

ladder.vec <- Sys.glob(file.path(dirname(file.vec), "*Ladder*"))

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

one.pair.mid.time <- one.pair.channels[2500<Time & Time<9000]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap(
    channel ~ sample,
    scales="free",
    labeller=label_both, ncol=3)+
  geom_line(aes(
    Time, RFU),
    data=one.pair.mid.time)
print(gg)

pdir <- "figure-one-pair-models-data"
if(FALSE){
  unlink(pdir, recursive=TRUE)
}
standard.after <- one.pair.mid.time[channel==105]
fit.dt <- standard.after[, {
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
  data.table(
    fit$loss,
    segs=list(fit$segments),
    others=list(fit$others),
    rfu=list(data.table(Time, RFU, Time0=Time)))
}, by=sample]

max.dt <- fit.dt[, {
  peak.dt <- segs[[1]][status=="peak"][order(chromStart), .(
    mean,
    peakStart=chromStart+0.5,
    peakEnd=chromEnd+0.5,
    bases=ip.sizes)]
  setkey(peak.dt, peakStart, peakEnd)
  rfu.dt <- rfu[[1]]
  setkey(rfu.dt, Time, Time0)
  rfu.in.peaks <- foverlaps(rfu.dt, peak.dt, nomatch=0L)
  rfu.in.peaks[, .SD[which.max(RFU)], by=bases]
}, by=sample]

one.pair.channels[, bases := {
  samp.max <- max.dt[sample, on="sample"]
  approx(samp.max$Time, samp.max$bases, Time)$y
}, by=sample]

min.bases <- 90
max.bases <- 370
bases.dt <- data.table(base=min.bases:max.bases)
one.pair.middle <- one.pair.channels[min.bases<bases & bases<max.bases]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_wrap(
    channel ~ sample,
    scales="free",
    labeller=label_both, ncol=3)+
  geom_line(aes(
    bases, RFU),
    data=one.pair.middle)
print(gg)

one.pair.middle[, relative.height := {
  m <- mean(RFU[bases>350])
  norm <- function(x)x-m
  norm(RFU)/norm(max(RFU))
}, by=.(channel, sample)]
sample.colors <- c(
  mixture="black",
  ## one.part="#A6CEE3", #blues
  ## four.parts="#1F78B4",
  ## one.part="#B2DF8A", #green
  ## four.parts="#33A02C",
  one.part="#FB9A99", #red
  four.parts="#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  geom_line(aes(
    bases, relative.height, color=sample),
    size=1, alpha=0.5,
    data=one.pair.middle)+
  scale_y_continuous(breaks=seq(0,1, by=0.2))
print(gg)

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_vline(aes(
    xintercept=base),
    color="grey",
    data=bases.dt)+
  facet_grid(channel ~ ., labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  geom_point(aes(
    bases, RFU, color=sample),
    size=1, alpha=0.5,
    data=one.pair.middle)
print(gg)

one.pair.middle[, milli.bases := as.integer(bases*1000)]
one.pair.middle[, diff.prev := c(NA, diff(milli.bases)), by=.(sample, channel)]
one.pair.middle[, offset := as.integer(
  c(diff.prev[2], diff.prev[-1])/2), by=.(sample, channel)]
one.pair.middle[, chromStart := milli.bases-offset]
one.pair.middle[, chromEnd := c(
  chromStart[-1], milli.bases[.N]+offset[.N]), by=.(sample, channel)]
one.pair.middle[, count := RFU-min(RFU), by=.(sample, channel)]
one.pair.middle[, Sample := factor(
  sample, c("mixture", "four.parts", "one.part"))]
ggplot()+
  geom_point(aes(
    bases, (chromEnd+chromStart)/2000),
    data=one.pair.middle)+
  coord_equal()

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both, scales="free")+
  scale_color_manual(values=sample.colors)+
  geom_step(aes(
    chromStart, count, color=sample),
    size=1, alpha=0.5,
    data=one.pair.middle)+
  coord_cartesian(ylim=c(0,100))
print(gg)

p <- function(channel, problemStart, problemEnd){
  data.table(channel, problemStart=problemStart-0.5, problemEnd=problemEnd+0.5)
}
problems <- rbind(
  p(4, 102, 108),
  p(4, 108, 114),
  p(4, 143, 147),
  p(4, 147, 151),
  p(4, 151, 155),
  p(4, 155, 159),
  p(4, 231, 235),
  p(4, 248, 252),
  p(4, 252, 256),
  p(4, 256, 260),
  p(3, 105, 111),
  p(3, 114, 118),
  p(3, 124, 128),
  p(3, 128, 132),
  ##p(3, 160, 166),
  p(3, 166, 170),
  p(3, 174, 178),
  p(3, 178, 182),
  p(3, 182, 186),
  p(3, 226, 232),
  p(3, 242, 248),
  p(3, 248, 260),
  p(3, 275, 279),
  p(3, 279, 283),
  p(3, 283, 287),
  p(3, 287, 300),
  p(3, 306, 312)
)

one.channel <- one.pair.middle[channel==2 & -Inf<bases & bases<Inf]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  scale_color_manual(values=sample.colors)+
  geom_line(aes(
    bases, count, color=Sample),
    size=1, alpha=0.5,
    data=one.channel)
print(gg)

problems[, problemMid := (problemEnd+problemStart)/2]
problems[, problem.i := 1:.N, by=channel]
setkey(problems, channel, problemStart, problemEnd)
one.pair.middle[, bases0 := bases]
setkey(one.pair.middle, channel, bases, bases0)
prob.data <- foverlaps(problems, one.pair.middle, nomatch=0L)
prob.data[, bases.from.middle := bases-problemMid]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel ~ problem.i, scales="free", space="free", labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(breaks=seq(-1000, 1000, by=2))+
  geom_line(aes(
    bases.from.middle, relative.height, color=Sample),
    size=1, alpha=0.5,
    data=prob.data)
print(gg)

seg.dt <- prob.data[, {
  one.prob <- data.table(
    sample.id=sample,
    sample.group=sample,
    chromStart, chromEnd, count)
  plist <- PeakSegJoint::ProfileList(one.prob)
  fit <- PeakSegJoint::PeakSegJointFaster(plist)
  segStart <- with(fit, c(data_start_end[1], peak_start_end))
  segEnd <- with(fit, c(peak_start_end, data_start_end[2]))
  with(fit, data.table(
    segStart,
    segEnd,
    mean=as.numeric(t(mean_mat)),
    status=c("background", "peak", "background"),
    sample=sub(".*/", "", rep(rownames(mean_mat), each=3))))
}, by=.(channel, problem.i)][problems, on=.(channel, problem.i)]

for(millibases.name in c("segStart", "segEnd")){
  millibases <- seg.dt[[millibases.name]]
  norm.name <- paste0("norm.", millibases.name)
  seg.dt[[norm.name]] <- millibases/1000-seg.dt$problemMid
}
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel + problem.i ~ ., scales="free", labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(breaks=seq(-1000, 1000, by=2))+
  geom_point(aes(
    bases.from.middle, count, color=Sample),
    size=1, alpha=0.5,
    data=prob.data)+
  geom_segment(aes(
    norm.segStart, mean,
    xend=norm.segEnd, yend=mean),
    color="green",
    size=1,
    data=seg.dt)
print(gg)


bg.means <- seg.dt[, .(
  mean=mean(mean)
), by=.(channel, problem.i, sample, status)]
bg.means.wide <- dcast(
  bg.means, channel + problem.i + sample ~ status, value.var="mean")
bg.means.wide[, log10.ratio := log10(peak/background)]
mixture <- bg.means.wide[sample=="mixture", .(
  channel, problem.i, mixture.log10.ratio=log10.ratio)]
single <- bg.means.wide[sample!="mixture", .(
  channel, problem.i, sample, single.log10.ratio=log10.ratio)
  ][mixture, on=.(channel, problem.i)]

ggplot()+
  coord_equal()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    mixture.log10.ratio, single.log10.ratio, color=sample),
    data=single)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(paste(
    "mixture log10(peak/background),",
    "larger values for bigger peaks in mixture sample"))+
  scale_y_continuous(paste(
    "single source log10(peak/background),",
    "larger values for bigger peaks in single source sample"))

single.wide <- dcast(
  bg.means.wide,
  channel + problem.i ~ sample,
  value.var="log10.ratio")
ggplot()+
  coord_equal()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    mixture.log10.ratio, single.log10.ratio, color=sample),
    data=single)+
  geom_segment(aes(
    mixture, one.part,
    xend=mixture, yend=four.parts),
    color="red",
    alpha=0.5,
    data=single.wide)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(paste(
    "mixture log10(peak/background),",
    "larger values for bigger peaks in mixture sample"))+
  scale_y_continuous(paste(
    "single source log10(peak/background),",
    "larger values for bigger peaks in single source sample"))

mixture.m <- bg.means.wide[sample=="mixture", .(
  channel, problem.i, mixture.mean=peak)]
single.m <- bg.means.wide[sample!="mixture", .(
  channel, problem.i, sample, single.mean=peak)
  ][mixture.m, on=.(channel, problem.i)]
single.m.wide <- dcast(
  bg.means.wide,
  channel + problem.i ~ sample,
  value.var="peak")
ggplot()+
  coord_equal()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    mixture.mean, single.mean, color=sample),
    data=single.m)+
  geom_segment(aes(
    mixture, one.part,
    xend=mixture, yend=four.parts),
    color="red",
    alpha=0.5,
    data=single.m.wide)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(paste(
    "mixture peak mean,",
    "larger values for bigger peaks in mixture sample"))+
  scale_y_continuous(paste(
    "single source peak mean,",
    "larger values for bigger peaks in single source sample"))

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel + problem.i ~ ., scales="free", labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  coord_cartesian(xlim=c(-3, 3))+
  geom_point(aes(
    bases.from.middle, count, color=Sample),
    size=1, alpha=0.5,
    data=prob.data)+
  geom_segment(aes(
    norm.segStart, mean,
    xend=norm.segEnd, yend=mean),
    color="green",
    size=1,
    data=seg.dt)
print(gg)

(peak.dt <- seg.dt[status=="peak"])
setkey(peak.dt, channel, problem.i, sample, segStart, segEnd)
setkey(prob.data, channel, problem.i, sample, chromStart, chromEnd)
data.in.peaks <- foverlaps(prob.data, peak.dt, nomatch=0L)
largest.in.peaks <- data.in.peaks[, .SD[
  which.max(count)
], by=.(channel, problem.i, sample)]

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel + problem.i ~ ., labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  coord_cartesian(xlim=c(-3, 3))+
  geom_line(aes(
    bases.from.middle, relative.height, color=Sample),
    size=1, alpha=0.5,
    data=prob.data)+
  geom_point(aes(
    bases.from.middle,
    color=Sample,
    relative.height),
    data=largest.in.peaks)
print(gg)

mixture.l <- largest.in.peaks[sample=="mixture", .(
  channel, problem.i, mixture=relative.height)]
single.l <- largest.in.peaks[sample!="mixture", .(
  channel, problem.i, sample, single=relative.height)
  ][mixture.l, on=.(channel, problem.i)]
single.l.wide <- dcast(
  largest.in.peaks,
  channel + problem.i ~ sample,
  value.var="relative.height")
ggplot()+
  coord_equal()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    mixture, single, color=sample),
    data=single.l)+
  geom_segment(aes(
    mixture, one.part,
    xend=mixture, yend=four.parts),
    color="red",
    alpha=0.5,
    data=single.l.wide)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(paste(
    "mixture peak relative height,",
    "larger values for bigger peaks in mixture sample"))+
  scale_y_continuous(paste(
    "single source peak relative height,",
    "larger values for bigger peaks in single source sample"))


gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(
    channel + problem.i ~ ., labeller=label_both)+
  scale_color_manual(values=sample.colors)+
  coord_cartesian(xlim=c(-3, 3))+
  geom_line(aes(
    bases.from.middle, relative.height, color=Sample),
    size=1, alpha=0.5,
    data=prob.data)+
  geom_point(aes(
    bases.from.middle,
    color=Sample,
    relative.height),
    data=largest.in.peaks)
print(gg)

mixture.l <- largest.in.peaks[sample=="mixture", .(
  channel, problem.i, mixture=relative.height)]
single.l <- largest.in.peaks[sample!="mixture", .(
  channel, problem.i, sample, single=relative.height)
  ][mixture.l, on=.(channel, problem.i)]
single.l.wide <- dcast(
  largest.in.peaks,
  channel + problem.i ~ sample,
  value.var="relative.height")
gg <- ggplot()+
  coord_equal()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    mixture, single, color=sample),
    data=single.l)+
  geom_segment(aes(
    mixture, one.part,
    xend=mixture, yend=four.parts),
    color="grey",
    alpha=0.5,
    data=single.l.wide)+
  scale_color_manual(values=sample.colors)+
  scale_x_continuous(paste(
    "mixture peak relative height,",
    "larger values for bigger peaks in mixture sample"))+
  scale_y_continuous(paste(
    "single source peak relative height,",
    "larger values for bigger peaks in single source sample"))
png("figure-one-pair-models-joint.png", width=8, height=8, units="in", res=100)
print(gg)
dev.off()

