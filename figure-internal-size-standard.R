source("packages.R")
fsa_db <- fread("fsa_db.csv")
fsa_present <- fsa_db[file.exists(fsa)]
fsa_present[, type := data.table::fcase(
  control != "", "control",
  sample.identifier != "", "single",
  identifiers.str != "", "mixture",
  default="unknown")]
fsa_present[, file.i := 1:.N]
mix.dt <- fsa_present[type=="mixture"]
##one mixture.
one.mix <- mix.dt[1]
ip.fsa <- fsa_present[kit=="IP" & type=="single"][2]

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

ip.dt <- fsa2dt(ip.fsa$fsa)

## looks like the first two are missing, so the first peak is 60.
ip.sizes <- c(
  ##20, 40,
  60, 80, 100,
  114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314,
  320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540,
  560, 580, 600)
length(ip.sizes)

ip.thresh <- 2400
ip.after <- ip.dt[Time>ip.thresh]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=ip.after)
print(gg)
## looks like 34 peaks.

ip.before <- ip.dt[Time<ip.thresh]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=ip.before)
print(gg)
##?

ip.standard.after <- ip.after[channel==105]
div <- 1000
ip.standard.after[, relative := Time %% div]
ip.standard.after[, panel := Time %/% div]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(panel ~ ., labeller=label_both)+
  geom_line(aes(
    relative, RFU),
    data=ip.standard.after)
print(gg)

## The Internal Lane Standard 600 (ILS 600)
## consists of 22 bands ranging in size from 60bp to 600bp. Fragments of
## 60–200bp are spaced at 20bp intervals, fragments of 200–500bp are
## spaced every 25 bases, and fragments of 500–600bp are spaced every 50
## bases. Fragments that are multiples of 100 bases have fluorescence
## intensities approximately twice that of other fragments to simplify
## size assignment.
ILS600 <- c(
  seq(60, 200, by=20),
  seq(200, 500, by=25)[-1],
  seq(500, 600, by=50)[-1])
n.peaks <- length(ILS600)
pp16.fsa <- fsa_present[kit=="PP16" & type=="single"][1]
pp16.dt <- fsa2dt(pp16.fsa$fsa)

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both)+
  geom_line(aes(
    Time, RFU),
    data=pp16.dt)
print(gg)


pp16.thresh <- 2500
pp16.after <- pp16.dt[Time>pp16.thresh]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both, scales="free")+
  geom_line(aes(
    Time, RFU),
    data=pp16.after)
print(gg)

pp16.after[, non.neg := RFU-min(RFU)]
pp16.after[, non.neg := ifelse(RFU<0, 0, RFU)]
range(pp16.after$non.neg)
str(pp16.after)
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(channel ~ ., labeller=label_both, scales="free")+
  geom_line(aes(
    Time, non.neg),
    data=pp16.after)
print(gg)

gg+xlim(5900, 6600)

gg+xlim(4000, 5000)

pdir <- "figure-internal-size-standard-data"
if(FALSE){
  unlink(pdir, recursive=TRUE)
}
dir.create(pdir)
pp16.standard.after <- pp16.after[channel==4]
rle.list <- rle(pp16.standard.after$non.neg)
table(rle.list$lengths)
edge.dt <- data.table(
  Time=cumsum(c(pp16.standard.after$Time[1]-1L, rle.list$lengths)))
coverage.dt <- edge.dt[, data.table(
  chrom = "chrUnknown",
  chromStart = Time[-.N],
  chromEnd = Time[-1],
  count = rle.list$values)]
coverage.bedGraph <- file.path(pdir, "coverage.bedGraph")
PeakSegDisk::writeBedGraph(coverage.dt, coverage.bedGraph)

fit <- PeakSegDisk::sequentialSearch_dir(pdir, n.peaks)

peak.dt <- fit$segments[status=="peak"][order(chromStart), .(
  mean,
  peakStart=chromStart+0.5,
  peakEnd=chromEnd+0.5,
  bases=ILS600)]
setkey(peak.dt, peakStart, peakEnd)
pp16.standard.after[, Time0 := Time]
setkey(pp16.standard.after, Time, Time0)
standard.in.peaks <- foverlaps(pp16.standard.after, peak.dt, nomatch=0L)
max.dt <- standard.in.peaks[, .SD[RFU==max(RFU)], by=bases]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_line(aes(
    Time, non.neg),
    data=pp16.standard.after)+
  geom_segment(aes(
    chromStart+0.5, mean,
    xend=chromEnd+0.5, yend=mean),
    color="green",
    data=fit$segments)+
  geom_text(aes(
    Time, RFU, label=bases),
    data=max.dt,
    vjust=-0.5,
    color="green")
print(gg)

ggplot()+
  geom_text_repel(aes(
    peakEnd-peakStart, mean, label=bases),
    data=peak.dt)+
  geom_point(aes(
    peakEnd-peakStart, mean),
    shape=1,
    data=peak.dt)+
  xlab("peak width")+
  ylab("peak height")

peak.dt[, peakBases := peakEnd-peakStart]
expand.bases <- median(peak.dt$peakBases)*1
peak.dt[, regionStart := peakStart-expand.bases]
peak.dt[, regionEnd := peakEnd+expand.bases]
peak.dt[, peak.i := 1:.N]
setkey(peak.dt, regionStart, regionEnd)
coverage.zoom <- foverlaps(peak.dt, pp16.standard.after, nomatch=0L)
segs.dt <- fit$segments[order(chromStart), data.table(
  segStart=chromStart+0.5,
  segEnd=chromEnd+0.5,
  mean)]
setkey(segs.dt, segStart, segEnd)
segs.zoom <- foverlaps(peak.dt, segs.dt, nomatch=0L)
segs.zoom[, showStart := ifelse(
  segStart < regionStart, regionStart, segStart)]
segs.zoom[, showEnd := ifelse(
  regionEnd < segEnd, regionEnd, segEnd)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  ##facet_grid(. ~ bases, scales="free", space="free")+
  facet_wrap("bases", scales="free_x")+
  geom_point(aes(
    Time, non.neg),
    shape=1,
    data=coverage.zoom)+
  geom_segment(aes(
    showStart, mean,
    xend=showEnd, yend=mean),
    color="green",
    data=segs.zoom)+
  geom_vline(aes(
    xintercept=showEnd),
    color="green",
    data=segs.zoom[regionEnd!=showEnd])+
  scale_x_continuous(breaks=pp16.standard.after[, seq(0, max(Time), by=20)])+
  geom_point(aes(
    Time, RFU),
    data=max.dt)
print(gg)
peak.dt[bases==200]

max.dt[, Time := chromEnd]
lin.fit <- lm(bases ~ Time, max.dt)
lin.fit.dt <- data.table(t(coef(lin.fit)))
pred.dt <- max.dt[, data.table(bases=seq(min(bases), max(bases)))]
m <- lin.fit.dt$Time
b <- lin.fit.dt[["(Intercept)"]]
## bases = m*Time + b so Time = (bases - b)/m
pred.dt[, Time := (bases-b)/m]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    Time, bases, color=what),
    data=data.table(max.dt, what="data"))+
  geom_abline(aes(
    slope=Time, intercept=`(Intercept)`, color=what),
    data=data.table(lin.fit.dt, what="model"))+
  geom_point(aes(
    Time, bases, color=what),
    shape=1,
    data=data.table(pred.dt, what="model"))
print(gg)

## log fit looks worse.
max.dt[, log.Time := log(Time)]
max.dt[, log.bases := log(bases)]
log.fit <- lm(log.bases ~ log.Time, max.dt)
log.fit.dt <- data.table(t(coef(log.fit)))
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(log.Time, log.bases), data=max.dt)+
  geom_abline(aes(
    slope=log.Time, intercept=`(Intercept)`),
    data=log.fit.dt)
print(gg)

## Poisson model.
pois.fit <- glm(bases ~ Time, max.dt, family=poisson(link="identity"))
max.dt[, pois.pred := predict(pois.fit)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    Time, bases, color=what),
    data=data.table(max.dt, what="data"))+
  geom_line(aes(
    Time, pois.pred, color=what),
    data=data.table(max.dt, what="model"))
print(gg)
