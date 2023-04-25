library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x


not_positioned_ags <- c("")
neut <- removeAntigens(neut, not_positioned_ags)
# Titer and antigen noise
neutBootTA <- bootstrapMap(
  neut,
  "bayesian",
  bootstrap_repeats = 500, # was 1000
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 1000,
  ag_noise_sd = 0.7,
  titer_noise_sd = 0.7,
  options = list(ignore_disconnected = TRUE)
)

save.acmap(neutBootTA, "./som/bootstrapping/neutBootTA_bayesian500_1000.ace")

# Titer noise
neutBootT <- bootstrapMap(
  neut,
  "bayesian",
  bootstrap_repeats = 500,
  bootstrap_ags = FALSE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 1000,
  ag_noise_sd = 0,
  titer_noise_sd = 0.7,
  options = list()
)

save.acmap(neutBootT, "./som/bootstrapping/neutBootT_bayesian.ace")

# Antigen noise
neutBootA <- bootstrapMap(
  neut,
  "bayesian",
  bootstrap_repeats = 500,
  bootstrap_ags = TRUE,
  bootstrap_sr = FALSE,
  reoptimize = TRUE,
  optimizations_per_repeat = 1000,
  ag_noise_sd = 0.7,
  titer_noise_sd = 0,
  options = list()
)

save.acmap(neutBootA, "./som/bootstrapping/neutBootA_bayesian.ace")

# do the blobs
neutBootTA <- read.acmap("./som/bootstrapping/neutBootTA_bayesian500_1000.ace")
neutBootTABlobs <- bootstrapBlobs(neutBootTA, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
neutBootT<- read.acmap("./som/bootstrapping/neutBootT_bayesian.ace")
neutBootTBlobs <- bootstrapBlobs(neutBootT, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
neutBootA <- read.acmap("./som/bootstrapping/neutBootA_bayesian.ace")
neutBootABlobs <- bootstrapBlobs(neutBootA, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)


# Plot the figure
png("som/bootstrapping/bayesian-bootstrap.png", width = 9, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1, 2, 3), ncol=3))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(neutBootTABlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'A', cex = 1.4)
plot(neutBootTBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'B', cex = 1.4)
plot(neutBootABlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'C', cex = 1.4)

dev.off()
