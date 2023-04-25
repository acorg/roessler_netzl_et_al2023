library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x 

labels <- data.frame(
  row.names = agNames(neut),
  val = LETTERS[c(1:length(agNames(neut)))]
)

paste(labels$val, ag_pretty$val, collapse = ", ")

ag_pretty <- data.frame(
  row.names = c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", 
                'CB.1', 'BR.3', 'CH.1.1', 'BA.5.3.2', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF'),
  val = c('D614G', 'alpha', 'alpha+E484K', 'gamma', 'beta', 'delta', 'BA.1 omicron', 'BA.2 omicron',
          'CB.1', 'BR.3', 'CH.1.1', 'BA.5 omicron', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF')
)

png("som/bootstrapping/bootstrapping-antigens.png", width = 10, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:nrow(labels)), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(ag in agNames(neut)){
  print(ag)
  
  newMap <- removeAntigens(neut, c(ag))
  newMap <- optimizeMap(
    map                     = newMap,
    number_of_dimensions    = 2,
    number_of_optimizations = 1000,
    minimum_column_basis    = "none",
    options = list(ignore_disconnected = TRUE)
  )
  
  newMap <- realignMap(newMap, neut)
  
  save.acmap(map = newMap, filename = paste0("./som/bootstrapping/wo_ag_",ag,".ace"))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = ag, cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[ag, ], cex = 1.2)
}

dev.off()

png("som/bootstrapping/bootstrapping-antigens.png", width = 10, height = 12, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:24), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(ag in agNames(neut)){
 
  newMap <- read.acmap(filename = paste0("./som/bootstrapping/wo_ag_",ag,".ace"))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = ag_pretty[ag,], cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.6, ylim_no_zoom[2]-0.6, labels[ag, ], cex = 1.2)
}

dev.off()
