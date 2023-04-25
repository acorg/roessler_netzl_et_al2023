library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

sr_groups <- as.character(unique(srGroups(neut)))

n_samples <- c(1:4)

for(n in n_samples) {

    for(i in 1:10){

        sera <- c()
        for(sg in sr_groups){
            sr_names <- srNames(neut)[srGroups(neut) == sg]

            n_sera <- length(sr_names)
            if(n_sera > n){
                 sera <- c(sera, sr_names[sample.int(n_sera, n)])
             } else {
                sera <- c(sera, sr_names)
             }
        }
        
        sub_map <- subsetMap(neut, sera = sera)
        sub_map <- optimizeMap(sub_map, 2, 1000, options = list(ignore_disconnected = TRUE))
        sub_map <- realignMap(sub_map, neut)
        save.acmap(sub_map, paste0("som/subset_serum_maps/n_serum_",n, "_run_", i, ".ace"))
    }

}


for(n in n_samples){

png(paste0("som/subset_serum_maps/",n,"_sera.png"), width = 8, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:10), ncol = 5, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(i in 1:10){

    p <- read.acmap(paste0("som/subset_serum_maps/n_serum_",n,"_run_", i, ".ace"))
    
  plot(procrustesMap(p, neut, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  
}

dev.off()

}

