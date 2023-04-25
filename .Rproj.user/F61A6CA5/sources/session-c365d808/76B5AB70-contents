# labelled map 
rm(list = ls())
library(Racmacs)

move_coords <- function(map, at = 2, by = -0.5){
  agCoords(map)[,at] <- agCoords(map)[,at] + by
  srCoords(map)[,at] <- srCoords(map)[,at] + by

  return(map)
}

# move all map y-coords down by 0.5
# read in map
map <- move_coords(read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace"))
map_positioned <- move_coords(read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace"))


lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
lims_no_zoom <- Racmacs:::mapPlotLims(map, sera = TRUE)

xlim_zoom <- round(lims$xlim)
ylim_zoom <- round(lims$ylim)

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

# Setup plotting function
doplot <- function(map, xlims, ylims, show_labels = TRUE) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  srSize(map) <- srSize(map) - 4
  agSize(map) <- agSize(map) - 4
  # Plot the regular map
  srOutlineWidth(map) <- 0.8
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.6)
  plot(map, xlim = xlims, 
       ylim =ylims, fill.alpha = 0.9)
  
  # Plot labels
  label_adjustments <- matrix(0, numAntigens(map), 2)
  rownames(label_adjustments) <- agNames(map)
  label_adjustments["B.1.351",] <- c(0.9, 0)
  label_adjustments["P.1.1",] <- c(-0.9, 0)
  label_adjustments["B.1.1.7+E484K",] <- c(0.9, -0.5)
  label_adjustments["BA.1",] <- c(-0.4, 0.7)
  label_adjustments["BA.2",] <- c(0, -0.6)
  label_adjustments["B.1.1.7",] <- c(0.3, -0.6)
  label_adjustments["D614G",] <- c(0, -0.5)
  label_adjustments["B.1.617.2",] <- c(0,-0.6)
  label_adjustments["BA.5.3.2",] <- c(0, 0.7)
  
  labels <- agNames(map)
  names(labels) <- agNames(map)
  labels["B.1.351"] <- "beta\n(B.1.351)"
  labels["P.1.1"] <- "gamma\n(P.1.1)"
  labels["BA.1"] <- "BA.1 omicron\n(B.1.1.529+BA.1)"
  labels["BA.2"] <- "BA.2 omicron\n(B.1.1.529+BA.2)"
  labels["B.1.617.2"] <- "delta\n(B.1.617.2)"
  labels["B.1.1.7"] <- "alpha\n(B.1.1.7)"
  labels["B.1.1.7+E484K"] <- "alpha + E484K\n(B.1.1.7+E484K)"
  labels["BA.5.3.2"] <- "BA.5.3.2 omicron\n(B.1.1.529+BA.5)"
  
  label_size <- rep(1, numAntigens(map))
  names(label_size) <- agNames(map)
  if(show_labels){
    text(
     agCoords(map) + label_adjustments,
     cex = label_size,
     label = labels,
     font = 1
    )
  }
  
  
}

png("figures/labelled_map/map_no_zoom_positioned.png", 7, 7, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map_positioned, xlim_no_zoom, ylim_no_zoom, FALSE)
dev.off()

png("figures/labelled_map/map_zoom_positioned.png", 5, 4, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map_positioned, xlim_zoom, ylim_zoom, FALSE)
dev.off()
