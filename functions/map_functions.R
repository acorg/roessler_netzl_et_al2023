makeMap <- function(table, baseMap = NULL, dimensions = 2,
                    nOptimisations = 500, mcb = "none", dilution_stepsize = 1, options = list()) {
  # Take a titer table, make a map. Optionally apply plotspec and re-align
  # map to an already existing map.
  m <- acmap(titer_table = table)
  
  dilutionStepsize(m) <- dilution_stepsize
  
  m <- optimizeMap(
    map                     = m,
    number_of_dimensions    = dimensions,
    number_of_optimizations = nOptimisations,
    minimum_column_basis    = mcb,
    options = options
  )
  
  agNames(m) <- rownames(table)
  srNames(m) <- colnames(table)
  
  if(!(is.null(baseMap))) {
    print("in here")
    m <- applyPlotspec(m, baseMap)
    m <- realignMap(m, baseMap)
  }
  
  ptDrawingOrder(m) <- rev(seq_len(numPoints(m)))
  
  return(m)
}

apply_color <- function(map, map_colors) {
  
  agFill(map) <- map_colors[agNames(map),]

  srOutline(map) <- map_colors[as.character(srGroups(map)),]

  return(map)
}

apply_style <- function(map){
  
  nag <- numAntigens(map)
  nsr <- numSera(map)
  N <- nsr + nag
  ptDrawingOrder(map) <- c(N:(nag+1), 1:nag)
  
  srOutlineWidth(map) <- 1
  srSize(map) <- 9
  agSize(map) <- 18
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.7)
  
  sublineages <- c("B.1.1.7+E484K", "BE.1.1", "BA.5.2.1", "CB.1", "BQ.1.3", "BQ.1.1",
    "BQ.1.18", "BR.3", "CH.1.1", "BF.7")
  agSize(map)[agNames(map) %in% sublineages] <- 15
  
  return (map)
  
}