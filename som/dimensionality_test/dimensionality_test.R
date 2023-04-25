rm(list = ls())
library(Racmacs)
library(ggplot2)
set.seed(100)
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace")


# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 1000,
  replicates_per_dimension = 100,
  options = list()
)

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))

saveRDS(df, "./som/dimensionality_test/dimension_test_result.rds")

# do the plot
df <- readRDS("./som/dimensionality_test/dimension_test_result.rds")
ggplot(data=df, aes(x=dimensions, y=rmse)) +
  geom_line()+
  geom_point() +
  theme_bw() +
  ylim(c(0,3.5))+
  xlab('Dimension') +
  ylab('Mean RMSE of detectable titers') +
  theme(strip.background = element_blank()) ->dp


png("./som/dimensionality_test/dimension_test.png", 3, 3, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
dp
dev.off()

# Plot a 3D map
map3D <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 3,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none",
  options = list(ignore_disconnected = TRUE)
)

map3D <- applyPlotspec(map3D, map)
map3D <- realignMap(map3D, map)

p <- procrustesMap(map3D, map, sera=FALSE)
agSize(p) <- 6
srSize(p) <- 4
Racmacs::view(p)
