# load required packages
rm(list = ls())
library(Racmacs)

# set seed
set.seed(100)

# ------------------------------------------ FULL MAP --------------------------------------
full_map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-full.ace")

reactivities <- rep(0, length(agNames(full_map)))
reactivities[agNames(full_map) == "P.1.1"] <- -1

full_map_p1_adj <- optimizeAgReactivity(full_map, fixed_ag_reactivities = reactivities, reoptimize = F)
full_map_p1_adj <- optimizeMap(full_map_p1_adj, number_of_dimensions = 2, number_of_optimizations = 1, options = list(ignore_disconnected = TRUE))
full_map_p1_adj <- realignMap(full_map_p1_adj, full_map)

save.acmap(full_map_p1_adj, "./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace")
##------------------------------- SINGLE EXPOSURE MAP ----------------------------------------
map_single <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")


# optimize reactivity adjustment
single_map_p1_adj <- optimizeAgReactivity(map_single, fixed_ag_reactivities =reactivities, reoptimize = F)
single_map_p1_adj <- optimizeMap(single_map_p1_adj, number_of_dimensions = 2, number_of_optimizations = 1000, options = list(ignore_disconnected = TRUE))
single_map_p1_adj <- realignMap(single_map_p1_adj, map_single)

save.acmap(single_map_p1_adj, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

# # make map with only reliably positioned sera
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")
reliable_map <- removeSera(map, srNames(map)[(grepl("G348|G353|G393|G25|G26|F647|G650", srNames(map)))])
reliable_map <- optimizeMap(reliable_map, 2, 1000, options = list(ignore_disconnected = TRUE))
reliable_map <- realignMap(reliable_map, map)
save.acmap(reliable_map, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace")

