# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

# load data
mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 
mapColors["B.1.1.7+E484K", "Color"] <- mapColors["B.1.1.7", "Color"]
mapColors["P.1.1", "Color"] <- mapColors["P.1", "Color"]

# load titer table
titer_table <- read.titerTable("./data/titer_data/titer_table.csv")
# ------------------------------------------ FULL MAP --------------------------------------
map <- makeMap(titer_table, nOptimisations = 1, dilution_stepsize = 0, options = list(ignore_disconnected = TRUE))

agFill(map) <- mapColors[agNames(map),]
agGroups(map) <- agNames(map)

sr_groups <- unlist(lapply(srNames(map), function(x) str_split(x, "_")[[1]][2]))
srGroups(map) <- factor(sr_groups, levels = srGroup_colors$SerumGroup)
srOutline(map) <- mapColors[as.character(srGroups(map)),]

map <- apply_style(map)
save.acmap(map, "./data/maps/map-OmicronI+II+III-thresholded-full.ace")

##------------------------------- SINGLE EXPOSURE MAP ----------------------------------------
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-full.ace")
single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.","BA.5 conv.","WT conv.", "CK.2.1.1 conv.")


single_exposure_sr <- srNames(map)[as.character(srGroups(map)) %in% single_exposure_sr_groups]

map_single <- subsetMap(map, sera = single_exposure_sr)

map_single <- optimizeMap(map_single, number_of_dimensions = 2, number_of_optimizations = 1000, 
                          options =  list(ignore_disconnected = TRUE))

# Copy of the same map for alignment
alignment_map <- read.acmap("../CoVAIL_trial/data/maps/alignment_map.ace")
map_single <- realignMap(map_single, alignment_map)

save.acmap(map_single, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")

