rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

# read in base map
map_single <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";", row.names = "SerumGroup")

sr_groups <- c("WT conv.", "delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT",
               "BA.1 conv." ,"BA.2 conv.", "BA.5 conv.", "CK.2.1.1 conv.")
# calculate GMTs using titertools
log_table <- adjustedTiterTable(map_single)

gmt_table <- sapply(sr_groups, function(x) {
  sub_table <- log_table[,as.character(srGroups(map_single)) == x]
 
  if(x == "BA.5 conv."){
    log2(sub_table/10)
  } else {
    apply(sub_table, 1, function(gt){
      titertools::gmt(gt, dilution_stepsize = 0)[1,1]
    })
  }
})

colnames(gmt_table) <- sr_groups

titer_gmt_table <- apply(gmt_table, 2, function(x){
  2^as.numeric(x)*10
})
rownames(titer_gmt_table) <- rownames(gmt_table)
titer_gmt_table[is.na(titer_gmt_table[,"BA.5 conv."]), "BA.5 conv."] <- "<16"
titer_gmt_table[is.na(titer_gmt_table[,"CK.2.1.1 conv."]), "CK.2.1.1 conv."] <- "<16"
titer_gmt_table[is.na(titer_gmt_table)] <- "*"

gmt_map <- makeMap(titer_gmt_table, baseMap = map_single, nOptimisations = 1000, dilution_stepsize = 0, options = list(ignore_disconnected = TRUE))
srGroups(gmt_map) <- srNames(gmt_map)
srOutline(gmt_map) <- srGroup_colors[as.character(srGroups(gmt_map)),]
gmt_map <- apply_style(gmt_map)

gmt_map <- realignMap(gmt_map, map_single)

save.acmap(gmt_map, "som/GMT_map/gmt_map.ace")
