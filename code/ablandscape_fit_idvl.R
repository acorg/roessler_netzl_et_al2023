#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(ablandscapes)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

fit_ags <- "all" # options are "all", "sub"
sub_ags <- c('CB.1', 'BR.3', 'CH.1.1','BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF')

source("./functions/remove_reactivity_bias.R")
source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")

figure_dir <- file.path("figures", "landscapes", "gmt_landscapes")
suppressWarnings(dir.create(figure_dir, recursive = T))

ags_to_exclude <- c("")
# Read the base map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace")
map <- removeAntigens(map, ags_to_exclude)
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)

if(fit_ags == "sub"){
  padding <- 1
  ags_to_fit_lndscp <- agNames(map)[!(agNames(map) %in% sub_ags)]
} else {
  ags_to_fit_lndscp <- agNames(map)
}

# read the full map
map_orig <- read.acmap(paste0("./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace"))

sr_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";", row.names = "SerumGroup")

# set the single exposure groups
# set the single exposure groups
single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT", "BA.1 conv." ,"BA.2 conv.","BA.5 conv.",
                               "CK.2.1.1 conv.")

vacc_groups <- c("BNT/BNT", "BNT/BNT/BNT", "WT conv.")
single_exposure_sr <- srNames(map_orig)[as.character(srGroups(map_orig)) %in% single_exposure_sr_groups]

# subset the map to only multi exposure sera
map_multi <- removeSera(map_orig, single_exposure_sr)
map_multi <- removeAntigens(map_multi, ags_to_exclude)

map_long <- long_map_info(map_multi)

map_long %>%
  select(titer, ag_name, sr_name, sr_group) -> titerdata

titerdata %>%
  group_by(
    sr_group
  ) -> titerdata

titerdata %>%
  group_map(
    get_titertable
  ) -> titertables

lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable[,ags_to_fit_lndscp],
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = map,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
    
  }
)

individual_lndscp_fits <- lapply(titertables, function(titertable){
  
  lapply(1:nrow(titertable), function(sr){

    ablandscape.fit(
      titers    = titertable[sr, ags_to_fit_lndscp, drop = F],
      acmap     = map,
      bandwidth = 1,
      degree    = 1,
      method = "cone",
      error.sd = 1,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
  }
)
 
})

individual_sr_colors <- list()
for(table_nr in 1:length(titertables)){
  individual_sr_colors[[table_nr]] <- sapply(rownames(titertables[[table_nr]]), function(sr){
    "grey80"
  })
}

titertables_groups <- group_data(titerdata)

# Add impulses
titerdata %>%
  group_by(
    sr_group,
    ag_name
  ) %>%
  summarize(gmt = titertools::gmt(titer, dilution_stepsize = 0)["mean", "estimate"]) %>%
  # manually set GMT's that are lower than that to LOD2
  mutate(gmt = ifelse(gmt < log2(0.8), log2(0.8), gmt))-> gmt_data


# angle for html page
angle <- list(
  rotation = c(-1.4427, 0.0100, -0.0263), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5
  # zoom = 1.1646 # higher is more zoomed out
)

unique(titertables_groups$sr_group)

lndscp_list <- list()
data3js <- base_plot_data3js(map, lndscp_fits, agNames(map), lims, agNames(map))

opacity_val <- 1
# plot landscapes
for(srg in 1:length(unique(titertables_groups$sr_group))){

        target_rows <- srg
        lndscp_fits_t <- lndscp_fits[target_rows]
        titertables_groups_t <- titertables_groups[target_rows,]
      
        # plot idvl landscapes
        lndscp_fits_idvl <- individual_lndscp_fits[[target_rows]]
        lndscp_colors_idvl <- individual_sr_colors[[target_rows]]
        
        data3js_idvl <- plot_idvl_landscapes_from_list(data3js, lndscp_fits_idvl, lndscp_colors_idvl)
        
        lndscp_3js <- plot_landscapes_from_list(data3js_idvl, titertables_groups_t, lndscp_fits_t, map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors)
    
        lndscp <-r3js(
          lndscp_3js,
          rotation = angle$rotation,
          zoom = angle$zoom
        )

        if(unique(titertables_groups$sr_group)[srg] %in% vacc_groups){
          lndscp <-r3js(
            lndscp_3js,
            rotation = c(-1.5836, 0.0100, -0.0131),
            zoom = 2
          )
        
        
        }
        
        srg_n <- titertables_groups$sr_group[target_rows]
        srg_n <- gsub("/", "_", srg_n)
        srg_n <- gsub(" ", "", srg_n)
        save_name <- file.path(figure_dir, paste0(srg_n, "idvl_gmt_landscapes"))
        plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)
  
  }

