homologous_ags <- c(
  "mRNA1273/mRNA1273" = "D614G",
  "BNT/BNT" = "D614G",
  "AZ/BNT" = "D614G",
  "AZ/AZ" = "D614G",
  "alpha/alpha+E484K conv." = "B.1.1.7",
  "beta conv." = "B.1.351",
  "WT conv." = "D614G",
  "delta conv." = "B.1.617.2",
  "BA.1 conv." = "BA.1",
  "BA.2 conv." = "BA.2",
  "BA.5 conv." = "BA.5.3.2",
  "CK.2.1.1 conv." = "BA.5.2.1"
)

sr_pretty <- data.frame(
  row.names = c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', "CK.2.1.1 conv.",
  'BNT/BNT/BNT', 'Vacc+BA.1', 'Vacc+BA.2', 'Vacc+BA.5', 'BA.1 reinf.', 'BA.2 reinf.', 'BA.2+BA.5 reinf.', 'Vacc+BA.1 reinf.', 'Vacc+BA.2 reinf.', 'AZ/AZ+delta', 'BNT/BNT+delta'),
  val = c('mRNA-1273/mRNA-1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"First wave conv.", 'alpha conv.', 'beta conv.', 'delta conv.', 'BA.1 omicron conv.', 'BA.2 omicron conv.', 'BA.5 omicron conv.', "CK.2.1.1 conv.",
  'BNT/BNT/BNT', 'Vacc+BA.1', 'Vacc+BA.2', 'Vacc+BA.5', 'BA.1 reinf.', 'BA.2 reinf.', 'BA.2+BA.5 reinf.', 'Vacc+BA.1 reinf.', 'Vacc+BA.2 reinf.', 'AZ/AZ+delta', 'BNT/BNT+delta')
)

sr_group_homologous_antigen <- function(data){
  homologous_ags <- unlist(lapply(levels(data$sr_group), function(x) {
    name <- strsplit(x, " ")[[1]][1]
    if(TRUE %in% grepl("WT|mRNA1273|BNT162b2|ChAdOx1-S", name)) {
      name <- "D614G"
    } else if(name == "B.1.526+E484K"){
      name <- "B.1.526_E484K"  
    } else {
      name
    }
  }))
  names(homologous_ags) <- levels(data$sr_group)
  
  return(homologous_ags)
}

# calculate fold change
calculate_fold_change_from_homolgous <- function(data) {
  #homologous_ags <- sr_group_homologous_antigen(data)
  
  data$fold_change <- unlist(lapply(1:nrow(data), function(x) {
    
    sr_group <- as.character(data[x, "sr_group"]$sr_group)
    

    hom_ag <- homologous_ags[sr_group]
    print(hom_ag)
    hom_gmt <- data[data$ag_name == hom_ag & data$sr_group == sr_group, "logtiter"]$logtiter
    hom_gmt <- (2^hom_gmt*10)
  
    this_titer <- (2^data[x, "logtiter"]$logtiter*10)
    
    if(length(hom_gmt) < 1 | length(this_titer) <1 ) {
      ""
    } else if(is.na(hom_gmt) | is.na(this_titer)) {
      ""
    } else if(this_titer > hom_gmt) {
      paste0("+", round(this_titer/hom_gmt,1))
    } else if(this_titer == hom_gmt) {
      "1"
    } else {
      paste0("-", round(hom_gmt/this_titer,1))
    }
    
  }))
  
  return(data)
}


# calculate relative value from different data frame
calculate_realtive_titer <- function(data_big, data_small, column_name = "gmt", new_col_name = "relative_titer"){
  
  data_big_gmts <- data_big[match(interaction(data_small$ag_name, data_small$sr_group), interaction(data_big$ag_name, data_big$sr_group)),column_name]
  data_big_gmts <- (2^data_big_gmts*10)$gmt
  data_small_gmts <- (2^data_small[,column_name]*10)$gmt
  
 # return(data_small_gmts/data_big_gmts)
  fold_change <- unlist(lapply(c(1:length(data_big_gmts)), function(x){
    
    if(is.na(data_small_gmts[x])){
      ""
    } else {
      paste0(round(data_small_gmts[x]/data_big_gmts[x],1))
    }
  }))
  
  data_big[new_col_name] <- ""
  data_big[match(interaction(data_small$ag_name, data_small$sr_group), interaction(data_big$ag_name, data_big$sr_group)),new_col_name] <- fold_change

  return(data_big)  
}

# gmt and fold change label 
get_fold_change_label <- function(data, ymax) {
  
  data$titer <- round(2^data$logtiter*10)
  data$titer[is.na(data$titer)] <- ""
  if(length(data$label) > 0) {
    data$label <- paste0(data$label, "\n", data$fold_change)
  }else {
    data$label <- paste0(data$titer, "\n", data$fold_change)
  }

  data$y <- ymax - 0.5
  
  return(data)
}

# gmt and fold change label 
get_gmt_label <- function(data, ymax) {
  
  data$titer <- round(2^data$logtiter*10)
  data$titer[is.na(data$titer)] <- ""
  data$y <- ymax - 0.5
  data$label <- paste0(data$titer)
  
  return(data)
}



sr_group_gmt_calc <- function(plotdata, thresh, half_thresh = FALSE) {
  plotdata$below_thresh <- plotdata$titer == paste0("<",thresh)
  
  if(!(half_thresh)) {
    # Get gmts by serum group
    plotdata %>%
      group_by(
        sr_group,
        ag_name
      ) %>%
      summarise(
        logtiter = meantiter::mean_titers(
          titers = titer,
          method = "bayesian",
          dilution_stepsize = 0
        )$mean,
        titer = 2^logtiter*10,
        all_below_thresh = !(FALSE %in% below_thresh)
      ) -> sr_group_gmt_plotdata
    
    sr_group_gmt_plotdata$titer <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
      if(sr_group_gmt_plotdata$all_below_thresh[x]) {
        paste0("<",thresh)
      } else {
        sr_group_gmt_plotdata$titer[x]
      }
    }))
    
    sr_group_gmt_plotdata$logtiter <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
      if(sr_group_gmt_plotdata$all_below_thresh[x]) {
        log2(thresh/10/2)
      } else {
        sr_group_gmt_plotdata$logtiter[x]
      }
    }))
  } else {
    plotdata$titer[plotdata$below_thresh] <- thresh/2
    
    plotdata %>%
      group_by(
        sr_group,
        ag_name
      ) %>%
      summarise(
        logtiter = meantiter::mean_titers(
          titers = titer,
          method = "bayesian",
          dilution_stepsize = 0
        )$mean,
        titer = 2^logtiter*10,
        all_below_thresh = !(FALSE %in% below_thresh)
      ) -> sr_group_gmt_plotdata
    
  }
 
  
  return(sr_group_gmt_plotdata)
}

do_titer_plot_fc_label <- function(map, facet_col = 3, thresh = 20, fc_label =TRUE, adj_titers = TRUE,
  remove_idvl_effects = FALSE,
  highlight_samples = NULL) {
  # Set the antigen x axis order
  ag_order <- agNames(map)
  
  # Get long info
  plotdata <- long_map_info(map)
  
  if(adj_titers){
    plotdata$logtiter <- plotdata$logtiter_adjusted
    plotdata$titer <- plotdata$titer_adjusted
  }
  
  ymin <- -0.5
  ymax <- 11.5
  if(remove_idvl_effects){
    plotdata <- adjust_individual_effect(plotdata)
    ymin <- -0.5
    ymax <- 10
  }
  # Get gmts by serum group
  sr_group_gmt_plotdata <- sr_group_gmt_calc(plotdata, thresh)
  
  # Clamp titers to <thresh and indicate they are <thresh
  plotdata$lessthan20 <- FALSE
  plotdata$lessthan20[plotdata$logtiter < log2(thresh/10)] <- TRUE
  plotdata$logtiter[plotdata$titertype == 2 | plotdata$logtiter < log2(thresh/10)] <- log2(thresh/10)/2
  
  if(fc_label) {
    # add label
    sr_group_gmt_plotdata %>% 
      calculate_fold_change_from_homolgous() %>%
      get_fold_change_label(.,ymax = 11.5) -> sr_group_gmt_plotdata
  } else {
    sr_group_gmt_plotdata %>% 
      get_gmt_label(.,ymax = 11.5) -> sr_group_gmt_plotdata
  }
  

  
  # Clamp sr group gmts to <20 and indicate they are <20
 # if(!remove_idvl_effects){
    sr_group_gmt_plotdata$lessthan20 <- FALSE
    sr_group_gmt_plotdata$lessthan20[sr_group_gmt_plotdata$logtiter < log2(thresh/10)] <- TRUE
    sr_group_gmt_plotdata$logtiter[sr_group_gmt_plotdata$logtiter < (log2(thresh/10)/2)] <- log2(thresh/10)/2
#  }
  
 # facet_labels <- paste(sr_pretty[unique(plotdata$sr_group), "vals"], "sera")
 # names(facet_labels) <- unique(plotdata$sr_group)
  
facet_labeller <- function(x){
    paste0(sr_pretty[x,"val"], " sera")
  }

#plotdata$ag_name <- factor(plotdata$ag_name, levels = 
#        agNames(map) <- c("D614G", "Alpha", "Alpha+E484K", "Gamma", "Delta","BA.2 omicron", "BA.5 omicron","Beta", "BA.1 omicron")
#)
  # Do the plot
  plotdata %>%
    ggplot(
      aes(
        x = ag_name,
        y = logtiter,
        color = sr_group,
        # shape = lessthan20,
        # size = lessthan20,
        fill = sr_group
      )
    ) + 
    geom_line(
      aes(group = sr_name),
      alpha = 0.4,
      size = 0.8
    ) + 
    geom_point(
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes(group = sr_group),
      size = 1.2,
      alpha = 1
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 2.5,
      alpha = 1
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(
      values = srGroupOutline(map)
    ) +
    scale_fill_manual(
      values = srGroupOutline(map)
    ) +
    scale_y_titer(
      ymin = 1,
      ymax = 11
    ) + 
    scale_x_discrete(
      limits = ag_order,
      expand = expansion(add = 0)
    ) +
    # scale_shape_manual(
    #   values = c(
    #     "FALSE" = 21,
    #     "TRUE" = 25
    #   )
    # ) +
    # scale_size_manual(
    #   values = c(
    #     "FALSE" = 1.2,
    #     "TRUE" = 0.8
    #   )
  # ) +
  facet_wrap(
    vars(sr_group),
    labeller =  as_labeller(facet_labeller),
    ncol = facet_col
  ) + 
    labs(
      x = "Antigen variant"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(ymin, ymax), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none"
    ) -> gp
  
 
  # Annotate region below detectable
  
  gp <- gp +
    annotate(
      "tile",
      x = agNames(map),
      y = 0,
      height = 1+log2(thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    )
  
  if(!is.null(highlight_samples)){
    sub_samples <- plotdata[grepl(paste(highlight_samples, collapse = "|"), plotdata$sr_name),]

    gp <- gp + 
      geom_line(data = sub_samples, aes(group = sr_name),
      alpha = 0.4,
      size = 0.8, color = "red")
  }
  
  
  # Annotate colors for each antigen
  for (n in seq_len(numAntigens(map))) {
    gp <- gp +
      annotate(
        "tile",
        x = agNames(map)[n],
        y = ymin-0.25,
        height = 1,
        fill = agFill(map)[n],
        color = NA
      ) 
  }
  
  # Label fold change from homologous
  if(fc_label){
    gp <- gp +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y-0.5, label = label),
              color = "black",
              size = 2
    )
  }
  
  return(gp)
  
}



do_ag_titer_plot_fc_label <- function(map, facet_col = 3, thresh = 20, fc_label =TRUE, adj_titers = TRUE) {
  # Set the antigen x axis order
  ag_order <- agNames(map)
  
  # Get long info
  plotdata <- long_map_info(map)
  
  if(adj_titers){
    plotdata$logtiter <- plotdata$logtiter_adjusted
    plotdata$titer <- plotdata$titer_adjusted
  }
 
  # Get gmts by serum group
  sr_group_gmt_plotdata <- sr_group_gmt_calc(plotdata, thresh)
  
  # Clamp titers to <thresh and indicate they are <thresh
  plotdata$lessthan20 <- FALSE
  plotdata$lessthan20[plotdata$logtiter < log2(thresh/10)] <- TRUE
  plotdata$logtiter[plotdata$titertype == 2] <- log2(thresh/10)/2
  
  
  if(fc_label) {
    # add label
    sr_group_gmt_plotdata %>% 
      calculate_fold_change_from_homolgous() %>%
      get_fold_change_label(.,ymax = 11.5) -> sr_group_gmt_plotdata
  } else {
    sr_group_gmt_plotdata %>% 
      get_gmt_label(.,ymax = 11.5) -> sr_group_gmt_plotdata
  }
  

  
  # Clamp sr group gmts to <20 and indicate they are <20
  sr_group_gmt_plotdata$lessthan20 <- FALSE
  sr_group_gmt_plotdata$lessthan20[sr_group_gmt_plotdata$logtiter < log2(thresh/10)] <- TRUE
  sr_group_gmt_plotdata$logtiter[sr_group_gmt_plotdata$logtiter < (log2(thresh/10)/2)] <- log2(thresh/10)/2
  
 # facet_labels <- paste(sr_pretty[unique(plotdata$sr_group), "vals"], "sera")
 # names(facet_labels) <- unique(plotdata$sr_group)
  
facet_labeller <- function(x){
    paste0(sr_pretty[x,"val"], " sera")
  }

#plotdata$ag_name <- factor(plotdata$ag_name, levels = 
#        agNames(map) <- c("D614G", "Alpha", "Alpha+E484K", "Gamma", "Delta","BA.2 omicron", "BA.5 omicron","Beta", "BA.1 omicron")
#)
 plotdata$ag_name <- factor(plotdata$ag_name, levels = ag_order)
  sr_group_gmt_plotdata$ag_name <- factor(sr_group_gmt_plotdata$ag_name, levels = ag_order)
  # Do the plot
  plotdata %>%
    ggplot(
      aes(
        x = sr_group,
        y = logtiter,
        color = sr_group,
        # shape = lessthan20,
        # size = lessthan20,
        fill = sr_group
      )
    ) + 
    geom_line(
     # aes(group = ag_n),
      alpha = 0.4,
      size = 0.8
    ) + 
    geom_point(
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
    #  aes(group = sr_group),
      size = 1.2,
      alpha = 1
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 2.5,
      alpha = 1
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
     scale_color_manual(
       values = srGroupOutline(map)
    ) +
    scale_fill_manual(
       values = srGroupOutline(map)
    ) +
    scale_y_titer(
      ymin = 1,
      ymax = 11
    ) + 
    # scale_x_discrete(
    #   limits = ag_order,
    #   expand = expansion(add = 0)
    # ) +
    # scale_shape_manual(
    #   values = c(
    #     "FALSE" = 21,
    #     "TRUE" = 25
    #   )
    # ) +
    # scale_size_manual(
    #   values = c(
    #     "FALSE" = 1.2,
    #     "TRUE" = 0.8
    #   )
  # ) +
  facet_wrap(
    vars(ag_name),
   # labeller =  as_labeller(facet_labeller),
    ncol = facet_col
  ) + 
    labs(
      x = "Serum group"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(0, 11.5), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none"
    ) -> gp
  
 
  # Annotate region below detectable
  gp <- gp +
    annotate(
      "tile",
      x = unique(as.character((srGroups(map)))),
      y = 0,
      height = 1+log2(thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    )
  
  sr_groups_char <- unique(as.character((srGroups(map))))
  # Annotate colors for each antigen
  # for (n in seq_len(sr_groups_char)) {
  #   gp <- gp +
  #     annotate(
  #       "tile",
  #       x = sr_groups_char[n],
  #       y = -0.75,
  #       height = 1,
  #       fill = agFill(map)[n],
  #       color = NA
  #     ) 
  # }
  
  # Label fold change from homologous
  gp <- gp +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = sr_group, y = y-1, label = label),
              color = "black",
              size = 2
    )
  
  return(gp)
  
}

# Function for accounting for a sample's indivdual affect
# calculate GMT per sample, subtract serum group GMT from that to get reactivity bias
# subtract sample's reactivity bias from sample titrations
adjust_individual_effect <- function(data, dil_stepsize = 0){
  
  # calculate gmt per sample
  data %>%
    group_by(sr_name) %>%
    mutate(gmt_sample = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  # calculate gmt per serum group
  titers_variation_adjusted %>%
    group_by(sr_group) %>%
    mutate(gmt_sr_group = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  
  # adjust serum reactivity bias
  titers_variation_adjusted %>%
    mutate(reactivity_bias = gmt_sample - gmt_sr_group,
           logtiter = logtiter - reactivity_bias,
          # logtiter = logtiter -gmt_sample,
           titer = 2^logtiter*10) -> titers_variation_adjusted
  
  return(titers_variation_adjusted)
  
}