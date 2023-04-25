rm(list = ls())

library(Racmacs)
library(tidyverse)
library(patchwork)
library(titertools)
set.seed(100)

# define titerplot function
source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")


sr_pretty <- data.frame(
  row.names = c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', "BA.5 conv.", "CK.2.1.1 conv."),
  val = c('2xmRNA-1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"Anc. virus conv.", 'alpha conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', "BA.5. conv.", "CK.2.1.1 conv.")
)

# create table of titrated sera
map_positioned <- read.acmap("data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace")
data_long <- long_map_info(map_positioned)

data_long %>%
  filter(titer == "*")

data_long %>%
  filter(titer != "*") %>%
  group_by(ag_name, sr_group) %>%
  summarize(nr_sera = length(unique(sr_name))) %>%
  pivot_wider(names_from = sr_group, values_from = nr_sera) -> nr_sera_table

ag_order <- c("D614G", "B.1.1.7", "B.1.351", "P.1.1", "B.1.617.2", 
              "BA.1", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.3.2",
              "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18",
              "XBB.1", "XBB.1.5", "XBF")

nr_sera_table$ag_name <- factor(nr_sera_table$ag_name, c("D614G", "B.1.1.7", "B.1.1.7+E484K", "B.1.351", "P.1.1", "B.1.617.2", 
                                                         "BA.1", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.3.2",
                                                         "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18",
                                                         "XBB.1", "XBB.1.5", "XBF"))

nr_sera_table <- nr_sera_table[order(nr_sera_table$ag_name),]

colnames(nr_sera_table) <- gsub("WT", "Anc. virus", colnames(nr_sera_table))
colnames(nr_sera_table) <- gsub("ag_name", "Variant", colnames(nr_sera_table))
write.csv(nr_sera_table, "data/titer_data/nr_sera_in_map.csv", row.names = FALSE)

#----------- start fold change calculation
map <- read.acmap("data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace")
data_long <- long_map_info(map)

# get the data to long format with the titertable estimates

sr_groups <- unique(data_long$sr_group)
ags <- unique(data_long$ag_name)

old_ags <- c("D614G", "BA.1", "BA.5.3.2")
new_ags <- c("B.1.1.7", "B.1.351", "B.1.617.2", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18", "XBB.1", "XBB.1.5", "XBF")

all_ags <- c(old_ags, new_ags)

calc_fc_multiple_ags <- function(data_long, sr_groups, old_ags, new_ags){

  fc_df <- data.frame(sr_group = c(),
  ag1 = c(),
  ag2 = c(),
  fc = c(),
  fc_lower = c(),
  fc_upper = c())

# calculate fold changes per sr group and ag
# have to optimise ag selection to not do everything double
for(sr in sr_groups){
  for(ag1 in old_ags){
    for(ag2 in new_ags){

      titers1 <- data_long %>% filter(sr_group == sr) %>%
        filter(ag_name == ag1) %>%
        pull(titer_adjusted)

       titers2 <- data_long %>% filter(sr_group == sr) %>%
        filter(ag_name == ag2) %>%
        pull(titer_adjusted)

      diff <- log2diff(titers1, titers2, dilution_stepsize = 0)

      temp_df <- data.frame(sr_group = sr,
        ag1 = ag1,
        ag2 = ag2,
        fc = diff["mean", "estimate"],
        fc_lower = diff["mean", "lower"],
        fc_upper = diff["mean", "upper"])
      
      fc_df <- rbind(fc_df, temp_df)
    }
  }
}

  fc_df$ag1 <- factor(fc_df$ag1, levels = ags)
  fc_df$ag2 <- factor(fc_df$ag2, levels = ags)

  return(fc_df)

}


fc_df <- calc_fc_multiple_ags(data_long, sr_groups, old_ags, all_ags)
write.csv(fc_df, "som/titer_fold_changes/230412_titertools_fold_changes_3_ags.csv")

fc_df <- read.csv("som/titer_fold_changes/230412_titertools_fold_changes_3_ags.csv")

sr_labeller <- function(x){
  sr_pretty[as.character(x), "val"]
}

target_groups <- as.character(sr_groups[grepl("3x", sr_groups)])

## FC boxplots
boxplot_colors <- srOutline(map)
names(boxplot_colors) <- as.character(srGroups(map))
boxplot_colors <- boxplot_colors[target_groups]

fc_df$sr_group <- factor(fc_df$sr_group, levels = sr_groups)

# get order of 2nd antigen by fold drop from D614G in biv BA.1
fc_order <- fc_df %>%
  filter(sr_group == "3x Vacc + bivalent BA.1 boost") %>%
  filter(ag1 == "D614G")

fc_order <- fc_order[order(fc_order$fc, decreasing = TRUE),]

# change to Janine's preferred order
fc_order$ag2 <- c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2", 
                  "BA.1", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.3.2",
                  "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18",
                  "XBB.1", "XBB.1.5", "XBF")

fc_df$ag2 <- factor(as.character(fc_df$ag2), levels = as.character(fc_order$ag2))


# get fraction below detection threshold
# get nr of below threshold values
data_long %>%
  group_by(sr_group, ag_name) %>%
  summarize(n_sera = length(titer),
            below_thresh = length(titer[grepl("<", titer)]),
            frac_thresh = below_thresh/n_sera) -> below_thresh

below_thresh$sr_group <- as.character(below_thresh$sr_group)
below_thresh$ag_name <- as.character(below_thresh$ag_name)


fc_df$below_thresh <- sapply(1:nrow(fc_df), function(x){
  
  srg <- fc_df$sr_group[x]
  ag <- fc_df$ag2[x]
  
  below_thresh$frac_thresh[below_thresh$sr_group == srg & below_thresh$ag_name == ag]
  
})

fc_df$ag1 <- factor(fc_df$ag1, levels = old_ags)

# change names of variants on x axis
ag_names <-  c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2", 
               "BA.1", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.3.2",
               "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18",
               "XBB.1", "XBB.1.5.1", "XBF.3")
names(ag_names) <- fc_order$ag2

fc_df$ag2 <- factor(ag_names[as.character(fc_df$ag2)], levels = ag_names)
fc_df %>%
  filter(sr_group %in% target_groups) %>%
  filter(ag1 %in% old_ags) %>%
  filter(!ag2 %in% c("B.1.1.7", "BE.1.1", "BA.5.2.1")) %>%
  ggplot(aes(x = ag2, y = fc, color = sr_group, group = sr_group, fill = sr_group)) + 
  geom_pointrange(aes(ymin = fc_lower, ymax = fc_upper), position = position_dodge(width = 0.5), alpha = 1) + 
  geom_line(alpha = 0.7, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = boxplot_colors, name = c(""),
                     labels = c("BA.1 biv. / N-", "BA.1 biv. / N+", "BA.4/5 biv. / N-", "BA.4/5 biv. / N+")) + 
  geom_point(aes(alpha = below_thresh), shape = 21, fill = "white", position = position_dodge(width = 0.5), size = 1.5) +
  scale_alpha(name = "Fraction <LOD") +
  guides(fill = "none") +
  facet_wrap(~ag1,
             nrow = 3) + 
  scale_y_continuous(
    breaks = seq(-20, 2, 2),
    labels = function(x) ifelse(x<0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")),
    name = "Fold change from "
  ) +
  xlab("Variant") +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") + 
  titerplot_theme() + 
  theme(legend.position = "right") -> fc_vacc_sr_groups

ggsave("som/titer_fold_changes/fc_3_vacc_sera.png", fc_vacc_sr_groups, width = 9, height = 8, dpi = 300)

#-------------------  now fold change for map sr groups
fc_map <- fc_df %>%
  filter(sr_group %in% rownames(sr_pretty)) %>%
  filter(ag1 == "D614G")

# get order 
fc_order_map <- fc_df %>%
  filter(sr_group == "BNT/BNT") %>%
  filter(ag1 == "D614G")
fc_order_map <- fc_order_map[order(fc_order_map$fc, decreasing = TRUE), ]

fc_order$ag2 <- c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2", 
                  "BA.1", "BA.2", "CB.1", "BR.3", "CH.1.1", "BA.5.3.2",
                  "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18",
                  "XBB.1", "XBB.1.5", "XBF")

fc_map$ag2 <- factor(as.character(fc_map$ag2), levels = as.character(fc_order_map$ag2))
fc_map$sr_group <- factor(fc_map$sr_group, levels = rownames(sr_pretty))

# colours
map_colors <- srOutline(map)
names(map_colors) <- as.character(srGroups(map))
map_colors <- map_colors[rownames(sr_pretty)]

fc_map %>%
  ggplot(aes(x = ag2, y = fc, color = sr_group, group = sr_group)) + 
  geom_pointrange(aes(ymin = fc_lower, ymax = fc_upper), position = position_dodge(width = 0.5), alpha = 1) + 
  geom_line(alpha = 0.7, position = position_dodge(width = 0.5)) + 
  geom_point(aes(alpha = below_thresh),shape =21, fill = "white", position = position_dodge(width = 0.5), size = 1.5) + 
  scale_color_manual(values = map_colors, name = "") + 
  facet_wrap(~sr_group,
             labeller = as_labeller(sr_labeller)) + 
  scale_alpha(name = "Fraction <LOD") +
  guides(color = "none") +
  scale_y_continuous(
    breaks = seq(-20, 12, 2),
    limits = c(-12,12),
    oob=scales::squish,
    labels = function(x) ifelse(x<0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")),
    name = "Fold change from D614G"
  ) +
  xlab("Variant") +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") + 
  titerplot_theme() + 
  theme(legend.position = "top") -> fc_map_sr_groups

ggsave("som/titer_fold_changes/fc_D614G_map_sera.png", fc_map_sr_groups, width = 13, height = 8, dpi = 300)



