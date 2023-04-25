# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ggplot2)

results <- readRDS("som/cross_validation/cross_validation_titer.rds")

# Read the map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-positioned.ace")

agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}


# Set detectable results subset
detectable_results <- filter(results, measured_titer_type == 1)

# do histogram of pred - measured
mean <- round(mean(detectable_results$residual, na.rm = T),2)
sd <- round(sd(scale(detectable_results$residual, scale = F), na.rm = T), 2)

hist_diff <- ggplot(detectable_results) +
  geom_histogram(aes(x = residual), fill = "grey50", alpha = 0.8, bins = 100) +
  xlim(c(-15, 15)) +
  geom_vline(xintercept = mean, linetype = "dashed") +
  labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
  theme_bw()

ggsave(plot = hist_diff, filename = "./som/cross_validation/histogram_residuals.png", width = 5, height = 4, dpi = 300)

# change names
ag_pretty <- data.frame(
  row.names = c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", 
                'CB.1', 'BR.3', 'CH.1.1', 'BA.5.3.2', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF'),
  val = c('D614G', 'alpha', 'alpha+E484K', 'gamma', 'beta', 'delta', 'BA.1 omicron', 'BA.2 omicron',
          'CB.1', 'BR.3', 'CH.1.1', 'BA.5 omicron', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF')
)

sr_pretty <- data.frame(
  row.names = c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', "CK.2.1.1 conv."),
  val = c('2xmRNA-1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"Anc. virus conv.", 'alpha conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', "CK.2.1.1 conv.")
)

detectable_results$ag_pretty <- factor(ag_pretty[as.character(detectable_results$ag_name),], levels = ag_pretty$val)
detectable_results$sr_pretty <- factor(sr_pretty[as.character(detectable_results$sr_group),], levels = sr_pretty$val)

    # Antigen and serum group tab
        detectable_results %>%
          ggplot(
            aes(
              x = predicted_logtiter,
              y = measured_logtiter,
              color = ag_name
            )
          ) +
          labs(x = "Predicted log2 titer",
               y = "Measured log2 titer") +
          # geom_smooth() +
          geom_point(
            alpha = 0.1
          ) +
          geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
          ) +
          scale_color_manual(
            # values = Racmacs:::srGroupOutline(map)
            values = agFillScale(map)
          ) +
          xlim(c(-10,10))+
          ylim(c(0, 8))+
          facet_grid(
            rows = vars(sr_pretty),
            cols = vars(ag_pretty)
          ) +
          theme_bw() +
          theme(legend.position = "none",
               strip.text.x = element_text(size = 6),
               strip.text.y = element_text(size = 6))-> gp

 
ggsave(plot = gp, filename = "./som/cross_validation/scatter_pred_vs_measured.png", width = 14, height = 10, dpi = 300)
 