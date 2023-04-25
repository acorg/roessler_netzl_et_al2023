
shrinkrange <- function(x, amount) {
  
  xn <- c(amount / 2, 1 - amount / 2)
  xn * diff(range(x)) + x[1]
  
}

scale_x_titer <- function(
  threshold = "<20",
  logthreshold = 0,
  axisname = "Titer"
) {
  
  scale_x_continuous(
    name = axisname,
    breaks = function(x) {
      logthreshold:ceiling(max(x))
    },
    labels = function(x) {
      output <- 2^x*10
      output[x == logthreshold] <- threshold
      output
    },
    minor_breaks = NULL
  )
  
}

scale_y_titer <- function(
  threshold = "<20",
  logthreshold = 0,
  axisname = "Titer",
  ymin = NULL,
  ymax = NULL,
  ...
  ) {
  
  scale_y_continuous(
    name = axisname,
    breaks = function(x) {
      if (is.null(ymin)) ymin <- ceiling(min(x))
      if (is.null(ymax)) ymax <- floor(max(x))
      ymin:ymax
    },
    labels = function(x) {
      output <- 2^x*10
      output[x == logthreshold] <- threshold
      output
    },
    minor_breaks = NULL,
    ...
  )
  
}

titerplot_theme <- function(){
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.background = element_rect(
      fill = NA,
      colour = "grey40"
    ),
    panel.grid.major = element_line(
      linetype = "solid",
      colour = rgb(0,0,0,0.05)
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold",
      size = 10
    )
  )
  
}

agFillScale <- function(map) {
  ag_fill <- agFill(map)
  names(ag_fill) <- agNames(map)
  ag_fill
}

srGroupOutline <- function(map) {
  sr_groups <- as.character(unique(srGroups(map)))
  sr_group_outlines <- srOutline(map)[match(sr_groups, srGroups(map))]
  names(sr_group_outlines) <- sr_groups
  sr_group_outlines
}
