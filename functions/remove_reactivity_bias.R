# remove reactivity bias from titertable
remove_reactivity_bias_logtiter <- function(table) {
  row_names <- rownames(table)
  table <- sapply(as.data.frame(table), function(x) {
    x - (mean(x, na.rm = T) - mean(table, na.rm = T))
  })
  
  rownames(table) <- row_names
  
  return(table)
}