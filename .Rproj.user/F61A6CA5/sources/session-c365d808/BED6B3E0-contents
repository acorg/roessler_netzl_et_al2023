# format titer table from excel sheet
rm(list = ls())
library(tidyverse)
library(stringr)

# vaccine status for sheet I
vacc_stat <- c("n" = "",
               "1" = "ChAdOx1-S/ChAdOx1-S", 
               "2" = "ChAdOx1-S/BNT162b2",
               "3" = "BNT162b2/BNT162b2",
               "4" = "mRNA-1273/mRNA-1273",
               "d" = "divers")

# function for sr_info

sr_info_function <- function(x, kimpel1) {
  stat <- kimpel1[[x, "Immune status"]]
  vacc <- kimpel1[[x,4]]
  var <- kimpel1[[x, "Variant"]]

  if(var == "firstwave"){
    var <- "WT"
  }

  if(stat == "convalescent" | stat == "unvaccinated") {
    res <- paste0(var, " conv.")
    
  } else if(stat == "3x vaccinated"){
    res <- vacc
  } else if(str_detect(stat, "con/vacc")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(var, "+", substring(vacc, 4, nchar(vacc)-1))
    } else {
      res <- paste0(var, "+", vacc_stat[vacc])
    }
    
  } else if(str_detect(stat, "vacc/con")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(substring(vacc, 4, nchar(vacc)-1), "+", var)
    } else {
      res <- paste0(vacc_stat[vacc],"+", var)
    }
    
  } else {
    res <- vacc_stat[vacc]
  }
  
  return(res)
}

shorten_vacc_names <- function(sr_group) {
  sr_group <- gsub("ChAdOx1-S", "AZ", sr_group)
  sr_group <- gsub("ChAdOx1", "AZ", sr_group)
  sr_group <- gsub("BNT162b2", "BNT", sr_group)
  sr_group <- gsub("BNT162b", "BNT", sr_group)
  sr_group <- gsub("Ad26.COV2.S", "J\\&J", sr_group)
  sr_group <- gsub("mRNA\\-1273", "mRNA1273", sr_group)
  
  sr_group <- gsub("B.1.617.2", "delta", sr_group)
  sr_group <- gsub("B.1.351", "beta", sr_group)
  sr_group <- gsub("B.1.1.7\\/B.1.1.7\\+E484K", "alpha\\/alpha\\+E484K", sr_group)
  
  return(sr_group)
}

set_threshold <- function(tab, thresh) {
  tab[as.numeric(tab) < as.numeric(thresh)] <- paste0("<", thresh)
  tab[is.na(tab)] <- "*"
  
  return(tab)
}


# =============================== Omicron I sheet, published in NEJM
kimpel1 <- readxl::read_excel("./data/titer_data/230411_Omicron V_update.xlsx", sheet = 1)
colnames(kimpel1) <- kimpel1[2,]
kimpel1 <- kimpel1[c(3:120),]

# remove all NA rows
kimpel1 <- kimpel1[rowSums(is.na(kimpel1)) != ncol(kimpel1),]

kimpel1$Variant <- gsub(" ", "", kimpel1$Variant)
kimpel1$Variant <- gsub("\\(BA.5\\)", "", kimpel1$Variant)


all_ags <- c("D614G", "B.1.1.7", "B.1.1.7+E484K", "B.1.351", "P.1.1", "B.1.617.2", "BA.1",
  "BA.2", "CB.1", "BR.3","CH.1.1","BA.5.3.2", "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18", "XBB.1", "XBB.1.5", "XBF")

# add sample info
# save sample ID, then sr group, then sr info including time since dose
kimpel1 <- kimpel1 %>%
  mutate(
    sr_info = unlist(lapply(1:nrow(kimpel1), function(x) sr_info_function(x, kimpel1)
  )),
  sr_group = shorten_vacc_names(sr_info)
  )

# set serum_info_column that contains all info for map
kimpel1 <- kimpel1 %>%
  mutate(sr_info_complete = paste(`sample ID`,sr_group, sr_info, sep = "_"))

# pivot wider to titer table
kimpel1 %>%
  column_to_rownames("sr_info_complete") %>%
  select(all_of(all_ags[all_ags %in% colnames(kimpel1)])) -> kimpel1_table

#======================================  Hybrid Immunity
kimpel <- readxl::read_excel("./data/titer_data/230411_Omicron V_update.xlsx", sheet = 2)
View(kimpel)
# add sample info
# save sample ID, then sr group, then sr info including time since dose
kimpel <- kimpel %>%
  mutate(
    sr_group = `Immune status`, 
  # account for infection sequence
  sr_info_complete = paste(`sample ID`,sr_group, sep = "_")
  )

# select titer columns
# pivot wider to titer table
kimpel %>%
  column_to_rownames("sr_info_complete") %>%
  select(all_of(all_ags[all_ags %in% colnames(kimpel)])) -> kimpel2_table

#======================================  Combine tables


kimpel_table <- plyr::rbind.fill(kimpel1_table, kimpel2_table)
rownames(kimpel_table) <- c(rownames(kimpel1_table), rownames(kimpel2_table))
kimpel_table_wide <- t(kimpel_table)
# set not titrated to "*"
kimpel_table_wide[kimpel_table_wide == "pending"] <- "*"
kimpel_table_wide[kimpel_table_wide == "-"] <- "*"
kimpel_table_wide[is.na(kimpel_table_wide)] <- "*"

# add lower threshold at the end
kimpel_table_wide <- set_threshold(kimpel_table_wide, 16)

# table with all samples
write.csv(kimpel_table_wide, "./data/titer_data/titer_table_all_samples.csv")

# remove samples from titer table 
kimpel_table_wide <- kimpel_table_wide[, !(grepl("C711|C860|G776|G780", colnames(kimpel_table_wide)))]

# remove double samples 
kimpel_table_wide <- kimpel_table_wide[, !(grepl("G935|G803|G898", colnames(kimpel_table_wide)))]
write.csv(kimpel_table_wide, "./data/titer_data/titer_table.csv")

