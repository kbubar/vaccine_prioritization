# Combine 5 year age-group contact matrix from Prem for 10 years age-groups 

# setwd to file with Prem contact matrices
setwd("~/Vaccine Strategy/Vaccine_Allocation_Project/Prem_contact_matrices")

# IMPORT ----
library("xlsx")
library("tidyverse")
library("RColorBrewer")
library(gplots)

# FUNCTION ---- 
# convert C to 10 year age-groups
convert_bins_5to10 <- function(C_byfives) {
  l <- length(C_byfives)
  C_bytens <- matrix(nrow = l/2, ncol = l/2)
  
  col_count <- 1
  for (col in seq(1,l,by = 2)){
    row_count <- 1
    for (row in seq(1,l, by = 2)){
      C_bytens[row_count, col_count] <- C_byfives[row, col] +
        C_byfives[row, col + 1] +
        C_byfives[row + 1, col] +
        C_byfives[row + 1, col + 1]
      row_count <- row_count + 1
    }
    col_count <- col_count + 1
  }
  
  colnames(C_bytens) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
  rownames(C_bytens) <- colnames(C_bytens)
  
  C_bytens
}

add_80bin <- function(C_bytens){
  # INPUT:
  #  C_bytens: Contact matrix with 10 year age bins, with last bin = 70-79
  #
  # OUTPUT:
  # C_new: Same C matrix with new row & col for 80+ age bin
  
  bin_70 <- dim(C_bytens)[1]
  bin_80 <- bin_70 + 1
  bin_60 <- bin_70 - 1
  
  C_new <- matrix(nrow = bin_80, ncol = bin_80)
  colnames(C_new) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  rownames(C_new) <- colnames(C_new)
  
  # fill rows for 80+ with same data from 70-79
  C_new[1:bin_70, 1:bin_70] <- C_bytens
  C_new[2:bin_80, bin_80] <- C_bytens[1:bin_70, bin_70]
  C_new[bin_80, 2:bin_80] <- C_bytens[bin_70, 1:bin_70]
  
  C_new[1, bin_80] <- C_new[1, bin_70]
  C_new[bin_80, 1] <- C_new[bin_70, 1]
  
  # Assumption: 80+ have similar (but less) contact rates as 70-79, with increased contact with 70-80+ (long term living facilities)
  # Implementation: Decrease contact for bins 0 - 69 for 80+  by 10% then split these contacts between 70-79 & 80+
  
  # col 80+ 
  to_decrease <- 0.1 * C_new[1:bin_60, bin_80]
  C_new[1:bin_60, bin_80] <- 0.9 * C_new[1:bin_60, bin_80]
  
  to_increase <- sum(to_decrease)/2
  C_new[bin_70, bin_80] <- C_new[bin_70, bin_80] + to_increase
  C_new[bin_80, bin_80] <- C_new[bin_80, bin_80] + to_increase
  
  # row 80+ 
  to_decrease <- 0.1 * C_new[bin_80, 1:bin_60]
  C_new[bin_80, 1:bin_60] <- 0.9 * C_new[bin_80, 1:bin_60]
  
  to_increase <- sum(to_decrease)/2
  C_new[bin_80, bin_70] <- C_new[bin_80, bin_70] + to_increase
  C_new[bin_80, bin_80] <- C_new[bin_80, bin_80] + to_increase
  
  C_new
}

# MAKE CONTACT MATRICES ----
# Read in contact matrix 
# Convert 5 years age bins to 10 year age bins
# Add 80+ age bin

#country <- "United States of America"
country <- "Belgium"
#* C for all locations (home, work, school & other) ----
C_byfives_all <- read.xlsx("MUestimates_all_locations_2.xlsx", country, header = FALSE)
C_bytens_all <- convert_bins_5to10(C_byfives_all)
heatmap(C_bytens_all, NA, NA, scale = "column", xlab = "Age of Individual", ylab = "Age of Contact")

C_bytens_all <- add_80bin(C_bytens_all)
heatmap(C_bytens_all-belgium, NA, NA, scale = "none", 
        #xlab = "Age of Individual", 
        #ylab = "Age of Contact", 
        cexRow = 1.5, 
        cexCol = 1.5)



#* C without school (home, work & other) ----
school <- read.xlsx("MUestimates_school_2.xlsx", country, header = FALSE)
C_byfives_noschool <- C_byfives_all - school
C_bytens_noschool <- convert_bins_5to10(C_byfives_noschool)
heatmap(C_bytens_noschool, NA, NA, scale = "none", xlab = "Age of Individual", ylab = "Age of Contact")

C_bytens_noschool <- add_80bin(C_bytens_noschool)
heatmap(C_bytens_noschool, NA, NA, scale = "column", xlab = "Age of Individual", ylab = "Age of Contact")

# SAVE FILE as .RData  ----
 saveRDS(C_bytens_all, "C_Belgium_bytens_all.RData")
# saveRDS(C_bytens_noschool, "C_USA_bytens_noschool.RData")

dif <- matrix(nrow = 9, ncol = 9)
# check symmetry with age demographics
for (i in 1:9){
  for (j in 1:9){
    dif[i, j] <- C_bytens_all[j,i]*age_demo[j] - C_bytens_all[i,j]*age_demo[i]
  }
}
max(abs(dif))
