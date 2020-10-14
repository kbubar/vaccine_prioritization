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
  l <- dim(C_byfives)[1]
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

# country codes:
# Belgium: BEL
# United States: USA
# India: IND
# Spain: ESP
# Zimbabwe: ZWE
# Brazil: BRA
# China: CHN
# South Africa: ZAF
# Poland: POL

country <- "POL"
setting <- "overall" #overall, rural or urban
#* C for all locations (home, work, school & other) ----
df <- read.csv("synthetic_contacts_2020.csv")

df <- df[df$iso3c == country,]
df <- df[df$location_contact == "all",]
df <- df[df$setting == setting,]

nage <- 16 # number of age groups

# construct C(i,j) = # of age j people that an age i person contacts each day
C <- matrix(ncol = nage, nrow = nage)

contactor <- c("0 to 4" , "5 to 9", "10 to 14" ,"15 to 19", "20 to 24", "25 to 29", 
               "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59",
               "60 to 64", "65 to 69", "70 to 74", "75+")
contactee <- contactor

rownames(C) <- as.list(contactor)
colnames(C) <- as.list(contactee)

for (i in contactor){
  for (j in contactee){
    temp <- df[df$age_contactor == i & df$age_cotactee == j,]
    C[i,j] = temp$mean_number_of_contacts
  }
}

C_bytens <- convert_bins_5to10(C)
heatmap(C_bytens, NA, NA, scale = "column", xlab = "Age of Individual", ylab = "Age of Contact")

C_bytens <- add_80bin(C_bytens)
heatmap(C_bytens, NA, NA, scale = "none", 
        #xlab = "Age of Individual", 
        #ylab = "Age of Contact", 
        cexRow = 1.5, 
        cexCol = 1.5)

saveRDS(C_bytens, paste0("C_", country,"_bytens_", setting, ".RData"))

