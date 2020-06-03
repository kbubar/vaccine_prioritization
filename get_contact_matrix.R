# Combine 5 year age-group contact matrix from Prem for 10 years age-groups 

library("xlsx")
library("tidyverse")

# setwd to file with Prem contact matrices
setwd("~/Vaccine Strategy/Prem_contact_matrices_152_countries")

# choose country of interest and read in the file
country <- "United States of America"

C_byfives <- read.xlsx("MUestimates_all_locations_2.xlsx", country, header = FALSE)

# convert C to 10 year age-groups
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

# save as .RData  
# saveRDS(C_bytens, "C_USA_bytens.RData")

heatmap(C_bytens, NA, NA, xlab = "Age of Individual", ylab = "Age of Contact")
   