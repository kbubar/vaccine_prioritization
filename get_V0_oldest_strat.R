# Get V_0 for strategy to vaccinate oldest first (oldest_strat)
# For a given population demographic and vax supply 0:50

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# country codes:
# Belgium: BEL
# United States: USA
# India: IND
# Spain: ESP
# Zimbabwe: ZWE
# Brazil: BRA
# China: CHN
# South Africa: ZAF

country <- "POL"
age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total <- age_demo[10]
age_demo <- age_demo[1:9]
N_i <-  pop_total*age_demo      
num_groups <- length(age_demo)


V_0 <- matrix(0, 51, 9)

for (i in seq(0, 50, by = 1)){
  nvax <- pop_total*i/100
  j <- num_groups # start at the oldest age group
  
  while (nvax > N_i[j]){
    V_0[i+1,j] <- N_i[j]
    nvax <- nvax - N_i[j]
    j <- j - 1
  }
  V_0[i+1,j] <- nvax
}

saveRDS(V_0, paste0("V0_oldest_strat_", country,".RData"))
        