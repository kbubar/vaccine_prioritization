# Combine 5 year age-group demographics from UN into 10 years age-groups 
setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

library(grid)
library(readxl)

# population in 1000s
df <- read_excel("WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx", "ESTIMATES")

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

countrycode <- "POL"
country <- "Poland"
popdata <- df[df$Country == country,]
popdata <- popdata[popdata$Year == 2020,]

pop_demo <- rep(0, 10)
col_num <- 9 # col of 0-4 y.o.

for (i in 1:8){
  pop_demo[i] <- as.numeric(popdata[col_num]) + as.numeric(popdata[col_num + 1])
  col_num <- col_num + 2
}

pop_demo[9] <- sum(as.numeric(popdata[col_num:29]))
total_pop <- sum(pop_demo)
pop_demo <- pop_demo/total_pop
pop_demo[10] <- total_pop*1000

saveRDS(pop_demo, paste0("age_demographics_", countrycode, ".RData"))
