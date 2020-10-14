# Get YLL vec for a specific country from WHO GHO data
# https://apps.who.int/gho/data/view.main.LT62160?lang=en

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

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# read in YLL data
df <- read_excel("WPP2019_MORT_F16_1_LIFE_EXPECTANCY_BY_AGE_BOTH_SEXES.xlsx", "ESTIMATES", skip = 16)
demo_df <- read_excel("WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx", "ESTIMATES")

country <- "Belgium"
countrycode <- "BEL"
data <- df[df$Type == "Country/Area",]
data <- data[data$`Region, subregion, country or area *` == country,]
data <- data[data$Period == "2015-2020",]                 

# read in pop demographics to weigh YLL

popdata <- demo_df[demo_df$Country == country,]
popdata <- popdata[popdata$Year == 2020,]

pop_demo <- rep(0, 20)
yll_data <- rep(0, 20)
col_num_pop <- 9 # col of 0-4 y.o.
col_num_yll <- 10

for (i in 1:21){
  pop_demo[i] <- as.numeric(popdata[col_num_pop])
  yll_data[i] <- as.numeric(data[col_num_yll])
  col_num_pop <- col_num_pop + 1
  col_num_yll <- col_num_yll + 1
}

# weighted avg yll binned into 10y age bins (0-9, 10-19, ..., 80+)
yll_vec <- rep(0, 9)

count <- 1
for (i in 1:8){
  temp_sum <- pop_demo[count] + pop_demo[count + 1]
  yll_vec[i] <- yll_data[count]*pop_demo[count]/temp_sum + yll_data[count+1]*pop_demo[count+1]/temp_sum
  count <- count + 2
}

temp_sum <- sum(pop_demo[count:21]) 
yll_vec[9] <- sum(yll_data[count:21]*pop_demo[count:21]/temp_sum) 

saveRDS(yll_vec, paste0("yll_vec_", countrycode, ".RData"))
