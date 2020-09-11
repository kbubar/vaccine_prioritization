setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________
library(tidyverse) 
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(ggrepel)

# cases <- readRDS("optimal_simplemodel_cases.RData")
# deaths <- readRDS("optimal_simple_initialtooldest2_deaths.RData")

cases_best <- matrix(NA, 51, 11)

#total cases: 50x
cases1 <- as.matrix(readRDS("optimal_BEL_cases1.RData")) # 3x
cases2 <- as.matrix(readRDS("optimal_BEL_cases2.RData")) # 8x
cases3 <- as.matrix(readRDS("optimal_BEL_cases3.RData")) # 8x
cases4 <- as.matrix(readRDS("optimal_BEL_cases4.RData")) # 8x
cases5 <- as.matrix(readRDS("optimal_BEL_cases5.RData")) # 8x
cases6 <- as.matrix(readRDS("optimal_BEL_cases6.RData")) # 8x
cases7 <- as.matrix(readRDS("optimal_BEL_cases7.RData")) # 7x
cases <- c(cases1, cases2, cases3, cases4, cases5, cases6, cases7)

# total mortalitly runs: 50x
deaths1 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_3x_1.RData")) # 3x
deaths2 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_3x_2.RData")) # 3x
deaths2 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_3x_3.RData")) # 3x
deaths3 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_8x_4.RData")) # 3x
deaths4 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_8x_5.RData")) # 3x
deaths5 <- as.matrix(readRDS("optimal_BEL_deaths_newIFR_8x_6.RData")) # 3x
deaths6 <- as.matrix(readRDS("optimal_BEL_deaths6.RData")) # 6x
deaths7 <- as.matrix(readRDS("optimal_BEL_deaths7.RData")) # 6x
deaths8 <- as.matrix(readRDS("optimal_BEL_deaths8.RData")) # 5x
deaths9 <- as.matrix(readRDS("optimal_BEL_deaths9.RData")) # 5x
deaths10 <- as.matrix(readRDS("optimal_BEL_deaths10.RData")) # 5x
deaths11 <- as.matrix(readRDS("optimal_BEL_deaths11.RData")) # 5x
cases <- c(deaths1, deaths2, deaths3, deaths4, deaths5) #, deaths6, deaths7,
           # deaths8, deaths9, deaths10, deaths11)

for (i in 0:51){
  index <- which.min(unlist(lapply(cases, function(x) min(x$Percent_infected[i]))))
  print(index)
  temp <- as.matrix(cases[[index]])
  cases_best[i,] <- temp[i,]
}

#cases_best <- list_optimal[[1]]
cases_best <- as.data.frame(cases_best)
colnames(cases_best) <- c("Percent_vax", "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+", "Percent_infected")

cases_best[cases_best > 100] <- 100
saveRDS(cases_best, "optimal_BEL_deaths_newIFR_x33.RData")
# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________

#cases_best <- readRDS("optimal_BEL_cases_50x.RData")
new_df <- cases_best %>%
  select("Percent_vax", "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+") %>%
  gather(key = "age_group", value = "num", -Percent_vax)

pal <- brewer.pal(9, "Greys")
# <- pal[2:6]

ggplot(new_df[new_df$Percent_vax < 43,], aes(x = Percent_vax, y = num, col = age_group)) +
  geom_line(size = 2, alpha = 0.9) +
  #geom_line(color = "black", size = 2.2) +
  xlab("Total vaccine supply (% of pop)") +
  ylab("% of age group vaccinated") +
  theme_classic(base_size = 26) +
  # geom_label(aes(x = new_df[50,]$Percent_vax + 2, y = new_df[50,]$num, 
  #                label = new_df[50,]$age_group, fill = new_df[50,]$age_group),
  #            size = 6, fontface = "bold") +
  # geom_label(aes(x = new_df[101,]$Percent_vax + 2, y = new_df[101,]$num,
  #                label = new_df[101,]$age_group, fill = new_df[101,]$age_group),
  #            size = 6, fontface = "bold") +
  # geom_label(aes(x = new_df[152,]$Percent_vax + 2, y = new_df[152,]$num,
  #                label = new_df[152,]$age_group, fill = new_df[152,]$age_group),
  #            size = 6, fontface = "bold") +
  # geom_label(aes(x = new_df[203,]$Percent_vax + 2, y = new_df[203,]$num - 1,
  #                label = new_df[203,]$age_group, fill = new_df[203,]$age_group),
  #            size = 6, fontface = "bold") +
  # geom_label(aes(x = new_df[254,]$Percent_vax + 2, y = new_df[254,]$num - 7,
  #                label = new_df[254,]$age_group, fill = new_df[254,]$age_group),
  #            size = 6, fontface = "bold", colour = "white") +
  # geom_label(aes(x = new_df[305,]$Percent_vax + 2, y = new_df[305,]$num,
  #                label = new_df[305,]$age_group, fill = new_df[305,]$age_group),
  #            size = 6, fontface = "bold", colour = "white") +
  # geom_label(aes(x = new_df[356,]$Percent_vax + 2, y = new_df[356,]$num,
  #                label = new_df[356,]$age_group, fill = new_df[356,]$age_group),
  #            size = 6, fontface = "bold", colour = "white") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  # ggtitle("Optimal allocation for minimizing cases") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 42), breaks = c(0, 10, 20, 30, 40))
