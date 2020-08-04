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

# cases1 <- as.matrix(readRDS("optimal_belgium_new_cases.RData"))
# cases2 <- as.matrix(readRDS("optimal_belgium_new_cases_2.RData"))

cases_best <- matrix(NA, 51, 11)

cases1 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases1.RData"))
cases2 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases2.RData")) # noise +/ 3%
cases3 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases3.RData"))
cases4 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases4.RData"))
cases5 <- as.matrix(readRDS("optimal_belgium_uvar_list_cases1.RData"))
cases6 <- as.matrix(readRDS("optimal_belgium_uvar_list_cases2.RData"))
cases7 <- as.matrix(readRDS("optimal_belgium_uvar_list_cases3.RData"))
cases8 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases5.RData"))
cases9 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases6.RData"))
cases10 <- as.matrix(readRDS("optimal_belgium_40_50_list_cases7.RData"))

# cases <- c(cases1, cases2, cases3, cases4, cases5, cases6)
cases <- c(cases1, cases2, cases3, cases5, cases6, cases7, cases8, cases9, cases10)
# deaths1 <- readRDS("optimal_belgium_list_deaths.RData")
# deaths2 <- readRDS("optimal_belgium_list_deaths2.RData")
# deaths3 <- readRDS("optimal_belgium_list_deaths3.RData")
# cases <- c(deaths1, deaths2, deaths3)

for (i in 0:51){
  index <- which.min(unlist(lapply(cases, function(x) min(x$Percent_infected[i]))))
  print(index)
  temp <- as.matrix(cases[[index]])
  cases_best[i,] <- temp[i,]
}

#cases_best <- list_optimal[[1]]
cases_best <- as.data.frame(cases_best)
colnames(cases_best) <- c("Percent_vax", "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+", "Percent_infected")

saveRDS(cases_best, "optimal_belgium_u_var_best.RData")
# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________

new_df <- cases_best %>%
  select(Percent_vax, "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79") %>% #, "80+") %>%
  gather(key = "age_group", value = "num", -Percent_vax)

pal <- brewer.pal(9, "Greys")
pal <- pal[2:9]

ggplot(new_df, aes(x = Percent_vax, y = num, col = age_group)) +
  geom_line(size = 2, alpha = 0.9) +
  #geom_line(color = "black", size = 2.2) +
  xlab("Vaccine available (% of total pop)") +
  ylab("% of age group vaccinated") +
  theme_classic(base_size = 20) +
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
  #ylim(0, 101) + 
  # xlim(0, 52) + 
  scale_y_continuous(expand = c(0,0), limit = c(0, 101)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 55), breaks = c(0, 10, 20, 30, 40, 50))
