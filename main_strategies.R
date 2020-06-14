# Vaccine Strategy Simple Model

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________

library(tidyverse) 
library(deSolve) # ode solver
library(xlsx)
library(gridExtra)
library(RColorBrewer)

# ColorBrewer Dark2
col_adults = "#1B9E77" # green
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple
col_kids = "#E7298A" # magenta


# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# Demographics
C <- readRDS("C_USA_bytens_all.RData")
C_noschool <- readRDS("C_USA_bytens_noschool.RData")

age_demo <- read.xlsx("USA_demographics.xlsx", "data_bytens", header = FALSE, row.names = TRUE)
age_demo <- as.vector(age_demo[[1]])/100

pop_total <- 100000
N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

IFR   <- c(0.001, 0.001, 0.007, 0.02, 0.06, 0.2, 0.9, 2.4, 10.1) # Ref: Salje
IFR   <- IFR/100 # as decimal

# susceptibility with R0 ~ 3 (R0 <- compute_R0(u))
u_constant     <- rep(0.02, 9) # constant
u_var     <- c(0.33, 0.37, 0.69, 0.81, 0.74, 0.8, 0.89, 0.77, 0.77)/30 # Ref: Davies

# vaccine efficacy
v_e_constant <- get_v_e(p = 1)
v_e_var <- get_v_e(p = 0.5)

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________
# strategy options: "no vax", "all", "kids", "adults", "elderly"
df_novax <- run_sim(C, 0, "no vax", u_constant, v_e_constant)
df_novax_u_var <- run_sim(C, 0, "no vax", u_var, v_e_constant)
df_novax_v_e_var <- run_sim(C, 0, "no vax", u_constant, v_e_var)
df_novax_C_noschool<- run_sim(C_noschool, 0, "no vax")

# * Simple model (everything constant besides IFR) ----
list_all      <- vector(mode = "list")
list_kids     <- vector(mode = "list")
list_adults   <- vector(mode = "list")
list_elderly  <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(C, percent_vax = j, strategy = "all", u = u_constant)
  list_kids[[paste0(i)]] <- run_sim(C, j, "kids", u_constant)
  list_adults[[paste0(i)]] <- run_sim(C, j, "adults", u_constant)
  list_elderly[[paste0(i)]] <- run_sim(C, j, "elderly", u_constant)
}

# * Varying susceptibility, u ----
list_all_u_var <- vector(mode = "list")
list_kids_u_var <- vector(mode = "list")
list_adults_u_var <- vector(mode = "list")
list_elderly_u_var <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){ 
  j <- i/100
  list_all_u_var[[paste0(i)]] <- run_sim(C, j, "all", u_var)
  list_kids_u_var[[paste0(i)]] <- run_sim(C, j, "kids", u_var)
  list_adults_u_var[[paste0(i)]] <- run_sim(C, j, "adults", u_var)
  list_elderly_u_var[[paste0(i)]] <- run_sim(C, j, "elderly", u_var)
}

# * Varying v_e ----
list_all_v_e_var <- vector(mode = "list")
list_kids_v_e_var <- vector(mode = "list")
list_adults_v_e_var <- vector(mode = "list")
list_elderly_v_e_var <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_v_e_var[[paste0(i)]] <- run_sim(C, j, "all", u_constant, v_e_var)
  list_kids_v_e_var[[paste0(i)]] <- run_sim(C, j, "kids", u_constant, v_e_var)
  list_adults_v_e_var[[paste0(i)]] <- run_sim(C, j, "adults", u_constant, v_e_var)
  list_elderly_v_e_var[[paste0(i)]] <- run_sim(C, j, "elderly", u_constant, v_e_var)
}

# * C no school ----
list_all_C_noschool <- vector(mode = "list")
list_kids_C_noschool <- vector(mode = "list")
list_adults_C_noschool <- vector(mode = "list")
list_elderly_C_noschool <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_C_noschool[[paste0(i)]] <- run_sim(C_noschool, j, "all")
  list_kids_C_noschool[[paste0(i)]] <- run_sim(C_noschool, j, "kids")
  list_adults_C_noschool[[paste0(i)]] <- run_sim(C_noschool, j, "adults")
  list_elderly_C_noschool[[paste0(i)]] <- run_sim(C_noschool, j, "elderly")
}

# _____________________________________________________________________



# RESULTS ----
# _____________________________________________________________________r

# * Plot for all ages over time: infected & recovered ----
# TODO: Update
compartment <- "I"

p1 <- plot_allages_onestrategy(list_adults_u_constant[[1]], "I") +
    ggtitle("No Vaccines")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- plot_allages_onestrategy(list_all_u_constant[[2]], "I") +
    theme(legend.position = "none") +
    ggtitle("all, 1% Vax, constant u")

p3 <- plot_allages_onestrategy(list_all_u_constant[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("all, 25% Vax, constant u")

p4 <- plot_allages_onestrategy(list_kids_u_constant[[2]], "I") +
  theme(legend.position = "none") +
  ggtitle("Kids, 1% Vax, constant u")

p5 <- plot_allages_onestrategy(list_kids_u_constant[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("Kids, 25% Vax, constant u")

p6 <- plot_allages_onestrategy(list_adults_u_constant[[2]], "I") +
  theme(legend.position = "none") +
  ggtitle("Adults, 1% Vax, constant u")

p7 <- plot_allages_onestrategy(list_adults_u_constant[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("Adults, 25% Vax, constant u")

#grid.arrange(p2, p3, p4, p5, p6, p7, ncol=2, widths=c(2.3, 2.3))
grid.arrange(p2, p3, p4, p5, p6, p7,legend, ncol=2, widths=c(2.3, 2.3))

# * Plot for one age group, different strats ----
# TODO: Update
p_kids_I <- plot_oneage_allstrategy(col_name = "I2", age_group_num = 2) +
  theme(legend.position = "none") +
  ggtitle("Ages 10-19")

p_young_I <- plot_oneage_allstrategy(col_name = "I4", age_group_num = 4) +
  theme(legend.position = "none") +
  ggtitle("Ages 30-39")

p_middle_I <- plot_oneage_allstrategy(col_name = "I6", age_group_num = 6) +
  theme(legend.position = "none") +
  ggtitle("Ages 50-59")

p_old_I <- plot_oneage_allstrategy(col_name = "I8", age_group_num = 8) +
  ggtitle("Ages 70+")

legend <- get_legend(p_old_I)
p_old_I <- p_old_I + theme(legend.position = "none")

grid.arrange(p_kids_I, p_young_I, p_middle_I, p_old_I, legend, ncol=5, widths=c(2.3, 2.3, 2.3, 2.3, 0.8))

# * Plot final time point ----
# Vax strategies from run_sim
percent_vax <- 0.25
nvax <- percent_vax*pop_total 
# propall
vax_propall <- nvax*age_demo
vax_propall[vax_propall > N_i] <- N_i[vax_propall > N_i]

# propkids strategy
nkids <- N_i[1] + N_i[2]
vax_dist_propkids <- rep(0, num_groups)
vax_dist_propkids[1] <- N_i[1]/nkids
vax_dist_propkids[2] <- N_i[2]/nkids

vax_propkids <- nvax*vax_dist_propkids
vax_propkids[vax_propkids > N_i] <- N_i[vax_propkids > N_i]

# propadults strategy
nadults <- N_i[3] + N_i[4] + N_i[5]
vax_dist_propadults <- rep(0, num_groups)
vax_dist_propadults[3] <- N_i[3]/nadults
vax_dist_propadults[4] <- N_i[4]/nadults
vax_dist_propadults[5] <- N_i[5]/nadults

vax_propadults <- nvax*vax_dist_propadults
vax_propadults[vax_propadults > N_i] <- N_i[vax_propadults > N_i]


# propelderly strategy
nelderly <- N_i[7] + N_i[8] + N_i[9]
vax_dist_propelderly <- rep(0, num_groups)
vax_dist_propelderly[7] <- N_i[7]/nelderly
vax_dist_propelderly[8] <- N_i[8]/nelderly
vax_dist_propelderly[9] <- N_i[9]/nelderly

vax_propelderly <- nvax*vax_dist_propelderly
vax_propelderly[vax_propelderly > N_i] <- N_i[vax_propelderly > N_i]

t <- seq(0,150,1) 

# Simple Model

# y <- "infections"
y <- "deaths"
p1 <- barplot_at_finalT(df_novax, IFR, "no vax", y) 
p2 <- barplot_at_finalT(list_all$'25', IFR, "all", y) 
p3 <- barplot_at_finalT(list_kids$'25', IFR, "kids", y) 
p4 <- barplot_at_finalT(list_adults$'25', IFR, "adults", y) 
p5 <- barplot_at_finalT(list_elderly$'25', IFR, "elderly", y) 

grid.arrange(p1, p2, p3, p4, p5, ncol=1, widths=c(2.3))

## Changing susceptibility
y <- "infections"
# y <- "deaths"
p1 <- barplot_at_finalT(df_novax_ui, IFR, "no vax", y) 
p2 <- barplot_at_finalT(list_all_u_var$'25', IFR, "all", y) 
p3 <- barplot_at_finalT(list_kids_u_var$'25', IFR, "kids", y) 
p4 <- barplot_at_finalT(list_adults_u_var$'25', IFR, "adults", y) 
p5 <- barplot_at_finalT(list_elderly_u_var$'25', IFR, "elderly", y) 

grid.arrange(p1, p2, p3, p4, p5, ncol=1, widths=c(2.3))

## Changing contact matrix
#y <- "infections"
 y <- "deaths"
p1 <- barplot_at_finalT(df_novax_C_noschool, IFR, "no vax", y) 
p2 <- barplot_at_finalT(list_all_C_noschool$'25', IFR, "all", y) 
p3 <- barplot_at_finalT(list_kids_C_noschool$'25', IFR, "kids", y) 
p4 <- barplot_at_finalT(list_adults_C_noschool$'25', IFR, "adults", y) 
p5 <- barplot_at_finalT(list_elderly_C_noschool$'25', IFR, "elderly", y) 

grid.arrange(p1, p2, p3, p4, p5, ncol=1, widths=c(2.3))

## 
y <- "vaccinated"
p2 <- barplot_at_finalT(list_all_u_var$'25', IFR, "all", y) 
p3 <- barplot_at_finalT(list_kids_u_var$'25', IFR, "kids", y) 
p4 <- barplot_at_finalT(list_adults_u_var$'25', IFR, "adults", y) 
p5 <- barplot_at_finalT(list_elderly_u_var$'25', IFR, "elderly", y) 

grid.arrange(p2, p3, p4, p5, ncol=1, widths=c(2.3))
# * Bar plots ----
barplot_totalcases()
barplot_totaldeaths()


# * Plot total cases and deaths ----
# TODO: put in helper_functions.R
outcome <- "deaths"

plot_over_vax_avail(outcome, "Susceptibility", df_novax_u_var, list_all_u_var, list_kids_u_var, list_adults_u_var, list_elderly_u_var)

plot_over_vax_avail(outcome, "Vaccine efficacy", df_novax_v_e_var, list_all_v_e_var, list_kids_v_e_var, list_adults_v_e_var, list_elderly_v_e_var)

red <- plot_over_vax_avail(outcome, "Contact Matrix", df_novax_C_noschool, list_all_C_noschool, list_kids_C_noschool, list_adults_C_noschool, list_elderly_C_noschool)



# * Other Plots ----
# # Plot susceptibility
# groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
# df <- data.frame(groups, u_var, u_constant)
# 
# theme_set(theme_bw(base_size = 15))
# 
# ggplot(df, aes(x = groups)) +
#   geom_line(aes(y = u_constant, group = 1), size = 1.3) +
#   geom_point(aes(y = u_constant, group = 1), size = 2) +
#   geom_line(aes(y = u_var, group = 1), linetype = "dashed", size = 1.3) +
#   geom_point(aes(y = u_var, group = 1), size = 2) +
#   ylab("Susceptibility") + 
#   xlab("Age Groups")

# Plot vaccine efficacy
# groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
# df <- data.frame(groups, v_e_var, v_e_constant)
# 
# theme_set(theme_bw(base_size = 15))
# 
# ggplot(df, aes(x = groups)) +
#   geom_line(aes(y = v_e_constant*100, group = 1), size = 1.3) +
#   geom_point(aes(y = v_e_constant*100, group = 1), size = 2) +
#   geom_line(aes(y = v_e_var*100, group = 1), linetype = "dashed", size = 1.3) +
#   geom_point(aes(y = v_e_var*100, group = 1), size = 2) +
#   ylab("Vaccine Efficacy (%)") + 
#   xlab("Age Groups") + 
#   ylim(0, 100)

# Plot IFR
# ggplot(data = NULL, aes(x = 1:9, y = IFR))+
#   theme_bw() +
#   geom_line(size = 1.3) + 
#   ylab("IFR (%)") +
#   geom_point(size = 2) + 
#   scale_x_discrete(name = "Age groups")

# barplot(age_demo, xlab = ("Age Groups"), ylab = "Percent",
#        names.arg = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))

