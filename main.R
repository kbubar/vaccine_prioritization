# Vaccine Strategy Simple Model

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________

library(tidyverse) 
library(deSolve) # ode solver
library(gridExtra)
library(RColorBrewer)
library(wesanderson)

# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)
my_pal <- c(col_kids, col_youngadults, col_elderly, col_adults, col_all)
# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# Demographics
country <- "BEL"

C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))

age_demo <- readRDS(paste0("age_demographics_", country,".RData"))

pop_total <- age_demo[10]
age_demo <- age_demo[1:9]

N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

#IFR   <- c(0.001, 0.001, 0.007, 0.02, 0.06, 0.2, 0.9, 2.4, 10.1) # Ref: Salje
IFR   <- c(0.001, 0.001, 0.005, 0.02, 0.05, 0.2, 0.7, 1.9, 8.3) # Ref: Salje
IFR   <- IFR/100 # as decimal

# susceptibility with R0 ~ 3 (R0 <- compute_R0(u, C))
u_constant     <- rep(0.022, 9) # constant # 0.02 for US, 0.022 for Belgium
#u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/32.8 # Ref: Davies
u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 # R0 = 2.6 for BEL
R0 <- compute_R0(u_var, C)

# vaccine efficacy
v_e_constant <- get_v_e(p = 1, y0 = 1, hinge_age = 50)
v_e_var <- get_v_e(p = 0.5, y0 = 1, hinge_age = 50)

# serology 
sero_none <- rep(0, 9) # no prior immunity
sero_belgium <- c(0.052, 0.052, 0.076, 0.055, 0.063, 0.071, 0.038, 0.042, 0.1) # Ref: Herzog

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________
# strategy options: "no vax", "all", "kids", "adults", "elderly"
# df_novax <- run_sim(C, 0, "no vax", u_constant, v_e_constant) DELETE LATER, MAY BE

# * Simple model (US: everything constant besides IFR) ----
list_all      <- vector(mode = "list")
list_kids     <- vector(mode = "list")
list_adults   <- vector(mode = "list")
list_elderly  <- vector(mode = "list")
list_twentyplus   <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(C, percent_vax = j, strategy = "all", u = u_var)
  list_kids[[paste0(i)]] <- run_sim(C, j, "kids", u_var)
  list_adults[[paste0(i)]] <- run_sim(C, j, "adults", u_var)
  list_elderly[[paste0(i)]] <- run_sim(C, j, "elderly", u_var)
  list_twentyplus[[paste0(i)]] <- run_sim(C, j, "20+", u_var)
}

# * Nontransmission blocking ----
list_all_ntb      <- vector(mode = "list")
list_kids_ntb      <- vector(mode = "list")
list_adults_ntb    <- vector(mode = "list")
list_elderly_ntb   <- vector(mode = "list")
list_twentyplus_ntb  <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_ntb[[paste0(i)]] <- run_sim_nontransmissionblocking(C, percent_vax = j, strategy = "all", u = u_var)
  list_kids_ntb[[paste0(i)]] <- run_sim_nontransmissionblocking(C, j, "kids", u_var)
  list_adults_ntb[[paste0(i)]] <- run_sim_nontransmissionblocking(C, j, "adults", u_var)
  list_elderly_ntb[[paste0(i)]] <- run_sim_nontransmissionblocking(C, j, "elderly", u_var)
  list_twentyplus_ntb[[paste0(i)]] <- run_sim_nontransmissionblocking(C, j, "20+", u_var)
  }

# * Varying susceptibility, u ----
list_all_u_var <- vector(mode = "list")
list_kids_u_var <- vector(mode = "list")
list_adults_u_var <- vector(mode = "list")
list_elderly_u_var <- vector(mode = "list")
list_twentyplus_u_var   <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){ 
  j <- i/100
  list_all_u_var[[paste0(i)]] <- run_sim(C, j, "all", u_var)
  list_kids_u_var[[paste0(i)]] <- run_sim(C, j, "kids", u_var)
  list_adults_u_var[[paste0(i)]] <- run_sim(C, j, "adults", u_var)
  list_elderly_u_var[[paste0(i)]] <- run_sim(C, j, "elderly", u_var)
  list_twentyplus_u_var[[paste0(i)]] <- run_sim(C, j, "20+", u_var)
}

# * Varying v_e ----
list_all_v_e_var <- vector(mode = "list")
list_kids_v_e_var <- vector(mode = "list")
list_adults_v_e_var <- vector(mode = "list")
list_elderly_v_e_var <- vector(mode = "list")
list_twentyplus_v_e_var <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_v_e_var[[paste0(i)]] <- run_sim(C, j, "all", u_var, v_e_var)
  list_kids_v_e_var[[paste0(i)]] <- run_sim(C, j, "kids", u_var, v_e_var)
  list_adults_v_e_var[[paste0(i)]] <- run_sim(C, j, "adults", u_var, v_e_var)
  list_elderly_v_e_var[[paste0(i)]] <- run_sim(C, j, "elderly", u_var, v_e_var)
  list_twentyplus_v_e_var[[paste0(i)]] <- run_sim(C, j, "20+", u_var, v_e_var)
}

# * Serology (Belgium data) with no serology testing ----
list_all_sero_notest      <- vector(mode = "list")
list_kids_sero_notest     <- vector(mode = "list")
list_adults_sero_notest   <- vector(mode = "list")
list_elderly_sero_notest  <- vector(mode = "list")
list_twentyplus_sero_notest   <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_sero_notest[[paste0(i)]] <- run_sim(C, percent_vax = j, strategy = "all", u_var, frac_age = age_demo ,
                                              N = N_i , sero =  sero_belgium )
  list_kids_sero_notest[[paste0(i)]] <- run_sim(C , j, "kids", u_var, frac_age = age_demo ,
                                                N = N_i , sero =  sero_belgium )
  list_adults_sero_notest[[paste0(i)]] <- run_sim(C , j, "adults", u_var, frac_age = age_demo ,
                                                  N = N_i , sero =  sero_belgium )
  list_elderly_sero_notest[[paste0(i)]] <- run_sim(C , j, "elderly", u_var, frac_age = age_demo ,
                                                   N = N_i , sero =  sero_belgium )
  list_twentyplus_sero_notest[[paste0(i)]] <- run_sim(C , j, "20+", u_var, frac_age = age_demo ,
                                                      N = N_i , sero =  sero_belgium )
}

# * Serology (Belgium data) WITH serology testing ----
list_all_sero_test      <- vector(mode = "list")
list_kids_sero_test     <- vector(mode = "list")
list_adults_sero_test   <- vector(mode = "list")
list_elderly_sero_test  <- vector(mode = "list")
list_twentyplus_sero_test   <- vector(mode = "list")

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all_sero_test[[paste0(i)]] <- run_sim(C , percent_vax = j, strategy = "all", u_var, frac_age = age_demo ,
                                             N = N_i , sero =  sero_belgium , sero_testing = TRUE)
  list_kids_sero_test[[paste0(i)]] <- run_sim(C , j, "kids", u_var, frac_age = age_demo ,
                                              N = N_i , sero =  sero_belgium , sero_testing = TRUE)
  list_adults_sero_test[[paste0(i)]] <- run_sim(C , j, "adults", u_var, frac_age = age_demo ,
                                                N = N_i , sero = sero_belgium , sero_testing = TRUE)
  list_elderly_sero_test[[paste0(i)]] <- run_sim(C , j, "elderly", u_var, frac_age = age_demo ,
                                                 N = N_i , sero =  sero_belgium , sero_testing = TRUE)
  list_twentyplus_sero_test[[paste0(i)]] <- run_sim(C , j, "20+", u_var, frac_age = age_demo ,
                                                    N = N_i , sero =  sero_belgium , sero_testing = TRUE)
}



# Run multiple R0 ----
scale_u <- c(49.5, 45, 43.1, 41.3, 39.6, 38.1, 36.7, 35.4, 34.2, 33, 31.9)
#scale_u <- c(49.5)

list_cases_R0      <- vector(mode = "list")
list_deaths_R0     <- vector(mode = "list")

count <- 1
ptm <- proc.time()
for (j in scale_u){
  u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/j
  
  list_all      <- vector(mode = "list")
  list_kids     <- vector(mode = "list")
  list_adults   <- vector(mode = "list")
  list_elderly  <- vector(mode = "list")
  list_twentyplus   <- vector(mode = "list")
  
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(C, percent_vax = j, strategy = "all", u = u_var)
    list_kids[[paste0(i)]] <- run_sim(C, j, "kids", u_var)
    list_adults[[paste0(i)]] <- run_sim(C, j, "adults", u_var)
    list_elderly[[paste0(i)]] <- run_sim(C, j, "elderly", u_var)
    list_twentyplus[[paste0(i)]] <- run_sim(C, j, "20+", u_var)
  }
  
  list_cases_R0[[count]] <- get_reduction_in_cases_df_novar()
  list_deaths_R0[[count]] <- get_reduction_in_deaths_df_novar()
  count <- count + 1
}
proc.time() - ptm
list_cases_R0 <- readRDS("list_BEL_overR0_cases.RData")
list_deaths_R0 <- readRDS("list_BEL_overR0_deaths.RData")

p1 <- plot_over_vax_avail_varyingR0(list_deaths_R0, "deaths", "2.1", 1) + labs(tag = "A ")
p2 <- plot_over_vax_avail_varyingR0(list_deaths_R0, "deaths", "2.6", 2) + labs(tag = "B")
p3 <- plot_over_vax_avail_varyingR0(list_deaths_R0, "deaths", "3.1", 3) + labs(tag = "C")

p4 <- plot_over_vax_avail_varyingR0(list_cases_R0, "cases", "2.1", 4) + labs(tag = "D ")
p5 <- plot_over_vax_avail_varyingR0(list_cases_R0, "cases", "2.6", 5) + labs(tag = "E")
p6 <- plot_over_vax_avail_varyingR0(list_cases_R0, "cases", "3.1", 6) + labs(tag = "F")

grid.arrange(arrangeGrob(p1, top = textGrob("R0 = 2.1", vjust = 1, gp = gpar(fontsize = 18)),
                         left = textGrob("Reduction in deaths (%)", rot = 90, hjust = 0.5,
                                         gp = gpar(fontsize = 20))),
             arrangeGrob(p2, top = textGrob("R0 = 2.6", vjust = 1, gp = gpar(fontsize = 18))),
             arrangeGrob(p3, top = textGrob("R0 = 3.1", vjust = 1, gp = gpar(fontsize = 18)),
                         right = ""),
             arrangeGrob(p4, left = textGrob("Reduction in infections (%)", rot = 90, hjust = 0.5, 
                                             gp = gpar(fontsize = 20))),
             arrangeGrob(p5),
             arrangeGrob(p6, right = ""),
             ncol=3, 
             widths=c(2.6, 2.3, 2.3),
             heights = c(2.4, 2.5),
             #left = textGrob("Reduction in infections (%)", rot = 90, vjust = 1, gp = gpar(fontsize = 18)), 
             bottom = textGrob("Total vaccine supply (% of pop)", vjust = 0.4, gp = gpar(fontsize = 20)))

# RESULTS ----
# _____________________________________________________________________

# * Paper fig 1 ----
# export 800*600
outcome <- "deaths"
pB <- plot_over_vax_avail(outcome, "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)
outcome <- "cases"
pC <- plot_over_vax_avail(outcome, "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)

# * * Plot vaccination strategies ----
p2 <- barplot_vax_strat("all")
p3 <- barplot_vax_strat("kids")
p4 <- barplot_vax_strat("adults")
p5 <- barplot_vax_strat("elderly")
p6 <- barplot_vax_strat("20+")

# export 1000 * 500
pA <- grid.arrange(arrangeGrob(p3,p4,p6,p5,p2,
                         ncol=1, 
                         widths=c(2.6), 
                         heights=c(2.3,2.3,2.3,2.3,3.1),
                         left = textGrob("Age-distribution of vaccines (%)", rot = 90, vjust = 0.5, gp = gpar(fontsize = 24)), 
                         bottom = textGrob("Age (years)", vjust = 0, gp = gpar(fontsize = 24))))


# * Paper fig 2
outcome <- "deaths"
plot_over_vax_avail(outcome, "Vaccine efficacy", list_all_v_e_var, list_kids_v_e_var, list_adults_v_e_var, list_elderly_v_e_var, list_twentyplus_v_e_var)

groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
groups <- c("0", "10", "20", "30", "40", "50", "60", "70", "80")
groups <- c(0,10,20,30,40,50,60,70,80)

df <- data.frame(groups, v_e_var, v_e_constant)
df$xend <- c(10,20,30,40,50,60,70,80,90)

theme_set(theme_minimal(base_size = 38))

#export as 800 * 900
ggplot(df, aes(x = groups, y = v_e_var*100, xend = xend, yend = v_e_var*100)) +
  geom_rect(aes(xmin=50, xmax=60, ymin=0, ymax=100), alpha = 0.02) +
  geom_segment(size = 1.5, linetype = "dashed")+
  #geom_point(aes(y = v_e_var*100, group = 1), size = 1.5) +
  geom_hline(yintercept = 100, size = 1.5) +
  ylab("Vaccine Efficacy (%)") +
  scale_y_continuous(expand = c(0,0), limit = c(0, 100.5)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 90), breaks = c(0,10,20,30,40,50,60,70,80,90),
                     labels = c(0,10,20,30,40,50,60,70,80,90)) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Age (years)") 

# * Paper Fig 3


# * Plot for all ages over time: infected & recovered ----
compartment <- "I"
N <- N_i

# NAS figs
p1 <- plot_allages_onestrategy(list_elderly_u_var[[1]], "I")
print(p1)

p2 <- plot_allages_onestrategy(list_elderly_u_var[[21]], "I")
print(p2)

# ggtitle("No Vaccines")

legend <- get_legend(p1)
#p1 <- p1 + theme(legend.position = "none")
p1 <- plot_allages_onestrategy(list_twentyplus_u_var[[20]], "I") +
  theme(legend.position = "none") 
p2 <- plot_allages_onestrategy(list_twentyplus_u_var[[25]], "I") +
  theme(legend.position = "none")
p3 <- plot_allages_onestrategy(list_twentyplus_u_var[[30]], "I") +
  theme(legend.position = "none")
p4 <- plot_allages_onestrategy(list_twentyplus_u_var[[35]], "I") +
  theme(legend.position = "none")
p5 <- plot_allages_onestrategy(list_twentyplus[[20]], "I")  +
  theme(legend.position = "none")
p6 <- plot_allages_onestrategy(list_twentyplus[[25]], "I")  +
  theme(legend.position = "none")
p7 <- plot_allages_onestrategy(list_twentyplus[[30]], "I")  +
  theme(legend.position = "none")
p8 <- plot_allages_onestrategy(list_twentyplus[[35]], "I")  +
  theme(legend.position = "none")

print(p1)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4, widths=c(2.3, 2.3, 2.3, 2.3))


p2 <- plot_allages_onestrategy(list_all_belgium[[2]], "I") +
  theme(legend.position = "none") +
  ggtitle("all, 1% Vax, constant u")

p3 <- plot_allages_onestrategy(list_all_belgium[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("all, 25% Vax, constant u")

p4 <- plot_allages_onestrategy(list_kids_belgium[[2]], "I") +
  theme(legend.position = "none") +
  ggtitle("Kids, 1% Vax, constant u")

p5 <- plot_allages_onestrategy(list_kids_belgium[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("Kids, 25% Vax, constant u")

p6 <- plot_allages_onestrategy(list_adults_belgium[[2]], "I") +
  theme(legend.position = "none") +
  ggtitle("Adults, 1% Vax, constant u")

p7 <- plot_allages_onestrategy(list_adults_belgium[[26]], "I") +
  theme(legend.position = "none") +
  ggtitle("Adults, 25% Vax, constant u")

#grid.arrange(p2, p3, p4, p5, p6, p7, ncol=2, widths=c(2.3, 2.3))
grid.arrange(p2, p3, p4, p5, p6, p7,legend, ncol=2, widths=c(2.3, 2.3))

# * Serology Comparisons ----
i <- 1
p1 <- plot_totalinfections_sero(list_all[[i]], list_all_sero_notest[[i]], list_all_sero_test[[i]], "I") +
  ggtitle("No vaccine")

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

i <- 6
p2 <- plot_totalinfections_sero(list_all[[i]], list_all_sero_notest[[i]], list_all_sero_test[[i]], "I") +
  ggtitle("Vaccines for 5% of pop") +
  theme(legend.position = "none")

i <- 11
p3 <- plot_totalinfections_sero(list_all[[i]], list_all_sero_notest[[i]], list_all_sero_test[[i]], "I") +
  ggtitle("Vaccines for 10% of pop") +
  theme(legend.position = "none")

i <- 26
p4 <- plot_totalinfections_sero(list_all[[i]], list_all_sero_notest[[i]], list_all_sero_test[[i]], "I") +
  ggtitle("Vaccines for 25% of pop") +
  theme(legend.position = "none")

grid.arrange(p1, p2, legend, p3, p4, ncol=3, widths=c(2.3, 2.3, 1))

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
t <- seq(0,400,1) 

# Simple Model
percent_vax <- 0.10
nvax <- percent_vax*pop_total

 y <- "infections"
#y <- "deaths"
p1 <- barplot_at_finalT(list_all$`0`, IFR, "no vax", y) 
p2 <- barplot_at_finalT(list_all$'25', IFR, "all", y, nvax = 0.25*pop_total) 
p3 <- barplot_at_finalT(list_kids$'25', IFR, "kids", y, nvax = 0.25*pop_total) 
p4 <- barplot_at_finalT(list_adults$'25', IFR, "adults", y, nvax = 0.25*pop_total) 
p5 <- barplot_at_finalT(list_elderly$'25', IFR, "elderly", y, nvax = 0.25*pop_total) 
p6 <- barplot_at_finalT(list_twentyplus$'25', IFR, "20+", y, nvax = 0.25*pop_total)

grid.arrange(p1, p2, p3, p4, p5, p6, ncol=1, widths=c(2.3))

p1 <- barplot_at_finalT(list_all$`0`, IFR, "no vax", y) 
p2 <- barplot_at_finalT(list_adults$'20', IFR, "adults", y, nvax = 0.20*pop_total) 
p3 <- barplot_at_finalT(list_adults$'30', IFR, "adults", y, nvax = 0.30*pop_total)
p4 <- barplot_at_finalT(list_adults$'40', IFR, "adults", y, nvax = 0.40*pop_total)
#grid.arrange(p1, p2, p3, ncol=1, widths=c(2.3))

p5 <- barplot_at_finalT(list_all_u_var$`0`, IFR, "no vax", y) 
p6 <- barplot_at_finalT(list_adults_u_var$'20', IFR, "adults", y, nvax = 0.20*pop_total) 
p7 <- barplot_at_finalT(list_adults_u_var$'30', IFR, "adults", y, nvax = 0.30*pop_total)
p8 <- barplot_at_finalT(list_adults_u_var$'40', IFR, "adults", y, nvax = 0.40*pop_total)
grid.arrange(p1, p5, p2, p6, p3, p7, p4, p8, ncol=2, widths=c(2.3, 2.3))


## 
# percent_vax <- 0.10
# nvax <- percent_vax*pop_total 
# 
# y <- "vaccinated"
# p2 <- barplot_at_finalT(list_all$'25', IFR, "all", y, nvax = percent_vax*pop_total) 
# p3 <- barplot_at_finalT(list_kids$'25', IFR, "kids", y, nvax = percent_vax*pop_total) 
# p4 <- barplot_at_finalT(list_adults$'25', IFR, "adults", y, nvax = percent_vax*pop_total) 
# p5 <- barplot_at_finalT(list_elderly$'25', IFR, "elderly", y, nvax = percent_vax*pop_total) 
# p6 <- barplot_at_finalT(list_twentyplus$'25', IFR, "20+", y, nvax = percent_vax*pop_total)




# * Bar plots ----
barplot_totalcases()
barplot_totaldeaths()


# * Plot total cases and deaths ----
outcome <- "cases"
outcome <- "deaths"

plot_over_vax_avail(outcome, "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)

plot_over_vax_avail(outcome, "Non-transmission blocking", list_all_ntb, list_kids_ntb, list_adults_ntb, list_elderly_ntb, list_twentyplus_ntb)

#plot_over_vax_avail(outcome, "Susceptibility", list_all_u_var, list_kids_u_var, list_adults_u_var, list_elderly_u_var, list_twentyplus_u_var)

plot_over_vax_avail(outcome, "Vaccine efficacy", list_all_v_e_var, list_kids_v_e_var, list_adults_v_e_var, list_elderly_v_e_var, list_twentyplus_v_e_var)

plot_over_vax_avail(outcome, "Contact Matrix", list_all_C_noschool, list_kids_C_noschool, list_adults_C_noschool, list_elderly_C_noschool)

plot_over_vax_avail(outcome, "Serology (Belgium)", list_all_sero_test, list_kids_sero_test, list_adults_sero_test, list_elderly_sero_test, list_twentyplus_sero_test)

outcome <- "cases"
outcome <- "deaths"
p0 <- plot_over_vax_avail(outcome, "Serology (Belgium)", list_all_sero_notest, list_kids_sero_notest, list_adults_sero_notest, list_elderly_sero_notest, list_twentyplus_sero_notest)
p1 <-  plot_sero_over_vax_avail(outcome, "kids", col_kids)
p2 <-  plot_sero_over_vax_avail(outcome, "young adults", col_youngadults)
p3 <-  plot_sero_over_vax_avail(outcome, "adults", col_adults)
p4 <-  plot_sero_over_vax_avail(outcome, "elderly", col_elderly)
p5 <-  plot_sero_over_vax_avail(outcome, "all", col_all)

grid.arrange(p0, p1, p2, p3, p4, p5, 
             ncol=6, widths=c(2.3, 2.3, 2.3, 2.3, 2.3, 2.3),
             left = textGrob("Reduction in Infections (%)", rot = 90, vjust = 1, gp = gpar(fontsize = 22)), 
             bottom = textGrob("Total vaccine supply (% of pop)", vjust = 0, gp = gpar(fontsize = 22)))

# * vaccine efficacy tipping point ----
v_e_bisection(1, 50)

# * Other Plots ----
# Plot susceptibility
groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
df <- data.frame(groups, u_var, u_constant)

theme_set(theme_bw(base_size = 15))

ggplot(df, aes(x = groups)) +
  geom_line(aes(y = u_constant, group = 1), size = 1.3) +
  geom_point(aes(y = u_constant, group = 1), size = 2) +
  geom_line(aes(y = u_var, group = 1), linetype = "dashed", size = 1.3) +
  geom_point(aes(y = u_var, group = 1), size = 2) +
  ylab("Susceptibility") +
  xlab("Age Groups")

# Plot vaccine efficacy
groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
groups <- c("0", "10", "20", "30", "40", "50", "60", "70", "80")
groups <- c(0,10,20,30,40,50,60,70,80)
df <- data.frame(groups, v_e_var, v_e_constant)

theme_set(theme_minimal(base_size = 22))

ggplot(df, aes(x = groups)) +
  #geom_line(aes(y = v_e_constant*100, group = 1), size = 1.3) +
  #geom_point(aes(y = v_e_constant*100, group = 1), size = 2) +
  #geom_line(aes(y = v_e_var*100, group = 1), linetype = "dashed", size = 1.3) +
  geom_step(y = v_e_var*100) +
  geom_point(aes(y = v_e_var*100, group = 1), size = 2) +
  #geom_vline(xintercept = 50, size = 1, linetype = "dashed") +
  ylab("Vaccine Efficacy (%)") +
  scale_y_continuous(expand = c(0,0), limit = c(0, 110)) +
  scale_x_continuous(expand = c(0,0), limit = c(-5, 81), breaks = c(0,10,20,30,40,50,60,70,80), 
                     labels=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")) +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank()) +
  xlab("Age (years)") 
  #ggtitle("Vaccine efficacy by age", subtitle = "")


lines(groups, v_e_var, type = "s")

# Plot serology
groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
groups <- c("0", "10", "20", "30", "40", "50", "60", "70", "80")
groups <- c(0,10,20,30,40,50,60,70,80)
df <- data.frame(groups, sero_belgium)

theme_set(theme_minimal(base_size = 22))

ggplot(df, aes(x = groups)) +
  geom_hline(yintercept = sum(sero_belgium*100*age_demo)) +
  #geom_line(aes(y = avg_sero, group = 1), size = 1, color = "black") +
  geom_point(aes(y = sero_belgium*100, group = 1), size = 2.5) +
  ylab("Seroprevalence (%)") +
  scale_y_continuous(expand = c(0,0), limit = c(0, 11)) +
  scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
                   labels=c("0-9    ", "10-19    ", "20-29    ", "30-39    ", "40-49    ", "50-59    ",
                            "60-69    ","70-79    ", "80+    ")) +
  xlab("") +
  theme(axis.text.x = element_text(angle=45),
        axis.ticks = element_line(size = 1),
        axis.line = element_line(size = 1, color = "black"))

# Plot IFR
# ggplot(data = NULL, aes(x = 1:9, y = IFR))+
#   theme_bw() +
#   geom_line(size = 1.3) + 
#   ylab("IFR (%)") +
#   geom_point(size = 2) + 
#   scale_x_discrete(name = "Age groups")

# Plot Age demographics
barplot(age_demo*100, xlab = ("Age Groups"), ylab = "Population (%)",
       names.arg = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
       ylim = c(0,25),
       cex.names  = 1.5, cex.axis = 1.5, cex.lab = 1.5)

# # Direct vs Indirect protection
# baseline_cases <- compute_total_cases(list_all$`0`)
# number_cases <- compute_total_cases(list_elderly$`15`)
# number_vac <- sum(list_elderly$'15'[1, 38:46])/pop_total * 100
# num_cases_prevented <- baseline_cases - number_cases
#   