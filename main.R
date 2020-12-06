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
library(egg)
library(foreach)
library(doParallel)
library(ggplot2)

# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green (20-49)
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple (60+)
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)

# myPal<- c(col_adults, col_all, col_youngadults,
#               col_elderly, col_kids)
# names(myPal) <- c("twentyplus", "all", "adults",
#                      "elderly", "kids")
# tippingpointPal <- c(col_adults, col_all, col_youngadults,
#                      "#E6E6E6", col_kids)
# names(tippingpointPal) <- c("twentyplus", "all", "adults",
#                   "elderly", "kids")
# colFill <- scale_fill_manual(name = "strats", values = tippingpointPal)

tippingpointPal2 <- c("#E6E6E6", "#FFEDA0", "#FEB24C", "#F03B20", "#000000")
names(tippingpointPal2) <- c("None", "0-25%", "25-50%",
                            "50-100%", "NA")

colFill2 <- scale_fill_manual(name = "Tipping point", values = tippingpointPal2)

nolabels_theme <- theme(axis.title.x =element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        plot.title = element_text(size = 12, face = "plain"),
                        legend.position = "none")
onlyx_theme <- theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")
onlyy_theme <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")

theme_set(theme_minimal(base_size = 12))
# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
country    <- "USA"

C          <- readRDS(paste0("C_", country, "_bytens_overallupdated.RData"))
C_low      <- C/2 # scale to an R_0 of 1.3

age_demo   <- readRDS(paste0("age_demographics_", country,".RData"))

pop_total  <- age_demo[10]
age_demo   <- age_demo[1:9]

N_i        <- pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

IFR        <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 
                4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # Ref: Levin
IFR        <- IFR/100 # as decimal
YLL_vec    <- readRDS(paste0("yll_vec_", country, ".RData"))

# susceptibility with R0 =2.6 for BEL, ref: Davies
#u_constant <- rep(0.0154, 9) 
u_var      <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/(39.80957/2)
R0         <- compute_R0(u_var, C)

this_v_e   <- get_v_e(p = 0.9, y0 = 0.9, hinge_age = 50)
v_e_var    <- get_v_e(p = 0.5, y0 = 0.9, hinge_age = 50)
v_e_type   <- "aorn"

sero_none <- rep(0, 9) # no prior immunity
sero_belgium <- c(0.036, 0.137, 0.130, 0.131, 0.133, 0.132, 0.133, 
                  0.067, 0.0501) # ref: Herzog, weighted average of 80 and 90 bins to get 80+
sero_CT <- c(0.039, 0.0382, 0.031, 0.031, 0.031, 0.037, 0.032, 0.027, 0.027)
sero_NY <- c(0.32, 0.3129, 0.249, 0.249, 0.264, 0.279, 0.2575, 0.2215, 0.207) # ref: NYC

#### temporary for debugging ####
#df1 <- run_sim_new(C, 0.39, "adults", num_perday, v_e_type, u_var, this_v_e, sero = sero_NY)
#df2<- run_sim_new(C, 0.40, "adults", num_perday, v_e_type, u_var, this_v_e)
# df41 <- run_sim_new(C, 0.41, "adults", num_perday, v_e_type, u_var, this_v_e)
# # 
# print(compute_total_cases_new(df_36))
# print(compute_total_deaths_new(df_36))
# # 
# # # calculate the total number of people vaccinated
# row <- 366
# sum(df_36[row,11:19] + df_36[row,38:46] + df_36[row,65:73] + df_36[row,92:100])
# 
# # 
# # # calculate the total number of infected
# # sum(df2_36[row,83:109])
# sum(list_elderly$`26`[row,11:19] + list_elderly$`26`[row,38:46] + list_elderly$`26`[row,65:73] + list_elderly$`26`[row,92:100])


# _____________________________________________________________________
# FIGURE 1 ####
# Scenario 1: R0 = 1.3, continuous rollout at num_perday = 0.01
# Scenario 2: R0 = 2.6, continuous rollout at num_perday = 0.01
# Scenario 3: R0 = 2.6, anticipatory rollout (all at t = 0)
# _____________________________________________________________________
ptm <- proc.time()

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

scenarios <- 1:3
fig1_plots <- foreach (p = scenarios) %dopar% {
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
  # specify parameters for each scenario
  if (p == 1){
    this_C <- C_low
    num_perday <- 0.01
  } else if (p == 2){
    this_C <- C
    num_perday <- 0.01
  } else {
    this_C <- C
    num_perday <- 1
  }

  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  } 
  p_mort <- plot_over_vax_avail_new("deaths", "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus) 
  p_infect <- plot_over_vax_avail_new("cases", "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)
  p_yll <- plot_over_vax_avail_new("YLL", "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)

  list(p_mort, p_infect, p_yll)
}
stopCluster(cl)
proc.time() - ptm

mort_1 <- fig1_plots[[1]][[1]] + 
  onlyy_theme + 
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25), "cm"))
mort_2 <- fig1_plots[[2]][[1]] + 
  nolabels_theme
mort_3 <- fig1_plots[[3]][[1]] + 
  nolabels_theme

infect_1 <- fig1_plots[[1]][[2]]  + 
  theme(axis.title.x = element_blank())
infect_2 <- fig1_plots[[2]][[2]] + 
  onlyx_theme + 
  theme(axis.title.x = element_blank())
infect_3 <- fig1_plots[[3]][[2]] +
  onlyx_theme  + 
  theme(axis.title.x = element_blank())

panel <- ggarrange(mort_1, mort_2, mort_3, infect_1, infect_2, infect_3,
                   nrow = 2,
                   labels = c('B', 'C', 'D', 'E', 'F', 'G'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1.1),
                   bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12), 
                                     vjust = 0.3, hjust = 0.36),
                   padding = unit(0.5, "line"))

p1 <- barplot_vax_strat("kids") + 
  theme(axis.title.y = element_blank())
  #ylab("Distribution\nof vaccines (%)")
p2 <- barplot_vax_strat("adults") + 
  theme(axis.title.y = element_blank())
p3 <- barplot_vax_strat("20+") + 
  theme(axis.title.y = element_blank())
p4 <- barplot_vax_strat("elderly") + 
  theme(axis.title.y = element_blank())
p5 <- barplot_vax_strat("all") + 
  theme(axis.title.y = element_blank())

strategy_panel <- ggarrange(p1, p2, p3, p4, p5,
                            nrow = 5, 
                            labels = c('A',  '', '', '', ''),
                            label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                              hjust=1.8, vjust = 1.1),
                            left = textGrob("Distribution of vaccines (%)", rot = 90, hjust = 0.5),
                            bottom = textGrob("Age (years)", vjust = 0))
# export as pdf 9.5x4"
grid.arrange(strategy_panel, panel,
             ncol = 2, widths = c(2, 7.5),
             padding = unit(1, "line"))

# _____________________________________________________________________
# FIGURE 2: Age-dep ve ####
# Make sure to update so solid line is thinner in plot_over_vax_avail_new (size = 0.75)
# _____________________________________________________________________
ptm <- proc.time()

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

scenarios <- 1:3
fig2_plots <- foreach (p = scenarios) %dopar% {
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
  # specify parameters for each scenario
  if (p == 1){
    this_C <- C_low
    num_perday <- 0.01
  } else if (p == 2){
    this_C <- C
    num_perday <- 0.01
  } else {
    this_C <- C
    num_perday <- 1
  }

  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
    list_all_var[[paste0(i)]] <- run_sim_new(this_C, j, "all", num_perday, v_e_type, v_e_var)
    list_kids_var[[paste0(i)]] <- run_sim_new(this_C, j, "kids", num_perday, v_e_type, v_e_var)
    list_adults_var[[paste0(i)]] <- run_sim_new(this_C, j, "adults", num_perday, v_e_type, v_e_var)
    list_elderly_var[[paste0(i)]] <- run_sim_new(this_C, j, "elderly", num_perday, v_e_type, v_e_var)
    list_twentyplus_var[[paste0(i)]] <- run_sim_new(this_C, j, "twentyplus", num_perday, v_e_type, v_e_var)
  }
  p_mort <- plot_over_vax_avail_new("deaths", "None", list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)

  p_infect <- plot_over_vax_avail_new("cases", "None", list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)

  p_yll <- plot_over_vax_avail_new("YLL", "None", list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)

  list(p_mort, p_infect, p_yll)
}
stopCluster(cl)
proc.time() - ptm

# ignore warning about missing values
leg_col <- color_legend()
leg_ve <- ve_legend() 

mort_1 <- fig2_plots[[1]][[1]] +
  theme(legend.position = "none",
        axis.title.x = element_blank())
mort_2 <- fig2_plots[[2]][[1]]+
  xlab("Total vaccine supply (% of population)") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
mort_3 <- fig2_plots[[3]][[1]]+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# age-dependent ve plot
groups <- c(0,10,20,30,40,50,60,70,80,90)
v_e_plot <- v_e_var
v_e_plot[10] <- v_e_var[9]

df <- data.frame(groups, v_e_plot)

age_dep_ve <- ggplot(df, aes(x = groups, y = v_e_plot*100)) +
  geom_step(size = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 90, size = 0.5) +
  geom_point(aes(x = 60, y = 90), size = 2) +
  geom_point(aes(x = 88, y = 50), size = 2) +
  ylab("Efficacy (%)") +
  scale_y_continuous(expand = c(0,0), limit = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 90.2), breaks = c(0,20,40,60,80),
                     labels = c(0,20,40,60,80)) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Age (years)") +
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_text(vjust = 0.2)) + 
  coord_fixed(9*1.25/10)

age_dep_ve <- ggarrange(age_dep_ve,
                        labels = c('A'),
                        label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                          hjust=0, vjust = 1))

panel <- ggarrange(mort_1, mort_2, mort_3,
          nrow = 1, widths = c(2, 2, 2),
          labels = c('B', 'C', 'D'),
          label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                            hjust=0, vjust = 1),
          padding = unit(0.5, "line"))

# export as 9.5x2
grid.arrange(age_dep_ve, panel,
             # layout_matrix = lay, 
             widths = c(2, 7.5))

# _____________________________________________________________________
# FIGURE 3: NYC seroprevalence ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________
compute_R_with_sero(u_var, C, sero_NY)

ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
num_perday <- 0.01

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_NY)
  list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_NY)
  list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY)
  list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY)
  list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_NY)
  list_all_var[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY, sp = 0.975, se = 0.967)
  list_kids_var[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY, sp = 0.975, se = 0.967)
  list_adults_var[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY, sp = 0.975, se = 0.967)
  list_elderly_var[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_NY, sp = 0.975, se = 0.967)
  list_twentyplus_var[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_NY, sp = 0.975, se = 0.967)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail_new("deaths", "None", list_all_var, list_kids_var, 
                                  list_adults_var, list_elderly_var, list_twentyplus_var)

p_infect <- plot_over_vax_avail_new("cases", "None", list_all_var, list_kids_var, 
                                    list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())
            
p_yll <- plot_over_vax_avail_new("YLL", "None", list_all_var, list_kids_var, 
                                 list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())

panel <- ggarrange(p_infect, p_mort, p_yll,
          nrow = 1,
          labels = c('A', 'B', 'C'),
          label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                            hjust=0, vjust = 1),
          padding = unit(0.5, "line"))

leg_sero <- sero_legend()

lay <- rbind(c(1, 3),
             c(2, 3),
             c(NA,3))

# export as pdf 9.5x2
grid.arrange(leg_sero, leg_col, panel,
             layout_matrix = lay, 
             widths = c(0.9, 8.6),
             heights = c(.9, .9, .2))
# _____________________________________________________________________
# FIGURE S12: Synthetic seroprevalence (40%) for US ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________
IC_syn <- readRDS("synthetic_sero_US.RData")
sero_syn <- IC_syn[82:90]/N_i
print(compute_R_with_sero(u_var, C, sero_syn))

ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
num_perday <- 0.01

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type,
                                       this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type,
                                        this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type,
                                          this_v_e, u = u_var,syn_sero_compartments = IC_syn)
  list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type,
                                           this_v_e, u = u_var,syn_sero_compartments = IC_syn)
  list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type,
                                              this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_all_var[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type,
                                           this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = 0.975, se = 0.967)
  list_kids_var[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, 
                                            this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = 0.975, se = 0.967)
  list_adults_var[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, 
                                              this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = 0.975, se = 0.967)
  list_elderly_var[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type,
                                               this_v_e, u = u_var, syn_sero_compartments = IC_syn, sp = 0.975, se = 0.967)
  list_twentyplus_var[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type,
                                                  this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = 0.975, se = 0.967)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail_new("deaths", "None", list_all_var, list_kids_var, 
                                  list_adults_var, list_elderly_var, list_twentyplus_var)

p_infect <- plot_over_vax_avail_new("cases", "None", list_all_var, list_kids_var, 
                                    list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())

p_yll <- plot_over_vax_avail_new("YLL", "None", list_all_var, list_kids_var, 
                                 list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())

panel <- ggarrange(p_infect, p_mort, p_yll,
                   nrow = 1,
                   labels = c('A', 'B', 'C'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1),
                   padding = unit(0.5, "line"))

leg_sero <- sero_legend()

lay <- rbind(c(1, 3),
             c(2, 3),
             c(NA,3))

# export as pdf 9.5x2
grid.arrange(leg_sero, leg_col, panel,
             layout_matrix = lay, 
             widths = c(0.9, 8.6),
             heights = c(.9, .9, .2))

# _____________________________________________________________________
# FIGURE S12: Low seroprevalence for US - Connecticut sero ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________
print(compute_R_with_sero(u_var, C, sero_CT))

ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
num_perday <- 0.01

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT)
  list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT)
  list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_all_var[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = 0.975, se = 0.967)
  list_kids_var[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = 0.975, se = 0.967)
  list_adults_var[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = 0.975, se = 0.967)
  list_elderly_var[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT, sp = 0.975, se = 0.967)
  list_twentyplus_var[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = 0.975, se = 0.967)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail_new("deaths", "None", list_all_var, list_kids_var, 
                                  list_adults_var, list_elderly_var, list_twentyplus_var)

p_infect <- plot_over_vax_avail_new("cases", "None", list_all_var, list_kids_var, 
                                    list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())

p_yll <- plot_over_vax_avail_new("YLL", "None", list_all_var, list_kids_var, 
                                 list_adults_var, list_elderly_var, list_twentyplus_var)+
  theme(axis.title.x = element_blank())

panel <- ggarrange(p_infect, p_mort, p_yll,
                   nrow = 1,
                   labels = c('A', 'B', 'C'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1),
                   padding = unit(0.5, "line"))

leg_sero <- sero_legend()

lay <- rbind(c(1, 3),
             c(2, 3),
             c(NA,3))

# export as pdf 9.5x2
grid.arrange(leg_sero, leg_col, panel,
             layout_matrix = lay, 
             widths = c(0.9, 8.6),
             heights = c(.9, .9, .2))

# _____________________________________________________________________
# FIGURE S?: Tipping point heatmaps ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________

ptm <- proc.time()
v_e_bisection(1, 50, 30, 1, "aorn")
proc.time() - ptm

# loop over vaccine supply
hinge_age_vec <- c(50, 60, 70) # age group that v_e begins decreasing after
baseline_vec <- seq(1, 0.3, by = -0.1)

tp_val <- matrix(NA, 8, 51)
tp_strat <- matrix(NA, 8, 51)

ptm <- proc.time()
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

v_e_type <- "aorn"
this_hinge_age <- 70
list_tp_leaky <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 50, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], this_hinge_age, i, num_perday = 1, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]

    #tp_val[j, i] <- temp[[1]]
    #tp_strat[j, i]<- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)
proc.time() - ptm

#saveRDS(list_tp_leaky, "tp_val_aorn_79_s3_IND.RData")

# list_tp_leaky[[1]] corresponds to baseline ve 100%
# list_tp_leaky[[8]] corresponds to baseline ve 30%

list_79_2 <- readRDS("tp_val_leaky_79_s2.RData")
list_69_2 <- readRDS("tp_val_leaky_69_s2.RData")
list_59_2 <- readRDS("tp_val_leaky_59_s2.RData")

list_79_3 <- readRDS("tp_val_leaky_79_s3.RData")
list_69_3 <- readRDS("tp_val_leaky_69_s3.RData")
list_59_3 <- readRDS("tp_val_leaky_59_s3.RData")
 
p79_2 <- plot_tipping_point_heatmap_2(list_79_2) + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) 
p69_2 <- plot_tipping_point_heatmap_2(list_69_2)+ 
  ylab("Baseline ve (%)") + 
  theme(legend.position = "none",
        axis.text.x = element_blank())
p59_2 <- plot_tipping_point_heatmap_2(list_59_2) + 
  ggtitle("Continuous rollout") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "plain"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank())

p79_3 <- plot_tipping_point_heatmap_2(list_79_3) + 
  theme(axis.title.y = element_blank(), 
        legend.position = "none",
        axis.text.y = element_blank())
p69_3 <- plot_tipping_point_heatmap_2(list_69_3) + 
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ,
        legend.position = "none")
p59_3 <- plot_tipping_point_heatmap_2(list_59_3) + 
  ggtitle("Anticipatory rollout") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "plain"))

leg <- get_legend(p59_3)

p59_3 <- p59_3 + theme(legend.position = "none")

#grid.arrange(p59_2, p59_3, p69_2, p69_3, p79_2, p79_3,
#             ncol = 2)

p <- ggarrange(p59_2, p59_3, p69_2, p69_3, p79_2, p79_3,
          nrow = 3,
          # labels = c("A", "B", "C", "D ", "E ", "F "),
          # label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
          #                   hjust=0, vjust = 1),
          bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12),
                            vjust = 0.3, hjust = 0.36),
          padding = unit(0.5,
                         "line"))

# save as 5x7
t1 <- textGrob("Hinge age\n59")
t2 <- textGrob("Hinge age\n69")
t3 <- textGrob("Hinge age\n79")

lay <- rbind(c(1, 4, 5),
             c(2, 4, 5),
             c(3, 4, 5))

# export as 8x6
g <- grid.arrange(t1, t2, t3, p, leg,
             layout_matrix = lay, 
             widths = c(1, 6, 1.5))

# _____________________________________________________________________
# FIGURE S10: Non transmission blocking heatmaps ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________

country <- "BEL"

v_e <- list(rep(0.3, 9),
            rep(0.4, 9),
            rep(0.5, 9),
            rep(0.6, 9),
            rep(0.7, 9),
            rep(0.8, 9),
            rep(0.9, 9),
            rep(1.0, 9))

this_v_e <- rep(0.9, 9)
v_e_type <- "aorn"

IFR     <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # Ref: Levin
IFR     <- IFR/100 # as decimal
YLL_vec <- readRDS(paste0("yll_vec_", country, ".RData"))

#u_var   <- rep(0.0154, 9) # susceptibility with R0 =2.6 for BEL, ref: Davies
u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 

C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
C_low <- C/2 # scale to an R_0 of 1.3
R0 <- compute_R0(u_var, C)

scale_C <- seq(0.5, 1, by = 0.5/13)

sero_none <- rep(0, 9) # no prior immunity
num_perday <- 0.01

rollouts <- seq(0.0025, 0.02, by = 0.0025)

# scale the susceptibility for each country to get R0 = 2.6 for scenario 2 & 3 and R0 = 1.3 
# for scenario 1
countries_scale_u_26 <- c()
unscaled_u <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)
count <- 1
for (i in countries){
  this_C <- readRDS(paste0("C_", i, "_bytens_overall.RData"))
  countries_scale_u_26[count] <- scale_u_for_R0(unscaled_u, this_C, 2.6)
  count <- count + 1
}

# _____________________________________________________________________
# RUN SIMS ----
# _____________________________________________________________________
ptm <- proc.time()

# clusters
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#HM <- foreach(this_scale_C = scale_C) %dopar% {
#HM <- foreach(this_v_e = v_e) %dopar% {
HM <- foreach(num_perday = rollouts) %dopar% {
  # HM <- foreach(country = countries) %dopar% {
  
  C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
  
  # when looping over countries
  # index <- match(country, countries)
  # u_var <- unscaled_u/countries_scale_u_26[index]
  
  # when looping over R0
  #C <- C*this_scale_C
  
  age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
  pop_total <- age_demo[10]
  age_demo <- age_demo[1:9]
  N_i <-  pop_total*age_demo      
  num_groups <- length(age_demo) # num age groups
  
  # * SCENARIO 3: r0 - 2.6, anticipatory rollout
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  num_perday <- 1
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  } 
  
  x_vec <- seq(1, 50, by = 1)
  I3 <- sapply(x_vec, FUN = function(x) get_best_strat_cases_new(x,
                                                                 list_all,list_kids,list_adults,
                                                                 list_elderly,list_twentyplus))
  M3 <- sapply(x_vec, FUN = function(x) get_best_strat_deaths_new(x,
                                                                  list_all,list_kids,list_adults,
                                                                  list_elderly,list_twentyplus))
  Y3 <- sapply(x_vec, FUN = function(x) get_best_strat_yll_new(x,
                                                               list_all,list_kids,list_adults,
                                                               list_elderly,list_twentyplus))
  # * SCENARIO 2: r0 - 2.6, continuous rollout
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  
  num_perday <- 0.01
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  }
  
  I2 <- sapply(x_vec, FUN = function(x) get_best_strat_cases_new(x,
                                                                 list_all,list_kids,list_adults,
                                                                 list_elderly,list_twentyplus))
  M2 <- sapply(x_vec, FUN = function(x) get_best_strat_deaths_new(x,
                                                                  list_all,list_kids,list_adults,
                                                                  list_elderly,list_twentyplus))
  Y2 <- sapply(x_vec, FUN = function(x) get_best_strat_yll_new(x,
                                                               list_all,list_kids,list_adults,
                                                               list_elderly,list_twentyplus))
  
  # * SCENARIO 1: r0 - 1.3, continuous rollout
  C_low <- C/2
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C_low, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C_low, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C_low, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C_low, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C_low, j, "twentyplus", num_perday, v_e_type, this_v_e)
  }
  
  x_vec <- seq(1, 50, by = 1)
  I1 <- sapply(x_vec, FUN = function(x) get_best_strat_cases_new(x,
                                                                 list_all,list_kids,list_adults,
                                                                 list_elderly,list_twentyplus))
  M1 <- sapply(x_vec, FUN = function(x) get_best_strat_deaths_new(x,
                                                                  list_all,list_kids,list_adults,
                                                                  list_elderly,list_twentyplus))
  Y1 <- sapply(x_vec, FUN = function(x) get_best_strat_yll_new(x,
                                                               list_all,list_kids,list_adults,
                                                               list_elderly,list_twentyplus))
  list(s1 = cbind(I1, M1, Y1),
       s2 = cbind(I2, M2, Y2),
       s3 = cbind(I3, M3, Y3))
}
stopCluster(cl)

#saveRDS(HM, "HM_BEL_rollout_speed.RData")
proc.time() - ptm
