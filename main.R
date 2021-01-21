# Vaccine Strategy Simple Model
#setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")
country <- "USA"
source("setup.R")

# _____________________________________________________________________
# FIGURE 1: Dynamics curves ####
# _____________________________________________________________________
this_C <- C/scale_15
ptm <- proc.time()  
for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
  list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
  list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
  list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
  list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
} 
proc.time() - ptm

p_mort <- plot_over_vax_avail("deaths")
p_infect <- plot_over_vax_avail("cases")

t_one <- 11
infect_10 <- plot_strat_overtime("I", list_all[[1]], list_all[[t_one]], list_adults[[t_one]], 
                                 list_kids[[t_one]], list_twentyplus[[t_one]], list_elderly[[t_one]], 0.1/num_perday) +
  onlyy_theme + 
  ggtitle("10% vaccine supply") + 
  theme(plot.title = element_text(color = "black"))
t_two <- 41
infect_30 <- plot_strat_overtime("I", list_all[[1]], list_all[[t_two]], list_adults[[t_two]], 
                                 list_kids[[t_two]], list_twentyplus[[t_two]], list_elderly[[t_two]], 0.3/num_perday) + 
  nolabels_theme +
  ggtitle("30% vaccine supply") + 
  theme(plot.title = element_text(color = "black"))

mort_10 <- plot_strat_overtime("D", list_all[[1]], list_all[[t_one]], list_adults[[t_one]], 
                               list_kids[[t_one]], list_twentyplus[[t_one]], list_elderly[[t_one]], 0.1/num_perday) 

mort_30 <- plot_strat_overtime("D", list_all[[1]], list_all[[t_two]], list_adults[[t_two]], 
                               list_kids[[t_two]], list_twentyplus[[t_two]], list_elderly[[t_two]], 0.3/num_perday) + 
  onlyx_theme

strategy_panel <- plot_vaxdist_hist()

sub_panel2 <- ggarrange(infect_10, infect_30, p_infect + onlyy_theme,
                        mort_10, mort_30, p_mort + xlab("Total vaccine supply\n(% of population)"),
                        nrow = 2)

# export as 9.5x4   
grid.arrange(strategy_panel, sub_panel2,
             widths = c(2, 7.5))

# _____________________________________________________________________
# FIGURE 1E and I: Reduction in infections & deaths ####
# Scenario 1: R0 = 1.15, continuous rollout at num_perday = 0.2%/day
# _____________________________________________________________________
ptm <- proc.time()

cores=detectCores()
cl <- makeCluster(cores[1]-1) # to not overload your computer
registerDoParallel(cl)

scenarios <- 1:2
fig1_plots <- foreach (p = scenarios) %dopar% {
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  # specify parameters for each scenario
  if (p == 1){
    this_C <- C/scale_115
  } else if (p == 2){
    this_C <- C/scale_15
  }
  
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  } 
  p_mort <- plot_over_vax_avail("deaths")
  p_infect <- plot_over_vax_avail("cases")
  p_yll <- plot_over_vax_avail("YLL")
  list(p_mort, p_infect, p_yll)
}
stopCluster(cl)

mort_1 <- fig1_plots[[1]][[1]]  + 
  theme(axis.title.x = element_blank()) +
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25), "cm"))
mort_2 <- fig1_plots[[2]][[1]]  + 
  onlyx_theme +
  theme(axis.title.x = element_blank())

infect_1 <- fig1_plots[[1]][[2]]  + 
  onlyy_theme
infect_2 <- fig1_plots[[2]][[2]] + 
  nolabels_theme

# export as pdf 5x4"
panel <- ggarrange(infect_1, mort_1,  
                   nrow = 2,
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1.1),
                   bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12), 
                                     vjust = 0.3, hjust = 0.36),
                   padding = unit(0.5, "line"))

proc.time() - ptm

# _____________________________________________________________________
# FIGURE 3: Age-dep ve ####
# Update plot_over_vax_avail so solid line is thinner than dashed line (size = 0.75)
# _____________________________________________________________________
ptm <- proc.time()

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

scenarios <- 1:2
fig2_plots <- foreach (p = scenarios) %dopar% {
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
  # specify parameters for each scenario
  if (p == 1){
    this_C <- C/scale_115
  } else if (p == 2){
    this_C <- C/scale_15
  }

  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
    list_all_var[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, v_e_var)
    list_kids_var[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, v_e_var)
    list_adults_var[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, v_e_var)
    list_elderly_var[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, v_e_var)
    list_twentyplus_var[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, v_e_var)
  }

  p_mort <- plot_over_vax_avail("deaths", TRUE)

  p_infect <- plot_over_vax_avail("cases", TRUE)

  p_yll <- plot_over_vax_avail("YLL", TRUE)

  list(p_mort, p_infect, p_yll)
}
stopCluster(cl)
proc.time() - ptm

mort_1 <- fig2_plots[[1]][[1]] +
  theme(legend.position = "none",
        axis.title.x = element_blank())
mort_2 <- fig2_plots[[2]][[1]]+
  onlyx_theme +
  theme(axis.title.x = element_blank())

age_dep_ve <- plot_age_dep_ve()

# export as 7.5x2
panel <- ggarrange(age_dep_ve, mort_1, mort_2,
          nrow = 1, widths = c(2, 2, 2),
          labels = c('A', 'B', 'C'),
          label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                            hjust=0, vjust = 1),
          padding = unit(0.5, "line"))

# _____________________________________________________________________
# FIGURE 4: NYC seroprevalence ####
# scenario 2 only: R0 = 1.5, continuous rollout with num_perday = 0.002
# _____________________________________________________________________
this_scale <- scale_u_for_R0(u_var, C, 2.038) # so that realized R is 1.5
print(compute_R_with_sero(u_var, C/this_scale, sero_NY))

ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")

this_C <- C/this_scale
this_sp <- 0.99
this_se <- 0.70

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e, sero = sero_NY)
  list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e, sero = sero_NY)
  list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e, sero = sero_NY)
  list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e, sero = sero_NY)
  list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e,  sero = sero_NY)
  list_all_var[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e, sero = sero_NY, sp = this_sp, se = this_se)
  list_kids_var[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e, sero = sero_NY, sp = this_sp, se = this_se)
  list_adults_var[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e, sero = sero_NY, sp = this_sp, se = this_se)
  list_elderly_var[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e, sero = sero_NY, sp = this_sp, se = this_se)
  list_twentyplus_var[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e, sero = sero_NY, sp = this_sp, se = this_se)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail("deaths", TRUE)

p_infect <- plot_over_vax_avail("cases", TRUE) +
  theme(axis.title.x = element_blank())
            
p_yll <- plot_over_vax_avail("YLL", TRUE)+
  theme(axis.title.x = element_blank())

panel <- ggarrange(p_infect, p_mort, p_yll,
          nrow = 1,
          labels = c('A', 'B', 'C'),
          label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                            hjust=0, vjust = 1),
          padding = unit(0.5, "line"))

# _____________________________________________________________________
# FIGURE S10: Synthetic seroprevalence (40%) for US ####
# scenario 2 only: R0 = 2.6, continuous rollout with num_perday = 0.01
# _____________________________________________________________________
IC_syn <- readRDS("synthetic_sero_US.RData")
sero_syn <- IC_syn[82:90]/N_i + IC_syn[100:108]/N_i
this_scale <- scale_u_for_R0(u_var, C, 2.6) # so that realized R is 1.43
print(compute_R_with_sero(u_var, C/this_scale, sero_syn))

this_C <- C/this_scale
ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
num_perday <- 0.002

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type,
                                       this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type,
                                        this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type,
                                          this_v_e, u = u_var,syn_sero_compartments = IC_syn)
  list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type,
                                           this_v_e, u = u_var,syn_sero_compartments = IC_syn)
  list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type,
                                              this_v_e, u = u_var, syn_sero_compartments = IC_syn)
  list_all_var[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type,
                                           this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = this_sp, se = this_se)
  list_kids_var[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, 
                                            this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = this_sp, se = this_se)
  list_adults_var[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, 
                                              this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = this_sp, se = this_se)
  list_elderly_var[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type,
                                               this_v_e, u = u_var, syn_sero_compartments = IC_syn, sp = this_sp, se = this_se)
  list_twentyplus_var[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type,
                                                  this_v_e, u = u_var,syn_sero_compartments = IC_syn, sp = this_sp, se = this_se)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail("deaths", TRUE)

p_infect <- plot_over_vax_avail("cases", TRUE) +
  theme(axis.title.x = element_blank())

p_yll <- plot_over_vax_avail("YLL", TRUE)+
  theme(axis.title.x = element_blank())

# export as pdf 9.5x2
panel <- ggarrange(p_infect, p_mort, p_yll,
                   nrow = 1,
                   labels = c('A', 'B', 'C'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1),
                   padding = unit(0.5, "line"))

# _____________________________________________________________________
# FIGURE S9: Low seroprevalence for US - Connecticut sero ####
# scenario 2 only: R0 = 1.5, continuous rollout with num_perday = 0.002
# _____________________________________________________________________
this_scale <- scale_u_for_R0(u_var, C, 1.5508) # so that realized R is 1.5
print(compute_R_with_sero(u_var, C/this_scale, sero_CT))

ptm <- proc.time()

list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
num_perday <- 0.002
this_C <- C/this_scale

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT)
  list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT)
  list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT)
  list_all_var[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = this_sp, se = this_se)
  list_kids_var[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = this_sp, se = this_se)
  list_adults_var[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = this_sp, se = this_se)
  list_elderly_var[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e, u = u_var, sero = sero_CT, sp = this_sp, se = this_se)
  list_twentyplus_var[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e, u = u_var,sero = sero_CT, sp = this_sp, se = this_se)
}
proc.time() - ptm
p_mort <- plot_over_vax_avail("deaths", TRUE)

p_infect <- plot_over_vax_avail("cases", TRUE) +
  theme(axis.title.x = element_blank())

p_yll <- plot_over_vax_avail("YLL", TRUE)+
  theme(axis.title.x = element_blank())

panel <- ggarrange(p_infect, p_mort, p_yll,
                   nrow = 1,
                   labels = c('A', 'B', 'C'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1),
                   padding = unit(0.5, "line"))

# _____________________________________________________________________
# FIGURE S8: Tipping point heatmaps ####
# scenario 2 only: R0 = 1.5, continuous rollout with num_perday = 0.002
# _____________________________________________________________________

# loop over vaccine supply
baseline_vec <- seq(1, 0.3, by = -0.1)

tp_val <- matrix(NA, 8, 51)
tp_strat <- matrix(NA, 8, 51)

this_C <- C_26

ptm <- proc.time()
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

v_e_type <- "aorn"
this_hinge_age <- 50 # 50, 60, or 70 : age group that v_e begins decreasing after
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

saveRDS(list_tp_leaky, "tp_aorn_59_A_ESP.RData")

# list_tp_leaky[[1]] corresponds to baseline ve 100%
# list_tp_leaky[[8]] corresponds to baseline ve 30%

list_79_2 <- readRDS("tp_aorn_79_C_IND.RData")
list_69_2 <- readRDS("tp_aorn_69_C_IND.RData")
list_59_2 <- readRDS("tp_aorn_59_C_IND.RData")

list_79_3 <- readRDS("tp_aorn_79_A_IND.RData")
list_69_3 <- readRDS("tp_aorn_69_A_IND.RData")
list_59_3 <- readRDS("tp_aorn_59_A_IND.RData")
 
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
# FIGURE S13-15: Dynamics curves ####
# scenario 2 only: R0 = 1.5, continuous rollout with num_perday = 0.002
# _____________________________________________________________________
num_perday <- 0.002
this_C <- C/scale_15

for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
  list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
  list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
  list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
  list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
} 
# export 6*5.5
infect <- plot_supp_dynamics_panel("I")
cumul_infect <- plot_supp_dynamics_panel("R")
cumul_mort <- plot_supp_dynamics_panel("D")

