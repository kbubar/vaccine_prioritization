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


# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green (20-49)
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple (60+)
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)

tippingpointPal2 <- c("#E6E6E6", "#FFEDA0", "#FEB24C", "#F03B20", "#000000")
names(tippingpointPal2) <- c("None", "0-25%", "25-50%",
                             "50-100%", "NA")

colFill2 <- scale_fill_manual(name = "Tipping point", values = tippingpointPal2)

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

C          <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
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

u_var      <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/39.80957
R0         <- compute_R0(u_var, C)

sero_none <- rep(0, 9) # no prior immunity

# _____________________________________________________________________
# FIND TP FOR HEATMAP ----
# looping over vax supply
# _____________________________________________________________________

tp_val_59_s3 <- tp_strat_59_s3 <- tp_val_69_s3 <- tp_strat_69_s3 <- tp_val_79_s3 <- tp_strat_79_s3 <- matrix(NA, 8, 51)
tp_val_59_s2 <- tp_strat_59_s2 <- tp_val_69_s2 <- tp_strat_69_s2 <- tp_val_79_s2 <- tp_strat_79_s2 <- matrix(NA, 8, 51)

baseline_vec <- seq(1, 0.3, by = -0.1)

# RUN AORN, SCENARIO 2 ####
v_e_type <- "aorn"
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

list_79_2 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 2, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 70, i, 0.01, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

###
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

list_69_2 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 50, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 60, i, 0.01, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

###
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

list_59_2 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 2, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 50, i, 0.01, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

# RUN AORN, SCENARIO 3 ####
v_e_type <- "aorn"
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

list_79_3 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 2, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 70, i, 1, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

###
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

list_69_3 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 50, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 60, i, 1, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

###
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

list_59_3 <- foreach (j = 1:length(baseline_vec)) %dopar% {
  tp_val   <- c()
  tp_strat <- c()
  for (i in seq(0, 2, by = 1)){
    temp <- v_e_bisection(baseline_vec[j], 50, i, 1, v_e_type)
    tp_val[i] <- temp[[1]]
    tp_strat[i] <- temp[[2]]
  }
  out <- list(tp_val, tp_strat)
}
stopCluster(cl)

# _____________________________________________________________________
# PLOT ----
# _____________________________________________________________________

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

p <- ggarrange(p59_2, p59_3, p69_2, p69_3, p79_2, p79_3,
               nrow = 3,
               bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12),
                                 vjust = 0.3, hjust = 0.36),
               padding = unit(0.5,
                              "line"))

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