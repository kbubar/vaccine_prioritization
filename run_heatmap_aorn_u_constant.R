# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________
library(tidyverse) 
library(deSolve) # ode solver
library(gridExtra)
library(RColorBrewer)
library(egg)
library(foreach)
library(doParallel)
library(gplots)

# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green (20-49)
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple (60+)
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)

# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# country codes:
country <- "USA"
v_e_type <- "aorn"

u_var     <- rep(0.01563, 9) # US R0 = 2.6

v_e <- list(rep(0.3, 9),
            rep(0.4, 9),
            rep(0.5, 9),
            rep(0.6, 9),
            rep(0.7, 9),
            rep(0.8, 9),
            rep(0.9, 9),
            rep(1.0, 9))

IFR     <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 
             4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # Ref: Levin
IFR     <- IFR/100 # as decimal
YLL_vec <- readRDS(paste0("yll_vec_", country, ".RData"))

C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
C_low <- C/2 # scale to an R_0 of 1.3
R0 <- compute_R0(u_var, C)

age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total <- age_demo[10]
age_demo <- age_demo[1:9]
N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

sero_none <- rep(0, 9) # no prior immunity
num_perday <- 0.01

ptm <- proc.time()

# clusters
cores=detectCores()
cl <- makeCluster(cores[1]-1) # to not overload your computer
registerDoParallel(cl)

HM <- foreach(this_v_e = v_e) %dopar% {
  
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

saveRDS(HM, "HM_USA_aorn_u_constant.RData")
proc.time() - ptm

# _____________________________________________________________________
# PLOT ----
# _____________________________________________________________________
param <- "ve"

num <- length(HM)

S1_I <- S2_I <- S3_I <- S1_M <- S2_M <- S3_M <- S1_Y <- S2_Y <- S3_Y <- matrix(NA, num, 50)

# Get data in right form for plotting
for (i in 1:num){
  S1_I[i ,] <- HM[[i]]$s1[, 1]
  S2_I[i ,] <- HM[[i]]$s2[, 1]
  S3_I[i ,] <- HM[[i]]$s3[, 1]
  
  S1_M[i ,] <- HM[[i]]$s1[, 2]
  S2_M[i ,] <- HM[[i]]$s2[, 2]
  S3_M[i ,] <- HM[[i]]$s3[, 2]
  
  S1_Y[i ,] <- HM[[i]]$s1[, 3]
  S2_Y[i ,] <- HM[[i]]$s2[, 3]
  S3_Y[i ,] <- HM[[i]]$s3[, 3]
}

S1_Inew <- S2_Inew <- S3_Inew <- S1_Mnew <- S2_Mnew <- S3_Mnew <- S1_Ynew <- S2_Ynew <- S3_Ynew <- matrix(NA, num, 50)

#adults = 1 (20-49)
#twentyplus = 2
#elderly = 3
#kids = 4
#all = 5
### I Matrices
for(i in 1){
  S1_Inew[S1_I == "adults"] = 1
  S1_Inew[S1_I == "twentyplus"] = 2
  S1_Inew[S1_I == "elderly"] = 3
  S1_Inew[S1_I == "kids"] = 4
  S1_Inew[S1_I == "all"] = 5
  
  S2_Inew[S2_I == "adults"] = 1
  S2_Inew[S2_I == "twentyplus"] = 2
  S2_Inew[S2_I == "elderly"] = 3
  S2_Inew[S2_I == "kids"] = 4
  S2_Inew[S2_I == "all"] = 5
  
  S3_Inew[S3_I == "adults"] = 1
  S3_Inew[S3_I == "twentyplus"] = 2
  S3_Inew[S3_I == "elderly"] = 3
  S3_Inew[S3_I == "kids"] = 4
  S3_Inew[S3_I == "all"] = 5
  
  ### M matrices
  S1_Mnew[S1_M == "adults"] = 1
  S1_Mnew[S1_M == "twentyplus"] = 2
  S1_Mnew[S1_M == "elderly"] = 3
  S1_Mnew[S1_M == "kids"] = 4
  S1_Mnew[S1_M == "all"] = 5
  
  S2_Mnew[S2_M == "adults"] = 1
  S2_Mnew[S2_M == "twentyplus"] = 2
  S2_Mnew[S2_M == "elderly"] = 3
  S2_Mnew[S2_M == "kids"] = 4
  S2_Mnew[S2_M == "all"] = 5
  
  S3_Mnew[S3_M == "adults"] = 1
  S3_Mnew[S3_M == "twentyplus"] = 2
  S3_Mnew[S3_M == "elderly"] = 3
  S3_Mnew[S3_M == "kids"] = 4
  S3_Mnew[S3_M == "all"] = 5
  
  ### Y matrices
  S1_Ynew[S1_Y == "adults"] = 1
  S1_Ynew[S1_Y == "twentyplus"] = 2
  S1_Ynew[S1_Y == "elderly"] = 3
  S1_Ynew[S1_Y == "kids"] = 4
  S1_Ynew[S1_Y == "all"] = 5
  
  S2_Ynew[S2_Y == "adults"] = 1
  S2_Ynew[S2_Y == "twentyplus"] = 2
  S2_Ynew[S2_Y == "elderly"] = 3
  S2_Ynew[S2_Y == "kids"] = 4
  S2_Ynew[S2_Y == "all"] = 5
  
  S3_Ynew[S3_Y == "adults"] = 1
  S3_Ynew[S3_Y == "twentyplus"] = 2
  S3_Ynew[S3_Y == "elderly"] = 3
  S3_Ynew[S3_Y == "kids"] = 4
  S3_Ynew[S3_Y == "all"] = 5
}

## ve
pS1_I <- plot_heatmap(reshape2::melt(S1_Inew), num, param) +
  ylab("Cumulative\nincidence\n\n") + 
  onlyy_theme +
  ggtitle("Scenario 1")  
pS2_I <- plot_heatmap(reshape2::melt(S2_Inew), num, param)  +
  nolabels_theme +
  ggtitle("Scenario 2")
pS3_I <- plot_heatmap(reshape2::melt(S3_Inew), num, param)  +
  nolabels_theme +
  ggtitle("Scenario 3")

pS1_M <- plot_heatmap(reshape2::melt(S1_Mnew), num, param)+
  ylab("Mortality\n\nVaccine efficacy (%)") + 
  onlyy_theme
pS2_M <- plot_heatmap(reshape2::melt(S2_Mnew), num, param) +
  nolabels_theme
pS3_M <- plot_heatmap(reshape2::melt(S3_Mnew), num, param) +
  nolabels_theme

pS1_Y <- plot_heatmap(reshape2::melt(S1_Ynew), num, param) +
  ylab("Years of\nlife lost\n\n")
pS2_Y <- plot_heatmap(reshape2::melt(S2_Ynew), num, param) +
  onlyx_theme
pS3_Y <- plot_heatmap(reshape2::melt(S3_Ynew), num, param) +
  onlyx_theme

# export as 9" x 6"
g <- ggarrange(pS1_I, pS2_I, pS3_I, pS1_M, pS2_M, pS3_M, pS1_Y, pS2_Y, pS3_Y,
               bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12)),
               padding = unit(0.5, "line"))

pdf("C:/Users/bubar/Documents/Vaccine Strategy/Plots/finaldraft_plots/heatmap_aorn_u_constant.pdf",
    height = 6, width = 8.5)
g
dev.off()