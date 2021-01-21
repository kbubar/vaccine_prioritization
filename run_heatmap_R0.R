# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")
country <- "USA"
source("setup.R")
# _____________________________________________________________________
# RUN SIM ----
# Create heatmap for a range of R0 values
# _____________________________________________________________________
# scale C to get range of R0
scale_C <- c()
R0 <- seq(1.1, 2, by = 0.05)
count <- 1
for (i in R0){
  scale_C[count] <- scale_u_for_R0(u_var, C, i)
  count <- count + 1
}

ptm <- proc.time()

# clusters
cores=detectCores()
cl <- makeCluster(cores[1]-1) # to not overload your computer
registerDoParallel(cl)

HM <- foreach(this_scale_C = scale_C) %dopar% {
  C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
  this_C <- C/this_scale_C
  
  # * Continuous rollout
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  
  num_perday <- 0.001
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  }
  x_vec <- seq(1, 50, by = 1)
  I1 <- sapply(x_vec, FUN = function(x) get_best_strat_cases_new(x, list_all,list_kids,list_adults,
                                                                 list_elderly,list_twentyplus))
  M1 <- sapply(x_vec, FUN = function(x) get_best_strat_deaths_new(x, list_all,list_kids,list_adults,
                                                                  list_elderly,list_twentyplus))
  Y1 <- sapply(x_vec, FUN = function(x) get_best_strat_yll_new(x, list_all,list_kids,list_adults,
                                                               list_elderly,list_twentyplus))
  
  # * Anticipatory rollout
  num_perday <- 1
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(this_C, j, "all", num_perday, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim(this_C, j, "kids", num_perday, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim(this_C, j, "adults", num_perday, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
  }
  
  x_vec <- seq(1, 50, by = 1)
  I2 <- sapply(x_vec, FUN = function(x) get_best_strat_cases_new(x, list_all,list_kids,list_adults,
                                                                 list_elderly,list_twentyplus))
  M2 <- sapply(x_vec, FUN = function(x) get_best_strat_deaths_new(x, list_all,list_kids,list_adults,
                                                                  list_elderly,list_twentyplus))
  Y2 <- sapply(x_vec, FUN = function(x) get_best_strat_yll_new(x, list_all,list_kids,list_adults,
                                                               list_elderly,list_twentyplus))
  list(s1 = cbind(I1, M1, Y1),
       s2 = cbind(I2, M2, Y2))
}
stopCluster(cl)

saveRDS(HM, "HM_USA_R0_001.RData")
proc.time() - ptm

# _____________________________________________________________________
# PLOT ----
# _____________________________________________________________________
param <- "R0"

num <- length(HM)

S1_I <- S2_I <- S1_M <- S2_M <- S1_Y <- S2_Y <- matrix(NA, num, 50)

# Get data in right form for plotting
for (i in 1:num){
  S1_I[i ,] <- HM[[i]]$s1[, 1]
  S2_I[i ,] <- HM[[i]]$s2[, 1]
  
  S1_M[i ,] <- HM[[i]]$s1[, 2]
  S2_M[i ,] <- HM[[i]]$s2[, 2]
  
  S1_Y[i ,] <- HM[[i]]$s1[, 3]
  S2_Y[i ,] <- HM[[i]]$s2[, 3]
}

S1_Inew <- S2_Inew <- S1_Mnew <- S2_Mnew <- S1_Ynew <- S2_Ynew <- matrix(NA, num, 50)

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
}

## R0
pS1_I <- plot_heatmap(reshape2::melt(S1_Inew), num, param) +
  ylab("Cumulative\nincidence\n\n") + 
  onlyy_theme +
  ggtitle("Continuous Rollout")  
pS2_I <- plot_heatmap(reshape2::melt(S2_Inew), num, param)  +
  nolabels_theme +
  ggtitle("Anticipatory Rollout")

pS1_M <- plot_heatmap(reshape2::melt(S1_Mnew), num, param)+
  ylab("Mortality\n\nR0") + 
  onlyy_theme
pS2_M <- plot_heatmap(reshape2::melt(S2_Mnew), num, param) +
  nolabels_theme

pS1_Y <- plot_heatmap(reshape2::melt(S1_Ynew), num, param) +
  ylab("Years of\nlife lost\n\n")
pS2_Y <- plot_heatmap(reshape2::melt(S2_Ynew), num, param) +
  onlyx_theme

# export as 8.5" x 6"
g <- ggarrange(pS1_I, pS2_I, pS1_M, pS2_M, pS1_Y, pS2_Y,
               bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12)),
               padding = unit(0.5, "line"))
