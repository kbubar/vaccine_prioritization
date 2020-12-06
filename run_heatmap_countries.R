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
library(gplots)

# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green (20-49)
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple (60+)
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)

theme_set(theme_minimal(base_size = 12))

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
# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("run_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# country codes:
countries <- c("BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL")

this_v_e <- rep(0.9, 9)
v_e_type <- "aorn"

IFR     <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 
             4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # Ref: Levin
IFR     <- IFR/100 # as decimal
YLL_vec <- readRDS(paste0("yll_vec_", country, ".RData"))

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

ptm <- proc.time()

# clusters
cores=detectCores()
cl <- makeCluster(cores[1]-1) # to not overload your computer
registerDoParallel(cl)

HM <- foreach(country = countries) %dopar% {
  
  C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
  
  # when looping over countries
  index <- match(country, countries)
  u_var <- unscaled_u/countries_scale_u_26[index]

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

saveRDS(HM, "HM_countries.RData")
proc.time() - ptm


