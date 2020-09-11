# Vaccine Strategy: Running multiple countries

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
# country codes:
# Belgium: BEL
# United States: USA
# India: IND
# Spain: ESP
# Zimbabwe: ZWE
# Brazil: BRA 
country <- "BEL"

C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
age_demo <- readRDS(paste0("age_demographics_", country,".RData"))

pop_total <- age_demo[10]
age_demo <- age_demo[1:9]

N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

IFR   <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 4.042049e-01, 1.355495e+00, 4.545632e+00,
           1.524371e+01) # Ref: Levin
IFR   <- IFR/100 # as decimal

# susceptibility
u_constant     <- rep(0.0154, 9) # constant # 0.02 for US, 0.022 for Belgium
u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 # R0 = 2.6 for BEL
R0 <- compute_R0(u_constant, C)

sero_none <- rep(0, 9) # no prior immunity

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________
v_e <- c(.75, 0.5)

ptm <- proc.time()
for (k in v_e){
  v_e_var <- get_v_e(p = k, y0 = k, hinge_age = 50)
  
  list_all      <- vector(mode = "list")
  list_kids     <- vector(mode = "list")
  list_adults   <- vector(mode = "list")
  list_elderly  <- vector(mode = "list")
  list_twentyplus   <- vector(mode = "list")
  
  for (i in seq(0, 50, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim(C, percent_vax = j, strategy = "all", u = u_constant, v_e_var)
    list_kids[[paste0(i)]] <- run_sim(C, j, "kids", u_constant, v_e_var)
    list_adults[[paste0(i)]] <- run_sim(C, j, "adults", u_constant, v_e_var)
    list_elderly[[paste0(i)]] <- run_sim(C, j, "elderly", u_constant, v_e_var)
    list_twentyplus[[paste0(i)]] <- run_sim(C, j, "20+", u_constant, v_e_var)
  }
  print(k)
  outcome <- "cases"
  p <- plot_over_vax_avail(outcome, "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)
  
  if (k == .75){
    p <- p + scale_y_continuous(expand = c(0,0), limit = c(0, 100))
  } else if (k == .50){
    p <- p + scale_y_continuous(expand = c(0,0), limit = c(0, 100))
  }
  print(p)
  
  outcome <- "deaths"
  p <- plot_over_vax_avail(outcome, "None", list_all, list_kids, list_adults, list_elderly, list_twentyplus)
  print(p)
}
proc.time() - ptm
