# Get synthetic seroprevalence from simulation
# for Supp Fig 12: US, R0 = 2.6, seroprevalence = 40%

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

source("run_sim.R")
source("helper_functions.R")

country <- "USA"
C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total <- age_demo[10]
age_demo <- age_demo[1:9]
N_i <-  pop_total*age_demo  

IFR     <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 
             4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # Ref: Levin
IFR     <- IFR/100 # as decimal

age_demo <- readRDS(paste0("age_demographics_", country,".RData"))

pop_total <- age_demo[10]
age_demo <- age_demo[1:9]

N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

u_var      <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/39.80957
v_e_constant <- get_v_e(p = 1, y0 = 1, hinge_age = 50)
sero_none <- rep(0, 9) # no prior immunity

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________

dat <- run_sim_new(C, percent_vax = 0, strategy = "all", num_perday = 0.01, v_e_type = "aorn")

R <- dat[,83:91]
R_tot <- rowSums(R)

R_tot <- floor((R_tot/pop_total)*100)

# store indices for R ~40%
val <- max(R_tot[R_tot <= 40])
index <- match(c(val), R_tot)

synthetic_sero <- dat[index, -1]
row.names(synthetic_sero) <- 1

saveRDS(synthetic_sero, "synthetic_sero_US.RData")
