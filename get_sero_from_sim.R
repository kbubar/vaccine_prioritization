# Get seroprevalence from simulation
# for serology paper
setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

country <- "USA"
C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total <- age_demo[10]
age_demo <- age_demo[1:9]
N_i <-  pop_total*age_demo  

age_demo <- readRDS(paste0("age_demographics_", country,".RData"))

pop_total <- age_demo[10]
age_demo <- age_demo[1:9]

N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 # R0 = 2.6 for BEL
v_e_constant <- get_v_e(p = 1, y0 = 1, hinge_age = 50)
sero_none <- rep(0, 9) # no prior immunity

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________
# strategy options: "no vax", "all", "kids", "adults", "elderly"
# df_novax <- run_sim(C, 0, "no vax", u_var, v_e_constant)
dat <- readRDS("US_sim_novax.RData")
R <- dat[,29:37]
R_tot <- rowSums(R)

R_tot <- floor((R_tot/pop_total)*100)

# store indices for R = 5,10,15...50%
indices <- match(c(5,10,15,20,25,30,35,40,45,50), R_tot)
indices[3] <- 101
indices[7] <- 109 # 35
indices[8] <- 111 # 40
indices[9] <- 113 # 45 

sero <- R[indices,]

for (i in 1:10){
  sero[i,] <- sero[i,]/N_i*100
}
rownames(sero) <- c("5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%", "45%", "50%")

sero <- round(sero, 2)
write.csv(sero, "US_sero_from_sim.csv")
