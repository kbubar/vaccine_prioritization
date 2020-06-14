setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________
library(tidyverse) 
library(deSolve) # ode solver
library(xlsx)
library(gridExtra)
library(RColorBrewer)
library(nloptr) # nonlinear optimization

# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________
source("optimize_sim.R")
source("helper_functions.R")

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# Demographics
C <- readRDS("C_USA_bytens_all.RData")
C_noschool <- readRDS("C_USA_bytens_noschool.RData")

age_demo <- read.xlsx("USA_demographics.xlsx", "data_bytens", header = FALSE, row.names = TRUE)
age_demo <- as.vector(age_demo[[1]])/100

pop_total <- 100000
N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

IFR   <- c(0.001, 0.001, 0.007, 0.02, 0.06, 0.2, 0.9, 2.4, 10.1) # Ref: Salje
IFR   <- IFR/100 # as decimal

# susceptibility with R0 ~ 3 (R0 <- compute_R0(u))
u_constant     <- rep(0.02, 9) 
u_var     <- c(0.33, 0.37, 0.69, 0.81, 0.74, 0.8, 0.89, 0.77, 0.77)/30 # Ref: Davies

# vaccine efficacy
v_e_constant <- get_v_e(p = 1)
v_e_var <- get_v_e(p = 0.5)

# _____________________________________________________________________
# RUN SIM ----
# _____________________________________________________________________
percent_vax <- 0.05
nvax <- percent_vax*pop_total
initial_vax <- rep(nvax/num_groups, 9) # start with uniform distribution

constraints <- function(x){
  # function defining the inequality constraints s.t. constraints >= 0 for all components
  h <- numeric(3)
  h[1] <- x[1]
  h[2] <- x[2]
  h[3] <- x[3]
  h[4] <- x[4]
  h[5] <- x[5]
  h[6] <- x[6]
  h[7] <- x[7]
  h[8] <- x[8]
  h[9] <- x[9]
  h[10] <- sum(x) - nvax
  h[11] <- nvax - sum(x)
  
  return(h)
}

to_minimize <- "deaths"

ptm <- proc.time()

optimal_vax <- cobyla(initial_vax, optimize_sim, hin = constraints,
            nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))

proc.time() - ptm

vax <- round(optimal_vax$par/N_i * 100, 2)
groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
test <- data.frame(groups, vax)

ggplot(test, aes(y=vax, x=groups)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Age group") + 
  ylab("Number vaccinated\n(% of age group)")+
  ylim(0, 100) + 
  ggtitle("Optimizing for least deaths")

test2 <- optimize_sim(optimal_vax$par)
p1 <- plot_allages_onestrategy(test2, "I")
print(p1)

#######

to_minimize <- "deaths"

ptm <- proc.time()

minimize_cases_df <- data.frame(matrix(ncol = 11, nrow = 50))
colnames(minimize_cases_df) <- c("Percent_vax", "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+", "Percent infected")
for (i in 1:50){
  percent_vax <- i/100
  nvax <- percent_vax*pop_total
  initial_vax <- rep(nvax/num_groups, 9)
 
  optimal_vax <- cobyla(initial_vax, optimize_sim, hin = constraints,
                        nl.info = FALSE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  minimize_cases_df[i, 2:10] <- round(optimal_vax$par/N_i * 100, 2)
  minimize_cases_df[i,1] <- i
  minimize_cases_df[i, 11] <- optimal_vax$value
}
proc.time() - ptm

saveRDS(minimize_cases_df, "optimal_simplemodel_deaths.RData")
