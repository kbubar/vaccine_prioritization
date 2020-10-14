# Comparing vaccinating oldest first vs infection minimizing distribution

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
library(wesanderson)

source("run_sim.R")
source("helper_functions.R")

# scenario initialization 
country <- "POL"
C <- readRDS(paste0("C_", country, "_bytens_overall.RData"))
age_demo <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total <- age_demo[10]
age_demo <- age_demo[1:9]
N_i <-  pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

IFR   <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 4.042049e-01, 1.355495e+00, 4.545632e+00,
           1.524371e+01)/100 # Ref: Levin

u_constant     <- rep(0.0154, 9) # constant # 0.0154 for Belgium
u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 # R0 = 2.6 for BEL
R0 <- compute_R0(u_var, C)
v_e_constant <- get_v_e(p = 1, y0 = 1, hinge_age = 50)
sero_none <- rep(0, 9) # no prior immunity


### Read in vax strategies
#min_cases <- readRDS("optimal_BEL_cases_50x.RData")
min_cases <- readRDS("optimal_POL_cases_best.RData")
V0_oldest <- readRDS("V0_oldest_strat_POL.RData")


list_mincases   <- vector(mode = "list")
list_minmort   <- vector(mode = "list")

ptm <- proc.time()
for (i in seq(0, 50, by = 1)){
  j <- i/100
  V_0_temp <- as.numeric(min_cases[1+i, 2:10] * N_i/100)
  list_mincases[[paste0(i)]] <- run_sim_WHOstrat(C, percent_vax = j, V_0 = V_0_temp, u = u_var)
  
  V_0_temp <- V0_oldest[1+i,]
  list_minmort[[paste0(i)]] <- run_sim_WHOstrat(C, percent_vax = j, V_0 = V_0_temp, u = u_var)
}
proc.time() - ptm

# Plotting

theme_set(theme_minimal(base_size = 26))
total_cases <- rep(NA, 102)
total_deaths <- rep(NA, 102)
count <- 1
for (i in list_mincases){
  total_cases[count] <- compute_total_cases(i)
  print(total_cases[count])
  total_deaths[count] <- compute_total_deaths(i)
  count <- count + 1
}
for (i in list_minmort){
  total_cases[count] <- compute_total_cases(i)
  total_deaths[count] <- compute_total_deaths(i)
  count <- count + 1
}

baseline_cases <- compute_total_cases(list_mincases$`0`)
baseline_cases <- c(rep(baseline_cases, 102))

baseline_deaths <- compute_total_deaths(list_minmort$`0`)
baseline_deaths <- c(rep(baseline_deaths, 102))

reduction_in_cases <- (1-(total_cases/baseline_cases))*100
reduction_in_deaths <- (1-(total_deaths/baseline_deaths))*100
vax_avail <- c((seq(0, 50, by = 1)), (seq(0, 50, by = 1)))
strat <- c(rep("min cases", 51), rep("min mort", 51))
df <- data.frame(vax_avail, strat, reduction_in_cases, reduction_in_deaths)

ggplot(df, aes(x = vax_avail, y = reduction_in_cases, col = strat, fill = strat)) +
  geom_line(size = 2, alpha = 1) +
  xlab("Total vaccine supply (% of pop)") +
  ylab("Reduction in infections (%)") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5, type = "discrete")[4:5],
                     name = "Allocation Strategy", labels = c("Minimize infections", "Oldest first")) +
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + # , breaks = c()) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 50)) + #, breaks = c()) +
  theme(legend.position = "none",
  ) +
  guides(colour = guide_legend(override.aes = list(size=6)))
 

ggplot(df, aes(x = vax_avail, y = reduction_in_deaths, col = strat, fill = strat)) +
  geom_line(size = 2, alpha = 1) +
  xlab("Total vaccine supply (% of pop)") +
  ylab("Reduction in deaths (%)") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5, type = "discrete")[4:5],
                     name = "Allocation Strategy", labels = c("Minimize infections", "Oldest first")) +
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + # , breaks = c()) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 50)) + #, breaks = c()) +
  theme(legend.position = "none",
  ) + 
  guides(colour = guide_legend(override.aes = list(size=6)))

######################################################
######################################################
vax_avail <- c(seq(0, 50, by = 1))
diff_cases <- df[df$strat == "min cases", ]$reduction_in_cases - df[df$strat == "min mort", ]$reduction_in_cases
diff_deaths <- df[df$strat == "min cases", ]$reduction_in_deaths - df[df$strat == "min mort", ]$reduction_in_deaths
  
diff_df <- data.frame(vax_avail, diff_cases, diff_deaths)

#700 * 550
ggplot(diff_df, aes(x = vax_avail, y = diff_cases)) + 
  geom_line(size = 1.5) + 
  xlab("Total vaccine supply (% of pop)") +
  ylab("Difference in percent\nreduction in infections") + 
  scale_y_continuous(expand = c(0,0), limit = c(0, 60), breaks = c(0, 20, 40, 60)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 50)) #, breaks = c()) +

ggplot(diff_df, aes(x = vax_avail, y = diff_deaths)) + 
  geom_line(size = 1.5) + 
  xlab("Total vaccine supply (% of pop)") +
  ylab("Difference in percent\nreduction in deaths") + 
  scale_y_continuous(expand = c(0,0), limit = c(-80, 5))#, breaks = c(0, 20, 40, 60)) +
  scale_x_continuous(expand = c(0,0), limit = c(0, 50)) #, breaks = c()) +
