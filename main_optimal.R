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

constraints <- function(x){
  # function defining the inequality constraints s.t. constraints >= 0 for all components
  h <- numeric(11)
  for (i in 1:9){
    h[i] <- x[i] # vax distribution > 0
    #h[i] <- x[i] - previous_step[i] # greedy
  }
  for (i in 1:9){
    h[i + 9] <- 100 - (x[i]/N_i[i]) # max vax is 100%
  }
  h[10] <- sum(x) - nvax
  h[11] <- nvax - sum(x)
  
  return(h)
}

get_initial_vax <- function(nvax){
  # vaccinate oldest groups first (used when minimizing deaths)
  initial_vax <- rep(0, num_groups)
  vax_left <- nvax
  count <- 9
  
  while (vax_left > 0){
    if (vax_left < N_i[count]){
      initial_vax[count] <- vax_left
      vax_left <- 0
    } else {
      initial_vax[count] <- N_i[count]
      vax_left <- vax_left - N_i[count]
      count <- count - 1
    }
  }
  initial_vax
}

# _____________________________________________________________________
# SET UP ----
# _____________________________________________________________________
# Demographics
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

# susceptibility with R0 ~ 3 (R0 <- compute_R0(u, C))
u_constant     <- rep(0.022, 9) # constant # 0.02 for US, 0.022 for Belgium
u_var     <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/38.1 # R0 = 2.6 for BEL
R0 <- compute_R0(u_var, C)

# vaccine efficacy
v_e_constant <- get_v_e(p = 1, y0 = 1, hinge_age = 50)
v_e_var <- get_v_e(p = 0.5, y0 = 1, hinge_age = 50)

# serology 
sero_none <- rep(0, 9) # no prior immunity
sero_belgium <- c(0.03760, 0.08981, 0.07008, 0.05616,0.02732, 0.03709, 0.02071, 0.02646,0.03477) # Ref: Herzog

# _____________________________________________________________________
# RUN OPTIMAL SIM ----
# _____________________________________________________________________
# * RUN OPTIMIZATION CODE FOR ONE PERCENT_VAX ----
to_minimize <- "deaths"
percent_vax <- 0.15
nvax <- percent_vax*pop_total
initial_vax <- rep(nvax/num_groups, 9) # start with uniform distribution
#initial_vax <- get_initial_vax(nvax) # start with vaccinating oldest first

ptm <- proc.time()
optimal_vax <- cobyla(initial_vax, optimize_sim, hin = constraints,
            nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
proc.time() - ptm

vax <- round(optimal_vax$par/N_i * 100, 2)
groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
test <- data.frame(groups, vax)

theme_set(theme_minimal(base_size = 20))

ggplot(test, aes(y=vax, x=groups)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Age group") + 
  ylab("Vaccinated\n(% of age group)")+
  ylim(0, 101) +
  scale_x_discrete(breaks=c("0-9","20-29", "40-49", "60-69", "80+"))
  #ggtitle("Optimizing for least deaths")

# * RUN OPTIMIZATION CODE OVER MULTIPLE PERCENT_VAX ----
to_minimize <- "deaths"

ptm <- proc.time()
minimize_cases_df <- data.frame(matrix(ncol = 11, nrow = 51))
colnames(minimize_cases_df) <- c("Percent_vax", "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+", "Percent_infected")

list_optimal <- vector(mode = "list")

for (j in 1:8){
  for (i in 0:42){
    percent_vax <- i/100
    nvax <- percent_vax*pop_total
    # uniform with noise
    noise <- runif(num_groups, -0.04*nvax, 0.04*nvax)
    initial_vax <- rep(nvax/num_groups, 9) + noise
    
    #initial_vax <- get_initial_vax(nvax)

    optimal_vax <- cobyla(initial_vax, optimize_sim, hin = constraints,
                          nl.info = FALSE, control = list(xtol_rel = 1e-8, maxeval = 5000))
    minimize_cases_df[i+1, 2:10] <- round(optimal_vax$par/N_i * 100, 2)
    minimize_cases_df[i+1, 1] <- i
    minimize_cases_df[i+1, 11] <- optimal_vax$value
  }
  print(j)
  list_optimal[[j]] <- minimize_cases_df
}
proc.time() - ptm

#saveRDS(list_optimal, "optimal_BEL_deaths_newIFR_8x_6.RData")

# _____________________________________________________________________
# ANALYZE OPTIMAL ACROSS VAX AVAIALBLE ####
# _____________________________________________________________________
simplemodel_cases <- readRDS("optimal_BEL_cases_50x.RData")

simplemodel_deaths <- readRDS("optimal_BEL_deaths_50x.RData")
simplemodel_deaths <- readRDS("optimal_BEL_deaths_newIFR_x33.RData")
# calculate percent reduction in cases
total_cases <- c(simplemodel_cases$Percent_infected)
baseline_cases_C <- simplemodel_cases$Percent_infected[1]

num_per_list <- 51
baseline_cases <- c(rep(baseline_cases_C, num_per_list))
reduction_in_cases <- (1-(total_cases/baseline_cases))*100

# calculate percent reduction in deaths
total_deaths <- c(simplemodel_deaths$Percent_infected)
baseline_deaths_C <- simplemodel_deaths$Percent_infected[1]

baseline_deaths <- c(rep(baseline_deaths_C, num_per_list))
reduction_in_deaths <- (1-(total_deaths/baseline_deaths))*100

vax_avail <- c(rep(seq(0, 50, by = 1), 1))
strat <- c(rep("optimal", num_per_list*1))
variable <- c(rep("var", num_per_list))#, rep("var", num_per_list))

# make dataframe and plot
optimal_df_C <- data.frame(vax_avail, strat, reduction_in_cases, variable)
optimal_df_C <- data.frame(vax_avail, strat, reduction_in_deaths, variable)

ggplot(optimal_df_C, aes(x = vax_avail, y = reduction_in_deaths, col = strat, fill = strat)) +
  geom_abline(slope = 1, intercept = 0, size = 3, alpha = 0.2) +
  geom_line(aes(linetype = variable), size = 2, alpha = 0.9) +
  xlab("Vaccine available (% of total pop)") +
  ylab("% Reduction in cases") +
  #ggtitle("Optimizing for least cases") +
  scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                     labels =  c("Optimal")) +
  scale_fill_brewer(palette = "Dark2", name = "Allocation Strategy",
                    labels =  c("Optimal")) +
  scale_linetype_discrete(name = "contact matrix")+
                          #labels = c("Constant", "Age-dependent"))+
  ylim(0, 100) +
  xlim(0,50) +
  theme(legend.position = "none")

### Bar chart for paper figure 4
# green: "#00BA38"
# blue: "#619CFF"
# gold: "#E6AB02"
ggplot(new_df[new_df$Percent_vax == 30,], aes(y=num, x=age_group, fill = as.factor(Percent_vax), group = as.factor(Percent_vax))) + 
  geom_bar(position = "stack", stat = "identity", fill = "#E6AB02") +
  #xlab("Age group") + 
  ylab("Vaccinated (%)")+
  theme_classic(base_size = 22) +
  scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                   labels=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")) +
  #scale_x_discrete(breaks = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =element_blank())
        #axis.title.y = element_blank())
#scale_fill_manual(pal)

dist_of_vax <- matrix(nrow = 51, ncol = 9)

for (i in 1:51){
  tot <- sum(N_i*cases_best[i, 2:10])
  dist_of_vax[i,] <- as.numeric(((N_i*cases_best[i, 2:10])/tot)*100)
}

dist_of_vax <- as.data.frame(dist_of_vax)
colnames(dist_of_vax) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

dist_of_vax$Percent_vax <- 0:50

new_df <- dist_of_vax %>%
  select(Percent_vax, "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+") %>%
  gather(key = "age_group", value = "num", -Percent_vax)

p1 <- ggplot(new_df[new_df$Percent_vax == 10,], aes(y=num, x=age_group, fill = as.factor(Percent_vax), group = as.factor(Percent_vax))) + 
  geom_bar(position = "stack", stat = "identity", fill = "#E6AB02") +
  #xlab("Age group") + 
  #ylab("Age distribution of vaccines (%)")+
  theme_classic(base_size = 26) +
  scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                   labels=c("", "", "", "", "", "", "", "", "")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 75), breaks = c(0,25,50,75)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =element_blank(),
        axis.title.y =element_blank())

p2 <- ggplot(new_df[new_df$Percent_vax == 20,], aes(y=num, x=age_group, fill = as.factor(Percent_vax), group = as.factor(Percent_vax))) + 
  geom_bar(position = "stack", stat = "identity", fill = "#E6AB02") +
  #xlab("Age group") + 
  #ylab("Age distribution of vaccines (%)")+
  theme_classic(base_size = 26) +
  scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                   labels=c("", "", "", "", "", "", "", "", "")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 75), breaks = c(0,25,50,75)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =element_blank(),
        axis.title.y =element_blank())

p3 <- ggplot(new_df[new_df$Percent_vax == 30,], aes(y=num, x=age_group, fill = as.factor(Percent_vax), group = as.factor(Percent_vax))) + 
  geom_bar(position = "stack", stat = "identity", fill = "#E6AB02") +
  #xlab("Age group") + 
  #ylab("Age distribution of vaccines (%)")+
  theme_classic(base_size = 26) +
  scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                   labels=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")) +
  #scale_x_discrete(breaks = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 75), breaks = c(0,25,50,75)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =element_blank(), 
        axis.title.y =element_blank())

# export as 600 * 1200
grid.arrange(p1, p2, p3, 
             ncol=1, widths=c(2.3),
             heights = c(2.3,2.3,2.8),
             left = textGrob("Age distribution of vaccines (%)", rot = 90, vjust = 0.5, hjust = 0.35, gp = gpar(fontsize = 26)),
             bottom = textGrob("Age (years)", vjust = 0, gp = gpar(fontsize = 26)))
