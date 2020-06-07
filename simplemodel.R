# Vaccine Strategy Simple Model

setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________

library(tidyverse) 
library(deSolve) # ode solver
library(xlsx)
library(gridExtra)

# _____________________________________________________________________
# FUNCTIONS ----
# _____________________________________________________________________

calculate_derivatives=function(t, x, parameters){
  # the parameters in the parameters list are:
  #    the probability of transmission on contact, beta
  #    the incubation period, nu
  #    the recovery period, gamma
  #    the contact matrix, C, that is the # contacts per day among age groups
  #
  # Note that x is a vector of length (#model compartment types)*(#age classes)
  # Thus, S, E, I and R are vectors, all of length nage
  ncompartment <- 4
  nage <- length(x)/ncompartment
  S    <- as.matrix(x[1:nage])
  E    <- as.matrix(x[(nage+1):(2*nage)])
  I    <- as.matrix(x[(2*nage+1):(3*nage)])
  R    <- as.matrix(x[(3*nage+1):(4*nage)])
  
  I[I<0] = 0
  with(as.list(parameters),{
    # note that because S, I and R are all vectors of length nage, so will N,
    # and dS, dI, and dR
    N = S+E+I+R
    dS = -as.matrix(S*beta)*(as.matrix(C)%*%as.matrix(I/N))
    dE = -dS - nu*as.matrix(E)
    dI = nu*as.matrix(E) - gamma*as.matrix(I)
    dR = +gamma*as.matrix(I)
    
    out=c(dS,dE,dI,dR)
    list(out)
  })
}

plot_all_ages_overtime=function(df, compartment){
  # INPUTS:
  # df: simulation to plot
  # compartment: character of compartment to plot i.e. "I", "R"
  #
  # OUTPUT:
  # p: ggplot object of plot
  
  # make new df for plotting
  col_to_gather <- vector(mode="character", length=nage)
  for (i in 1:9) {col_to_gather[i] <- paste0(compartment,i)}
  
  new_df <- df %>%
    select(time, col_to_gather) %>%
    gather(key = "age_group", value = "num", -time)
  
  # Add col for percent infected & recovered
  new_df$percent <- new_df$num
  for (i in 1:length(new_df$num)){
    if (new_df$age_group[i] == col_to_gather[1]) {new_df$percent[i] <- new_df$percent[i]/N[1]}
    else if (new_df$age_group[i] == col_to_gather[2]) {new_df$percent[i] <- new_df$percent[i]/N[2]}
    else if (new_df$age_group[i] == col_to_gather[3]) {new_df$percent[i] <- new_df$percent[i]/N[3]}
    else if (new_df$age_group[i] == col_to_gather[4]) {new_df$percent[i] <- new_df$percent[i]/N[4]}
    else if (new_df$age_group[i] == col_to_gather[5]) {new_df$percent[i] <- new_df$percent[i]/N[5]}
    else if (new_df$age_group[i] == col_to_gather[6]) {new_df$percent[i] <- new_df$percent[i]/N[6]}
    else if (new_df$age_group[i] == col_to_gather[7]) {new_df$percent[i] <- new_df$percent[i]/N[7]}
    else if (new_df$age_group[i] == col_to_gather[8]) {new_df$percent[i] <- new_df$percent[i]/N[8]}
    else if (new_df$age_group[i] == col_to_gather[9]) {new_df$percent[i] <- new_df$percent[i]/N[9]}
  }
  
  # Plot
  theme_set(theme_minimal(base_size = 15))
  
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    geom_line(aes(color = age_group), size = 2) +
    xlab("Time") +
    scale_color_brewer(palette = "Spectral", name = "Age Group", 
                       labels =  c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
  
  if (compartment == "I") {
    p <- p + ylab("Percent Infected") + ylim(0,0.4)
  } else if (compartment == "R") {
    p <- p + ylab("Percent Recovered") + ylim(0,1)
  }
}

plot_one_age_overtime = function(col_name, age_group_num) {
  # INPUTS:
  # col_name: "I2" == infected, age group 2
  # age_group_num: age group of interest, 1 = 0-9, 2 = 10-18, ...
  #
  # OUTPUT:
  # p: ggplot object plotting all strats for one age group & one compartment
  
  nstrat <- 3
  new_df <- data.frame(matrix(ncol = nstrat+1, nrow = dim(df_novax)[1]))
  colnames(new_df) <- c("time", "no_vax", "everyone", "kids")
  new_df$time <- df_novax$time
  
  # Make sure col of interest are stored as vectors
  new_df$no_vax <- unlist(as.data.frame(df_novax)[col_name]/N[age_group_num])
  new_df$everyone <- unlist(as.data.frame(df_prop)[col_name]/N[age_group_num])
  new_df$kids <- unlist(as.data.frame(df_propkids)[col_name]/N[age_group_num])
  
  new_df <- new_df %>%
    select(time, no_vax, everyone, kids) %>%
    gather(key = "strat", value = "percent", -time)
  
  theme_set(theme_minimal(base_size = 15))
  
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    geom_line(aes(color = strat), size = 2) +
    xlab("Time") +
    scale_color_brewer(palette = "Dark2", name = "Scenario", 
                       labels =  c("Everyone", "Just Kids", "No vaccines"))
  
  if (substr(col_name, 1, 1) == "I"){
    p <- p + ylab("Percent Infected") + ylim(0, 0.4)
  }
}  
  
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
# _____________________________________________________________________
# SETUP ----
# _____________________________________________________________________
# Demographics
C <- readRDS("C_USA_bytens_all.RData")

frac_age <- read.xlsx("USA_demographics.xlsx", "data_bytens", header = FALSE, row.names = TRUE)
frac_age <- as.vector(frac_age[[1]])

#barplot(frac_age, xlab = ("Age Groups"), ylab = "Percent", 
#        names.arg = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))

npop <- 1000000
N <-  npop*frac_age/100      
nage <- length(frac_age)

# Initialize simulation with 1 infected in each age group 
I_0    <- rep(1,nage)
S_0    <- N-I_0
E_0    <- rep(0,nage)
R_0    <- rep(0,nage)

inits_novax <- c(S=S_0,E=E_0,I=I_0,R=R_0)

# Disease Tranmission
nu    <- 1/3 # Davies
gamma <- 1/5        
R0    <- 1.5       
beta  <- 0.05 # TODO: double check how Towers calculates beta using R0

# _____________________________________________________________________
# VACCINE STRATEGIES ----
#     prop: distribute proportionally to each age group
#     propkids: distribute proportionally to age groups < 20
#     propadults: distribute proportionally to age groups 20-49
# _____________________________________________________________________
nvax <- 0.10*npop 

# prop strategy
vax_prop <- nvax*frac_age/100

I_0    <- rep(1,nage)
S_0    <- N-I_0-vax_prop
E_0    <- rep(0,nage)
R_0    <- rep(0,nage)

inits_prop <- c(S=S_0,E=E_0,I=I_0,R=R_0)

# propkids strategy
nkids <- N[1] + N[2]
vax_dist_propkids <- rep(0, nage)
vax_dist_propkids[1] <- N[1]/nkids
vax_dist_propkids[2] <- N[2]/nkids

vax_propkids <- nvax*vax_dist_propkids

I_0    <- rep(1,nage)
S_0    <- N-I_0-vax_propkids
E_0    <- rep(0,nage)
R_0    <- rep(0,nage)

inits_propkids <- c(S=S_0,E=E_0,I=I_0,R=R_0)

# propadults strategy
nadults <- N[3] + N[4] + N[5]
vax_dist_propadults <- rep(0, nage)
vax_dist_propadults[3] <- N[3]/nadults
vax_dist_propadults[4] <- N[4]/nadults
vax_dist_propadults[5] <- N[5]/nadults

vax_propadults <- nvax*vax_dist_propadults

I_0    <- rep(1,nage)
S_0    <- N-I_0-vax_propadults
E_0    <- rep(0,nage)
R_0    <- rep(0,nage) 

inits_propadults <- c(S=S_0,E=E_0,I=I_0,R=R_0)
# _____________________________________________________________________
# NUMERICALLY SOLVE ----
# _____________________________________________________________________
parameters <- c(beta=beta, nu=nu, gamma=gamma, C=C)

t <- seq(0,80,1) 

df_novax <- as.data.frame(lsoda(inits_novax, t, calculate_derivatives, parameters)) %>% as_tibble()
df_prop <- as.data.frame(lsoda(inits_prop, t, calculate_derivatives, parameters)) %>% as_tibble()
df_propkids <- as.data.frame(lsoda(inits_propkids, t, calculate_derivatives, parameters)) %>% as_tibble()
df_propadults <- as.data.frame(lsoda(inits_propadults, t, calculate_derivatives, parameters)) %>% as_tibble()

# _____________________________________________________________________
# RESULTS ----
# _____________________________________________________________________

# cat("The fraction of kids that were infected is ",max(mymodel_results$R1)/N[1],"\n")
# cat("The fraction of adults that were infected is ",max(mymodel_results$R2)/N[2],"\n")
# cat("The total final size is ",max(mymodel_results$R1+mymodel_results$R2)/npop,"\n")

# * Plot for all ages over time: infected & recovered ----
compartment <- "I"

p_novaxI <- plot_all_ages_overtime(df_novax, compartment) +
  theme(legend.position = "none") +
  ggtitle("No Vaccines")

p_propI <- plot_all_ages_overtime(df_prop, compartment) +
  theme(legend.position = "none") +
  ggtitle("Proportional to all ages")

p_propkidsI <- plot_all_ages_overtime(df_propkids, compartment) +
  ggtitle("Proportional to just kids (<20)")

legend <- get_legend(p_propkidsI)
p_propkidsI <- p_propkidsI + theme(legend.position = "none")

p_propadultsI <- plot_all_ages_overtime(df_propadults, compartment) +
  theme(legend.position = "none") +
  ggtitle("Proportional to adults (20-49)")

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
compartment <- "R"
p_novaxR <- plot_all_ages_overtime(df_novax, compartment) +
  theme(legend.position = "none") +
  ggtitle("No Vaccines")

p_propR <- plot_all_ages_overtime(df_prop, compartment) +
  theme(legend.position = "none") +
  ggtitle("Proportional to all ages")

p_propkidsR <- plot_all_ages_overtime(df_propkids, compartment) +
  theme(legend.position = "none") +
  ggtitle("Proportional to just kids (<20)")

grid.arrange(p_novaxI, p_propI, p_propkidsI, legend, p_novaxR, p_propR, p_propkidsR, ncol=4, widths=c(2.3, 2.3, 2.3, 0.8))
grid.arrange(p_novaxI, p_propI, legend, p_propkidsI,  p_propadultsI, ncol=3, widths=c(2.3, 2.3, 0.8))


# * Plot for one age group, different strats ----

p_kids_I <- plot_one_age_overtime(col_name = "I2", age_group_num = 2) +
  theme(legend.position = "none") +
  ggtitle("Ages 10-19")

p_young_I <- plot_one_age_overtime(col_name = "I4", age_group_num = 4) +
  theme(legend.position = "none") +
  ggtitle("Ages 30-39")

p_middle_I <- plot_one_age_overtime(col_name = "I6", age_group_num = 6) +
  theme(legend.position = "none") +
  ggtitle("Ages 50-59")

p_old_I <- plot_one_age_overtime(col_name = "I8", age_group_num = 8) +
  ggtitle("Ages 70+")

legend <- get_legend(p_old_I)
p_old_I <- p_old_I + theme(legend.position = "none")

grid.arrange(p_kids_I, p_young_I, p_middle_I, p_old_I, legend, ncol=5, widths=c(2.3, 2.3, 2.3, 2.3, 0.8))
